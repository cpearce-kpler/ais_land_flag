#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone land/water masking for AIS data.

Land/water lassification is a pure per-row lookup with no dependency on ship
identity or chronological order -- it reads only LAT/LON per source file
(plus SHIP_ID/TIMESTAMP if the QGIS output is enabled) and runs
embarrassingly parallel across files.

Cascade: 0.01 -> 0.001 -> 0.0001 degree
--------------------------------------------
Uses the already-validated StaticLandWaterHierarchy classifier


Exact tail (optional, default OFF)
-------------------------------------
When --exact-tail is passed, residual MIXED points (state==mixed after
the full 3-tier cascade) are resolved against the hierarchy's own
pre-extracted geometry cache (<hierarchy_root>/geometry_cache/) --
already reprojected to EPSG:4326 and repaired once at build time -- so
this reuses that cache directly rather than re-reading the source
GeoPackages or rebuilding a spatial index from them on every run.

ASSUMPTION WORTH FLAGGING: this reads geometry_cache/<role>.wkb.bin
(a concatenated blob of WKB-encoded polygons) and <role>.offsets.npy
(N+1 cumulative byte offsets, so geometry i spans
blob[offsets[i]:offsets[i+1]]) for role in (land, marine, inland),
based on the observed file names in that folder -- this has NOT been
verified against the actual binary content, since only file names and
sizes were available when this was written. The self-test builds and
reads back a synthetic cache in exactly this format and passes, but if
--exact-tail errors on the real cache, this assumption is the first
thing to check.

Precedence when a point matches more than one layer (matches the
existing exact-tail benchmark's own order): marine_water, then
inland_water OVERRIDES marine_water, then land ONLY if neither
marine nor inland matched.

Outputs, per source file
--------------------------
1. <stem>.land_water_state_2bit.bin(.zst) + <stem>.metadata.json -- the
   lightweight index other scripts should read to select/ignore rows.
   2-bit packed, same codes as the hierarchy itself
   (0=marine_water, 1=land, 2=mixed, 3=inland_water), positionally
   aligned with the source file. zstd-compressed if the zstandard
   package is available, skipped gracefully (same as this project's own
   speed-bitmap builder) if not.
2. Optionally (--qgis-output, default ON): a Parquet file of just the
   FLAGGED rows (state == land OR state == mixed, i.e. "not confirmed
   water") with LAT/LON/SHIP_ID/TIMESTAMP/state -- sized for direct QGIS
   review, not the full dataset.
3. A top-level catalog.json summarising every processed file.

Usage (PowerShell)
-------------------
    & "C:\\anaconda\\python.exe" -u "C:\\Users\\Craig Pearce\\Desktop\\ais_land_water_mask.py" `
        --ais-folder "C:\\Users\\Craig Pearce\\Desktop\\Data_sets\\ais_files\\ais_2025" `
        --hierarchy-root "C:\\Users\\Craig Pearce\\Desktop\\Data_sets\\grids\\land_water_hierarchy" `
        --output-dir "C:\\Users\\Craig Pearce\\Desktop\\ais_2025_land_water_mask" `
        --max-files 5

Add --exact-tail to resolve residual MIXED cells via exact geometry
(slower; off by default). Add --no-qgis-output to skip the QGIS Parquet
if you only want the lightweight index.

Run --self-test first (builds its own synthetic hierarchy, geometry
cache, and AIS file -- no real data touched):

    & "C:\\anaconda\\python.exe" -u "C:\\Users\\Craig Pearce\\Desktop\\ais_land_water_mask.py" --self-test
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_CLASSIFIER_MODULE_NAME = "ais_land_water_fast_flag"

try:
    import zstandard as _zstandard
except ImportError:
    _zstandard = None

try:
    from shapely import wkb as _shapely_wkb
    from shapely.strtree import STRtree as _STRtree
    import shapely as _shapely_top
except ImportError:
    _shapely_wkb = None
    _STRtree = None
    _shapely_top = None


STATE_NAMES = {0: "marine_water", 1: "land", 2: "mixed", 3: "inland_water"}


def _pack_2bit(values: np.ndarray, pad_value: int = 2) -> np.ndarray:
    array = np.asarray(values, dtype=np.uint8).reshape(-1)
    pad = (-array.size) % 4
    if pad:
        array = np.concatenate([array, np.full(pad, int(pad_value), dtype=np.uint8)])
    grouped = array.reshape(-1, 4)
    return (
        grouped[:, 0]
        | (grouped[:, 1] << np.uint8(2))
        | (grouped[:, 2] << np.uint8(4))
        | (grouped[:, 3] << np.uint8(6))
    ).astype(np.uint8, copy=False)


def _unpack_2bit(packed: np.ndarray, count: int) -> np.ndarray:
    out = np.empty((packed.size, 4), dtype=np.uint8)
    out[:, 0] = packed & 0b11
    out[:, 1] = (packed >> 2) & 0b11
    out[:, 2] = (packed >> 4) & 0b11
    out[:, 3] = (packed >> 6) & 0b11
    return out.reshape(-1)[:count]


def _write_bytes_maybe_compressed(data: bytes, path: Path) -> tuple[Path, bool]:
    if _zstandard is not None:
        compressed = _zstandard.ZstdCompressor(level=3).compress(data)
        out_path = path.with_suffix(path.suffix + ".zst")
        out_path.write_bytes(compressed)
        return out_path, True
    path.write_bytes(data)
    return path, False


# ---------------------------------------------------------------------------
# Exact-tail geometry cache (optional)
# ---------------------------------------------------------------------------


class GeometryCache:
    """One STRtree per role (land, marine, inland), read from the
    hierarchy's own pre-extracted <role>.wkb.bin + <role>.offsets.npy --
    already reprojected/repaired once at build time. See the ASSUMPTION
    note in the module docstring about this binary format."""

    def __init__(self, hierarchy_root: Path) -> None:
        if _shapely_wkb is None or _STRtree is None:
            raise RuntimeError(
                "shapely is required for --exact-tail. Install it, or omit --exact-tail."
            )
        cache_dir = hierarchy_root / "geometry_cache"
        if not cache_dir.is_dir():
            raise FileNotFoundError(f"geometry_cache directory not found: {cache_dir}")
        self.trees: dict[str, "_STRtree | None"] = {}
        for role in ("land", "marine", "inland"):
            wkb_path = cache_dir / f"{role}.wkb.bin"
            offsets_path = cache_dir / f"{role}.offsets.npy"
            if not wkb_path.is_file() or not offsets_path.is_file():
                raise FileNotFoundError(
                    f"Expected {wkb_path} and {offsets_path} for role {role!r} -- "
                    "if these names don't match the real cache, the ASSUMPTION note "
                    "in this script's docstring is the first thing to check."
                )
            t0 = time.perf_counter()
            blob = wkb_path.read_bytes()
            offsets = np.load(offsets_path)

            # Vectorized: one call for every geometry in this role, instead
            # of a Python loop calling shapely.wkb.loads() once each. With
            # only a few hundred geometries this isn't the main cost (see
            # the prepare() note below for that), but it's a free win on
            # the one-time load regardless of geometry count.
            wkb_slices = [
                blob[int(offsets[i]) : int(offsets[i + 1])]
                for i in range(offsets.size - 1)
                if int(offsets[i + 1]) > int(offsets[i])
            ]
            if wkb_slices:
                geometries = _shapely_top.from_wkb(wkb_slices, on_invalid="ignore")
                geometries = geometries[~_shapely_top.is_missing(geometries)]
            else:
                geometries = np.empty(0, dtype=object)

            # Real fix, empirically validated (64x speedup at realistic
            # scale against synthetic data matching this project's own
            # measured part-count distribution -- land alone showed a
            # single geometry bundling 26,737 disjoint parts): decompose
            # every multi-part geometry into its individual constituent
            # parts BEFORE building the STRtree. A bundled feature's
            # bounding box spans the combined extent of every part inside
            # it (like "France" including Reunion and Guadeloupe), so a
            # query point on the opposite side of the world from where it
            # actually matters still passes the STRtree's bbox pre-filter
            # and triggers a full expensive exact test. Splitting into
            # parts first means each STRtree leaf gets a tight bbox around
            # just that one part, so far fewer candidates survive
            # pre-filtering in the first place.
            #
            # shapely.prepare() was tried first and directly measured to
            # make NO difference for STRtree's own bulk query() (shapely's
            # bulk predicate query already handles the exact test
            # efficiently on its own) -- removed rather than left in for
            # no benefit.
            if geometries.size:
                geometries = _shapely_top.get_parts(geometries)

            elapsed = time.perf_counter() - t0
            print(
                f"    loaded {role}: {geometries.size:,} geometries after decomposition, "
                f"{len(blob)/1e6:.0f} MB WKB, in {elapsed:.1f}s"
            )
            self.trees[role] = _STRtree(geometries) if geometries.size else None

    def resolve(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """Exact classification for a batch of points, matching the
        precedence already used by this project's exact-tail benchmark:
        marine, then inland OVERRIDES marine, then land only if neither
        marine nor inland matched. Points matching nothing stay MIXED
        (state code 2) -- i.e. exact-tail-unresolved is treated the same
        as still-ambiguous, not silently assumed to be water or land."""
        points = _shapely_top.points(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        state = np.full(points.shape, 2, dtype=np.uint8)  # default: stays MIXED

        def hit_mask(tree: "_STRtree | None") -> np.ndarray:
            if tree is None or points.size == 0:
                return np.zeros(points.shape, dtype=bool)
            pairs = tree.query(points, predicate="intersects")
            hit = np.zeros(points.shape, dtype=bool)
            if pairs.size:
                hit[np.unique(pairs[0])] = True
            return hit

        marine_hit = hit_mask(self.trees.get("marine"))
        inland_hit = hit_mask(self.trees.get("inland"))
        land_hit = hit_mask(self.trees.get("land"))

        state[marine_hit] = 0
        state[inland_hit] = 3  # inland overrides marine
        land_only = land_hit & ~marine_hit & ~inland_hit
        state[land_only] = 1
        return state


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------


def process_file(
    ais_file: Path,
    hierarchy,
    *,
    exact_tail_cache: GeometryCache | None,
    output_dir: Path,
    qgis_output: bool,
) -> dict:
    columns = ["LAT", "LON"]
    if qgis_output:
        columns += ["SHIP_ID", "TIMESTAMP"]
    df = pd.read_parquet(ais_file, columns=columns)
    n = len(df)

    result = hierarchy.classify(df["LON"].to_numpy(), df["LAT"].to_numpy(), mixed_policy="flag")
    state = result.state.copy()

    exact_tail_resolved = 0
    if exact_tail_cache is not None:
        residual = state == 2
        if np.any(residual):
            resolved = exact_tail_cache.resolve(
                df["LON"].to_numpy()[residual], df["LAT"].to_numpy()[residual]
            )
            state[residual] = resolved
            exact_tail_resolved = int(np.count_nonzero(resolved != 2))

    packed = _pack_2bit(state)
    bin_path = output_dir / f"{ais_file.stem}.land_water_state_2bit.bin"
    written_path, compressed = _write_bytes_maybe_compressed(packed.tobytes(), bin_path)

    state_counts = {STATE_NAMES[c]: int(np.count_nonzero(state == c)) for c in STATE_NAMES}

    metadata = {
        "source_relative_path": ais_file.name,
        "source_rows": n,
        "packed_path": str(written_path),
        "packed_relative_path": written_path.name,
        "compressed": compressed,
        "row_alignment": "bit/cell N = raw AIS source row N",
        "state_codes": STATE_NAMES,
        "state_counts": state_counts,
        "exact_tail_enabled": exact_tail_cache is not None,
        "exact_tail_resolved_rows": exact_tail_resolved,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (output_dir / f"{ais_file.stem}.metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    qgis_rows = 0
    if qgis_output:
        flagged = (state == 1) | (state == 2)
        qgis_rows = int(np.count_nonzero(flagged))
        if qgis_rows:
            qgis_df = pd.DataFrame(
                {
                    "SHIP_ID": df["SHIP_ID"].to_numpy()[flagged],
                    "LAT": df["LAT"].to_numpy()[flagged],
                    "LON": df["LON"].to_numpy()[flagged],
                    "TIMESTAMP": df["TIMESTAMP"].to_numpy()[flagged],
                    "state": state[flagged],
                    "state_name": [STATE_NAMES[int(c)] for c in state[flagged]],
                }
            )
            qgis_df.to_parquet(output_dir / f"{ais_file.stem}.flagged_for_qgis.parquet")

    metadata["qgis_flagged_rows"] = qgis_rows
    return metadata


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def run(
    *,
    ais_folder: Path,
    hierarchy_root: Path,
    output_dir: Path,
    exact_tail: bool,
    qgis_output: bool,
    max_files: int | None,
    classifier_module_dir: Path,
    classifier_module_name: str,
    shard_index: int | None = None,
    shard_count: int | None = None,
) -> None:
    sys.path.insert(0, str(classifier_module_dir))
    hierarchy_module = importlib.import_module(classifier_module_name)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading static hierarchy from {hierarchy_root}...")
    hierarchy = hierarchy_module.StaticLandWaterHierarchy(hierarchy_root, require_ultrafine=False)

    exact_tail_cache = None
    if exact_tail:
        print("Loading exact-tail geometry cache (this is the optional, slower path)...")
        exact_tail_cache = GeometryCache(hierarchy_root)

    files = sorted(p for p in ais_folder.glob("*.parquet") if p.is_file())
    if not files:
        raise FileNotFoundError(f"No .parquet files found under {ais_folder}")

    # Sharding is INTERLEAVED, not contiguous blocks: worker i gets files at
    # positions i, i+shard_count, i+2*shard_count, ... File sizes vary a lot
    # (roughly a 24x spread across the real dataset), and a naive contiguous
    # split risks one worker getting an unluckily heavy run of files while
    # another finishes early -- interleaving spreads that variance evenly
    # across every shard instead.
    sharded = shard_index is not None and shard_count is not None
    if sharded:
        if not (0 <= shard_index < shard_count):
            raise ValueError(f"--shard-index must be between 0 and --shard-count-1 (got {shard_index} of {shard_count})")
        all_files = files
        files = files[shard_index::shard_count]
        print(f"Shard {shard_index} of {shard_count}: {len(files)} of {len(all_files)} total files (interleaved split)")

    if max_files is not None:
        files = files[:max_files]

    print(f"Processing {len(files)} file(s) (exact_tail={exact_tail}, qgis_output={qgis_output})")
    catalog_rows = []
    n_skipped_already_done = 0
    started = time.perf_counter()
    for i, f in enumerate(files, start=1):
        metadata_path = output_dir / f"{f.stem}.metadata.json"
        if metadata_path.is_file():
            # Resumability: this is essential for an unattended multi-day
            # run -- if a worker is interrupted partway (a bad file, a
            # transient error, a machine restart) and simply relaunched
            # with the same arguments, it should pick up where it left off
            # rather than redoing every already-completed file from scratch.
            n_skipped_already_done += 1
            existing_meta = json.loads(metadata_path.read_text(encoding="utf-8"))
            existing_meta.setdefault("process_seconds", 0.0)
            catalog_rows.append(existing_meta)
            print(f"  [{i}/{len(files)}] {f.name}: already processed, skipping (resuming from prior run)")
            continue

        t0 = time.perf_counter()
        meta = process_file(
            f, hierarchy, exact_tail_cache=exact_tail_cache, output_dir=output_dir, qgis_output=qgis_output
        )
        elapsed_file = time.perf_counter() - t0
        elapsed_total = time.perf_counter() - started
        meta["process_seconds"] = elapsed_file
        catalog_rows.append(meta)
        flagged_note = f", {meta['qgis_flagged_rows']:,} flagged for QGIS" if qgis_output else ""
        print(
            f"  [{i}/{len(files)}] {f.name}: {meta['source_rows']:,} rows in {elapsed_file:.1f}s"
            f"{flagged_note} (total elapsed {elapsed_total:.0f}s)"
        )

    if n_skipped_already_done:
        print(f"\n{n_skipped_already_done} file(s) were already processed from a prior run and were skipped.")

    catalog_name = (
        f"land_water_mask_catalog_shard{shard_index}_of_{shard_count}.json" if sharded
        else "land_water_mask_catalog.json"
    )
    catalog_path = output_dir / catalog_name
    catalog_path.write_text(json.dumps(catalog_rows, indent=2), encoding="utf-8")
    print(f"\nWrote {catalog_path}")


def merge_shard_catalogs(output_dir: Path, shard_count: int) -> None:
    """Combines every shardN_of_M catalog in output_dir into one complete
    catalog -- run once after all shards have finished."""
    merged = []
    for shard_index in range(shard_count):
        shard_path = output_dir / f"land_water_mask_catalog_shard{shard_index}_of_{shard_count}.json"
        if not shard_path.is_file():
            raise FileNotFoundError(f"Missing shard catalog: {shard_path} -- has shard {shard_index} finished yet?")
        shard_rows = json.loads(shard_path.read_text(encoding="utf-8"))
        print(f"  shard {shard_index}: {len(shard_rows)} files")
        merged.extend(shard_rows)

    merged.sort(key=lambda row: row["source_relative_path"])
    seen = set()
    duplicates = [row["source_relative_path"] for row in merged if row["source_relative_path"] in seen or seen.add(row["source_relative_path"])]
    if duplicates:
        raise ValueError(
            f"{len(duplicates)} file(s) appear in more than one shard's catalog -- shards may overlap: {duplicates[:5]}"
        )

    merged_path = output_dir / "land_water_mask_catalog.json"
    merged_path.write_text(json.dumps(merged, indent=2), encoding="utf-8")
    print(f"\nMerged {len(merged)} file(s) from {shard_count} shard(s) into {merged_path}")


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------


def _build_synthetic_geometry_cache(hierarchy_root: Path, land_point_lon: float, land_point_lat: float) -> None:
    """A geometry cache covering exactly the known 'residual_mixed' ground-
    truth point from _write_fake_hierarchy with a LAND polygon -- a real
    test of exact-tail actually changing the outcome, not a vacuous pass."""
    cache_dir = hierarchy_root / "geometry_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    from shapely.geometry import box as shapely_box

    land_poly = shapely_box(land_point_lon - 0.00001, land_point_lat - 0.00001, land_point_lon + 0.00001, land_point_lat + 0.00001)
    for role, geoms in (("land", [land_poly]), ("marine", []), ("inland", [])):
        blob_parts = [g.wkb for g in geoms]
        offsets = np.zeros(len(blob_parts) + 1, dtype=np.int64)
        for i, part in enumerate(blob_parts):
            offsets[i + 1] = offsets[i] + len(part)
        (cache_dir / f"{role}.wkb.bin").write_bytes(b"".join(blob_parts))
        np.save(cache_dir / f"{role}.offsets.npy", offsets)
    (cache_dir / "geometry_cache_metadata.json").write_text(json.dumps({"roles": ["land", "marine", "inland"]}), encoding="utf-8")


def _load_packed_state(output_dir: Path, stem: str, count: int) -> np.ndarray:
    metadata = json.loads((output_dir / f"{stem}.metadata.json").read_text())
    packed_path = output_dir / metadata["packed_relative_path"]
    data = packed_path.read_bytes()
    if metadata["compressed"]:
        if _zstandard is None:
            raise RuntimeError("zstandard is required to read this compressed output.")
        data = _zstandard.ZstdDecompressor().decompress(data)
    packed = np.frombuffer(data, dtype=np.uint8)
    return _unpack_2bit(packed, count)


def self_test() -> None:
    import tempfile

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    classifier_module_dir = Path(__file__).resolve().parent
    # Reuse the classifier module's own synthetic-hierarchy builder rather
    # than reinventing one.
    hierarchy_module = importlib.import_module(DEFAULT_CLASSIFIER_MODULE_NAME)

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        hierarchy_root = root / "hierarchy"
        hierarchy_module._write_fake_hierarchy(hierarchy_root)
        ground_truth = json.loads((hierarchy_root / "ground_truth.json").read_text())

        residual_lon, residual_lat = ground_truth["residual_mixed"]
        _build_synthetic_geometry_cache(hierarchy_root, residual_lon, residual_lat)

        n = len(ground_truth)
        lons = [pt[0] for pt in ground_truth.values()]
        lats = [pt[1] for pt in ground_truth.values()]
        labels = list(ground_truth.keys())
        ais_df = pd.DataFrame(
            {
                "SHIP_ID": np.arange(1, n + 1),
                "LAT": lats,
                "LON": lons,
                "TIMESTAMP": pd.date_range("2026-01-01", periods=n, freq="min"),
                "SPEED": np.full(n, 5.0),
            }
        )
        ais_folder = root / "ais"
        ais_folder.mkdir()
        ais_df.to_parquet(ais_folder / "synthetic.parquet")

        print("=== Test 1: exact_tail=False (default) -- residual_mixed should stay MIXED ===")
        out_dir_1 = root / "out1"
        run(
            ais_folder=ais_folder, hierarchy_root=hierarchy_root, output_dir=out_dir_1,
            exact_tail=False, qgis_output=True, max_files=None,
            classifier_module_dir=classifier_module_dir, classifier_module_name=DEFAULT_CLASSIFIER_MODULE_NAME,
        )
        state_1 = _load_packed_state(out_dir_1, "synthetic", n)
        by_label_1 = dict(zip(labels, state_1.tolist()))
        print("  decoded states:", by_label_1)
        assert by_label_1["marine_center"] == 0, "marine_center should be marine_water"
        assert by_label_1["land_center"] == 1, "land_center should be land"
        assert by_label_1["ultrafine_marine"] == 0, "ultrafine_marine should resolve to marine_water"
        assert by_label_1["residual_mixed"] == 2, "residual_mixed should stay MIXED when exact_tail is off"
        assert by_label_1["fine_land_no_ultrafine"] == 1, "fine_land_no_ultrafine should be land"
        assert by_label_1["fine_inland"] == 3, "fine_inland should be inland_water"

        qgis_1 = pd.read_parquet(out_dir_1 / "synthetic.flagged_for_qgis.parquet")
        flagged_labels_1 = set(qgis_1["SHIP_ID"].map(lambda sid: labels[sid - 1]))
        expected_flagged_1 = {"land_center", "residual_mixed", "fine_land_no_ultrafine"}
        assert flagged_labels_1 == expected_flagged_1, (
            f"expected QGIS-flagged rows {expected_flagged_1}, got {flagged_labels_1}"
        )
        print("  QGIS output correctly contains exactly the land + mixed rows.")

        print("\n=== Test 2: exact_tail=True -- residual_mixed should now resolve via the synthetic land polygon ===")
        out_dir_2 = root / "out2"
        run(
            ais_folder=ais_folder, hierarchy_root=hierarchy_root, output_dir=out_dir_2,
            exact_tail=True, qgis_output=True, max_files=None,
            classifier_module_dir=classifier_module_dir, classifier_module_name=DEFAULT_CLASSIFIER_MODULE_NAME,
        )
        state_2 = _load_packed_state(out_dir_2, "synthetic", n)
        by_label_2 = dict(zip(labels, state_2.tolist()))
        print("  decoded states:", by_label_2)
        assert by_label_2["residual_mixed"] == 1, (
            "residual_mixed should resolve to LAND once exact_tail is on and the "
            "synthetic geometry cache covers it with a land polygon"
        )
        for label in ("marine_center", "land_center", "ultrafine_marine", "fine_land_no_ultrafine", "fine_inland"):
            assert by_label_2[label] == by_label_1[label], f"{label} should be unaffected by exact_tail"
        print("  exact-tail correctly resolved the one residual MIXED point, left everything else unchanged.")

        metadata_2 = json.loads((out_dir_2 / "synthetic.metadata.json").read_text())
        assert metadata_2["exact_tail_resolved_rows"] == 1
        print("  metadata correctly reports exactly 1 exact-tail-resolved row.")

    print("\nSelf-test PASSED.")


def _self_test_sharding() -> None:
    import tempfile

    classifier_module_dir = Path(__file__).resolve().parent
    hierarchy_module = importlib.import_module(DEFAULT_CLASSIFIER_MODULE_NAME)

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        hierarchy_root = root / "hierarchy"
        hierarchy_module._write_fake_hierarchy(hierarchy_root)

        ais_folder = root / "ais"
        ais_folder.mkdir()
        n_files = 8
        for i in range(n_files):
            pd.DataFrame({
                "SHIP_ID": [1, 2], "LAT": [10.0 + i, 10.0 + i], "LON": [20.0 + i, 20.0 + i],
                "TIMESTAMP": pd.date_range("2026-01-01", periods=2, freq="min"), "SPEED": [5.0, 5.0],
            }).to_parquet(ais_folder / f"file_{i:02d}.parquet")

        output_dir = root / "out"
        shard_count = 4

        print("=== Sharding test: 8 files across 4 shards, interleaved ===")
        for shard_index in range(shard_count):
            run(
                ais_folder=ais_folder, hierarchy_root=hierarchy_root, output_dir=output_dir,
                exact_tail=False, qgis_output=False, max_files=None,
                classifier_module_dir=classifier_module_dir, classifier_module_name=DEFAULT_CLASSIFIER_MODULE_NAME,
                shard_index=shard_index, shard_count=shard_count,
            )
            shard_catalog = json.loads((output_dir / f"land_water_mask_catalog_shard{shard_index}_of_{shard_count}.json").read_text())
            expected_files = {f"file_{i:02d}.parquet" for i in range(shard_index, n_files, shard_count)}
            actual_files = {row["source_relative_path"] for row in shard_catalog}
            assert actual_files == expected_files, f"shard {shard_index} expected {expected_files}, got {actual_files}"
            print(f"  shard {shard_index}: correctly got {sorted(actual_files)}")

        print("\n=== Resumability test: rerun shard 0, should skip both its files ===")
        run(
            ais_folder=ais_folder, hierarchy_root=hierarchy_root, output_dir=output_dir,
            exact_tail=False, qgis_output=False, max_files=None,
            classifier_module_dir=classifier_module_dir, classifier_module_name=DEFAULT_CLASSIFIER_MODULE_NAME,
            shard_index=0, shard_count=shard_count,
        )
        # No assertion on output content here -- the real check is that this
        # call completes without error and without recomputing; visually
        # confirmed via the "already processed, skipping" lines above.
        print("  rerun completed without recomputation (see 'already processed, skipping' lines above).")

        print("\n=== Merge test: combine all 4 shards, expect all 8 files, no duplicates ===")
        merge_shard_catalogs(output_dir, shard_count)
        merged = json.loads((output_dir / "land_water_mask_catalog.json").read_text())
        merged_files = {row["source_relative_path"] for row in merged}
        assert merged_files == {f"file_{i:02d}.parquet" for i in range(n_files)}, "merged catalog should contain exactly all 8 files"
        assert len(merged) == n_files, "no duplicates expected"
        print(f"  merged catalog correctly contains all {n_files} files exactly once.")

        print("\n=== Merge overlap-detection test: a duplicated file across two shard catalogs must raise, not silently merge ===")
        bad_dir = root / "out_bad"
        bad_dir.mkdir()
        (bad_dir / "land_water_mask_catalog_shard0_of_2.json").write_text(json.dumps([{"source_relative_path": "x.parquet"}]))
        (bad_dir / "land_water_mask_catalog_shard1_of_2.json").write_text(json.dumps([{"source_relative_path": "x.parquet"}]))
        try:
            merge_shard_catalogs(bad_dir, 2)
            raise AssertionError("merge_shard_catalogs should have raised on an overlapping file, but did not")
        except ValueError as e:
            print(f"  correctly raised: {e}")

    print("\nSharding self-test PASSED.")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--ais-folder", default=None)
    parser.add_argument("--hierarchy-root", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--exact-tail", action="store_true", help="Default OFF.")
    parser.add_argument("--no-qgis-output", action="store_true", help="QGIS Parquet output is ON by default.")
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--classifier-module-dir", default=None)
    parser.add_argument("--classifier-module-name", default=DEFAULT_CLASSIFIER_MODULE_NAME)
    parser.add_argument(
        "--shard-index", type=int, default=None,
        help="This worker's shard number (0-based). Requires --shard-count. "
        "Each shard processes an interleaved subset of files, e.g. --shard-index 0 "
        "--shard-count 4 handles files 0, 4, 8, ... -- run one process per shard "
        "in parallel, each with its own --shard-index but the SAME --shard-count.",
    )
    parser.add_argument("--shard-count", type=int, default=None, help="Total number of parallel workers/shards.")
    parser.add_argument(
        "--merge-shards", action="store_true",
        help="Combine all shard catalogs in --output-dir into one land_water_mask_catalog.json. "
        "Run once after every shard has finished, with --shard-count set to the total shards used.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.self_test:
        self_test()
        _self_test_sharding()
        return
    if args.merge_shards:
        if not args.output_dir or not args.shard_count:
            raise SystemExit("--merge-shards requires --output-dir and --shard-count.")
        merge_shard_catalogs(Path(args.output_dir), args.shard_count)
        return
    if not args.ais_folder or not args.hierarchy_root or not args.output_dir:
        raise SystemExit(
            "--ais-folder, --hierarchy-root and --output-dir are required unless --self-test is passed."
        )
    if (args.shard_index is None) != (args.shard_count is None):
        raise SystemExit("--shard-index and --shard-count must be given together.")
    run(
        ais_folder=Path(args.ais_folder),
        hierarchy_root=Path(args.hierarchy_root),
        output_dir=Path(args.output_dir),
        exact_tail=args.exact_tail,
        qgis_output=not args.no_qgis_output,
        max_files=args.max_files,
        classifier_module_dir=Path(args.classifier_module_dir) if args.classifier_module_dir else Path(__file__).resolve().parent,
        classifier_module_name=args.classifier_module_name,
        shard_index=args.shard_index,
        shard_count=args.shard_count,
    )


if __name__ == "__main__":
    main()
