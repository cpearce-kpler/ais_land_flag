# ais_land_water_mask.py
Flags AIS data that is erroneously on land. These flags can then be used as a mask in other processes. 


## Dependencies
- Downloaded AIS data sorted by latitude and longitude.
- High-resolution land/water/inland water masks plus related cascade grids.

The high-resolution masks and related grids are available via the geospatial S3 bucket:
https://eu-west-1.console.aws.amazon.com/s3/buckets/kp-maritime-assets-geospatial-shared-dev-main?region=eu-west-1&prefix=land_water_mask/&showversions=false


## Run in PowerShell
Here is some example PowerShell code to run the script.
>>     & "C:\anaconda\python.exe" -u "C:\Users\John Doe\Desktop\ais_land_water_mask.py" `
>>     --ais-folder "C:\Users\Craig Pearce\Desktop\Data_sets\ais_files\ais_2025" `
>>     --hierarchy-root "C:\Users\Craig Pearce\Desktop\global_land_water_hierarchy_new" `
>>     --output-dir "C:\Users\Craig Pearce\Desktop\ais_2025_land_water_mask_full" `
>>     --exact-tail


## Purpose
The script provides a standalone land/water classification layer for AIS observations by assigning each source row to one of four states: marine_water, land, mixed, or inland_water. It uses the validated static land/water hierarchy as a three-level lookup cascade at 0.01°, 0.001° and 0.0001° resolution, so classification is based only on each observation's latitude and longitude and does not depend on vessel identity or chronological track order.   ais_land_water_mask.

The primary output is a compact, positionally aligned 2-bit state index for each source AIS file, allowing other processing stages to identify or exclude land and unresolved observations without rereading the full classification geometry. An optional exact-tail operation can resolve remaining mixed points against the hierarchy's pre-extracted geometry cache, while an optional QGIS output writes the land and mixed observations as a smaller review dataset.


## Logic Summary
For each AIS source file, the script reads the required LAT and LON values and passes them through the validated 0.01° → 0.001° → 0.0001° hierarchy cascade. The hierarchy resolves points to marine_water, land, inland_water or mixed; when a point matches more than one exact layer, the precedence is marine water first, inland water overriding marine water, and land only where neither water layer matches.

By default, unresolved mixed points remain mixed, but --exact-tail can optionally resolve them against the hierarchy's cached land, marine and inland geometries using STRtrees (i.e. the high-resolution masks). 


## Output Summary
The resulting state is packed into a 2-bit-per-row binary index aligned exactly with the original AIS row order, with optional QGIS output containing only land or mixed rows for visual review.

A companion JSON metadata file records the source file, row count, state definitions, state counts, compression status and exact-tail processing information. When QGIS output is enabled, the script also writes a smaller Parquet file containing only land and mixed observations, with LAT, LON, SHIP_ID, TIMESTAMP, state and state name for visual review. A top-level catalog.json summarises the processed source files.


## Code Summary
- Imports, configuration and state definitions — Lines 1–104: Loads the required Python libraries and defines the four land/water state codes used throughout the script: marine_water, land, mixed and inland_water. It also identifies the validated StaticLandWaterHierarchy classifier and optional Zstandard/Shapely dependencies.

- 2-bit state packing and output compression — Lines 107–137: Converts the four-state classification into a compact 2-bit-per-row binary representation and provides matching unpacking logic for validation/read-back. The binary state file can optionally be compressed with Zstandard.

- Optional exact-tail geometry cache — Lines 140–246: Defines an optional exact spatial fallback for points that remain mixed after the hierarchy cascade. It loads pre-extracted land, marine and inland-water geometries, decomposes multipart geometries before building STRtrees, and applies the defined marine/inland/land precedence rules.

- Per-file land/water classification — Lines 249–321: Reads the required AIS coordinate columns, applies the three-tier hierarchy classifier, optionally resolves residual mixed points with exact geometry, packs the resulting states into the row-aligned binary index, and writes per-file metadata. When enabled, it also writes a smaller QGIS Parquet containing land and mixed observations.

- Main processing and resumability — Lines 325–420: Loads the static hierarchy once, discovers the AIS Parquet files, optionally divides files into interleaved processing shards, skips files already represented by metadata checkpoints, and processes each remaining file independently. A catalogue records the results of all processed files.

- Shard catalogue merging — Lines 422–445: Combines the catalogues produced by separate file-processing shards into one complete catalogue and explicitly checks for overlapping files so that duplicated processing is not silently accepted.

- Synthetic validation and self-tests — Lines 447–638: Builds synthetic hierarchy and geometry-cache data to test the default three-tier classification, exact-tail resolution, QGIS output and state packing. Additional tests verify interleaved sharding, resumability and duplicate detection when shard catalogues are merged.

- Command-line interface and execution — Lines 639–698: Defines the PowerShell arguments for input/output paths, exact-tail processing, QGIS output, file limits and optional parallel shards, and then dispatches either the normal processing workflow, shard-catalogue merge or self-test.


## Performance and optimisation
The script is designed as a lightweight per-row classification process, using the pre-built land/water hierarchy so that each AIS observation can be classified from its latitude/longitude without reconstructing spatial indexes or rereading the source geometries. Processing is independent by source file, allowing files to be distributed across parallel shards, while only the coordinate columns are read unless QGIS output is enabled. Classification results are packed into a compact 2-bit-per-row binary index, substantially reducing the storage required for the land/water results. The optional exact-tail stage uses a pre-extracted geometry cache and STRtrees, with multipart geometries decomposed beforehand to reduce unnecessary spatial candidates.


## Validation
The script includes both input/process validation and synthetic self-tests to verify that the land/water classification behaves as intended. The self-tests check the three-tier hierarchy results, confirm that residual mixed points are resolved correctly when the optional exact-tail stage is enabled, and verify that QGIS output contains the expected land and mixed observations. Additional tests validate interleaved file sharding, resumability of previously processed files, and detection of duplicate files when shard catalogues are merged.

The outputs were subjectively assessed visually and objectively by using the zones dataset, where flagged points were tested for spatial overlap with berths.

