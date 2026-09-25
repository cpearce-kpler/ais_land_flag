# ais_land_flag
Flags AIS data that is erroneously on land. These flags can then be used as a mask in other processes. 


# Dependencies
- Downloaded AIS data sorted by latitude and longitude.
- Masks and relate land/water/inland water masks plus related grids.

Masks and grids are available via the geospatial S3 bucket.


# Run in PowerShell
Here is some example PowerShell code to run the script.
>>     & "C:\anaconda\python.exe" -u "C:\Users\John Doe\Desktop\ais_land_water_mask.py" `
>>     --ais-folder "C:\Users\Craig Pearce\Desktop\Data_sets\ais_files\ais_2025" `
>>     --hierarchy-root "C:\Users\Craig Pearce\Desktop\global_land_water_hierarchy_new" `
>>     --output-dir "C:\Users\Craig Pearce\Desktop\ais_2025_land_water_mask_full" `
>>     --exact-tail
