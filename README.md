# XCS Source.csv Generator

This code is a stripped down version of the MSL pipeline used by the XMM Cluster Survey (XCS). It is a lightweight version, and as such completes quicker, but without producing ancillary files such as Exposure.csv, Detections.csv etc.
___
## Prerequisites
This code requires python 3.7.2+ . For installing on the Apollo2 cluster at the University of Sussex:
```bash
module load python/3.7.2-GCCcore-8.2.0-Extra/3.7.2-GCC8.2.0-Extra

```
### Python Libraries
This module requires some additional python packages:

* astropy
* pandas
* tqdm
* numpy
* PyYAML

The exact version requirements can be found in *requirements.txt*.

To install the packages globally:

```bash
sudo -H pip3 install -r requirements.txt
```

Or to install packages locally, replace package_directory as needed:
```bash
 pip3 install --target=<package_directory> -r requirements.txt
```
___
## Configuring
The entire package is controlled by the *MSL_config.yaml*. There are only three variables you need to change - please do not change the header info files as these are set for the XCS data structure and format.

* data_path: This should point to the directory which contains the XMM observation directories where the folder is the 10 character XMM obsid.
* id_results_loc: This should point to the id_results directory which contain the outputs from the xapa runs.
* fits_path: This is the path from data_path which points to the keVmerged_img. Instead of the obsid in the path you **must** use {obs}.
___
## Running

```bash
python3 make_msl.py
```

___

## Expected output
The script will display multiple loading bars which denote the progress of the code.
The script will create three output files:

* all_finalsigs.dat: This is the aggregate of all the final_sigs files for each obsid.
* all_srcsums.dat: This is the aggregate of all the srcsums files for each obsid.
* Source.csv: The MSL - this is your piece de resistance!!


## Using the module
Most if not all methods have a docstring describing the process. This code base currently needs updating to provide help() functionality.

In the meantime these are the steps required if integrating the function into another script.


```python
msl = MSL(debug=True)
msl.process_idresults_files()
msl.aggregate_srcsums()
msl.aggregate_finalsigs()
msl.create_master_df()
msl.set_column_float(column = "SOURCE_COUNTS")
msl.sort_master_by(column = "SOURCE_COUNTS")
msl.set_master_ID()
msl.flag_source_type()
msl.remove_duplicates()
msl.sort_master_by(column = "OBSID")
msl.log_sig()
msl.make_source_csv()
msl.write_master_to_csv("./MSL.csv")
```

**You must keep debug=True**, non debug support is not ready - this just runs the code verbose unles you supply debug_size=<int>, in which case the script will run and reduce the number of sources parsed on loading to <int>.
