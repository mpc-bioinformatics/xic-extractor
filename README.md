# XIC-Extractor

A simply XIC-Extractor workflow written in Python and nextflow extracting XICs from `*.raw` and `*.d` files and returning hdf5 files. `n`-many XICs are extracted by providing a m/z-window and a retention time window in CSV-format. From there all peaks at a retention time (`m`-many) are aggregated to the maximum (basepeak), and saved to three arrays in the hdf5 file: `retention_times` (`(n, m)`),  `intensities` (`(n, m)`) and `labels` (`(n,)`). The XIC extraction has been implemented, so that thousands of XICs can be extracted in a short time frame.

Additinally this workflow provides an optional visualization step, plotting all XICs across the files and optionally, adjusts the retention time if retention time transformation information is provided (`*.trafoXML`).

This workflow has been written in DSL=2. It can be included into other workflow.

## How to run this workflow

For this workflow, you need `docker`, `git` and `nextflow`.  Follow these commands: 

```shell
git pull https://github.com/mpc-bioinformatics/xic-extractor.git
cd xic-extractor

nextflow run main.nf <params>
```

to get started. Check further below for required and optional parameters.

Alternatively you can execute this repository directly via nextflows git capabilities. Executing:

```shell
nextflow run mpc-bioinformatics/xic-extractor -r main <params>
```

downloads this repository automatically and starts the workflow directly. Checkout [the nextflow documentation](https://www.nextflow.io/docs/latest/sharing.html) for more information how to run workflows directly from GitHub.

### Docker container 

Make sure to have `docker buildx` installed. The docker can be build via the following command:

```shell
docker build -t xic-extractor:local . -f docker/Dockerfile
```

## Parameters

More information about parameters can be found in the source code (especially for python files). Further below are the required and optional parameters for the importable workflow as well es for the standalone workflow.

### Importing the workflow `extract_xics`

This workflow takes: 

* raw_files --> A channel containing only raw_files (e.g. Created via `Channel.fromPath` on a directory)
* xics_to_extract --> A channel containing only a single csv file (e.g. Created via `Channel.fromPath` on a single file)
* trafoXMLs --> Corresponding trafoXML files for the raw_files. Either a channel of .trafoXMLs or `Channel.empty()`

and emits:
* hdf5_files --> generated HDF5 files (as described above) for each added raw_file as a list and
* trafoXMLs --> The added trafoXML files as list (either appended with `NULL___(random)` or the real file if provided).

The trafoXML parameter is optional and can be either empty or could contain the corresponding trafoXML files. This workflow matches these files by the filename to the raw_files.

### Standalone Workflow

**Required parameters**:

* `--raw_spectra` (The path to the raw_spectra folder containing either `.raw` (Thermo) or `.d` (Bruker) files)
* `--extraction_csv` (The csv which needs to contain specific columns for extracting XICs. See further below)

**Optional parameters**:

* `--export_visualization` (A bool, either `true` or `false` (default) to also visualize the extracted XICs)
* `--trafoxmls` (A path to a folder, containing `.trafoXML` files, if RT should be adjusted in visualization)
* `--group_regex` (Visualize groups in the plot with the same color. The first captured regex group is used for grouping. Defaults to `(.*)`, the whole filename)
* `--outdir` (The output path where the HDF5 files and visualizations should be saved)

#### extraction_csv:

The CSV file needs to have the following columns present:

* `identifier` (An identifier for the XIC)
* `mz_start` (The m/z lower/left value, either fill this or provided `mz` and `ppm`)
* `mz_end` (The m/z upper/right value, either fill this or provided `mz` and `ppm`)
* `rt_start` (The RT value in minutes. Leave empty if you want to extract XIC across the whole run)
* `rt_end` (The RT value in minutes. Leave empty if you want to extract XIC across the whole run)
* `mz` (The m/z value where a tolerance interval should be caluclated. Either fill this or provided `mz_start` and `mz_end`)
* `ppm` (The ppm error across `mz`. Either fill this or provided `mz_start` and `mz_end`)
* `ms_level` (Currently can be `ms` for MS1 or `ms2` for MS2. If only Thermo files you can also provided `ms3`)


The `rt_start` and `rt_end` value is optional and does not need to be provided. **BUT** a m/z value need to be provided to extract a XIC. You can either provide `mz_start` and `mz_end` (if the left and right interval boundaries are known) or you can provide `mz` and `ppm`, where the interval is then calculated. An example file is provided for reference:  [`example_extraction.csv`](https://github.com/mpc-bioinformatics/xic-extractor/blob/main/example_extraction.csv)
