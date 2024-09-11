# XIC-Extractor

TBD

## Run This script

As long as you have docker/git installed, execute the following:

```shell
nextflow run mpc-bioinformatics/xic-extractor -r main <params>
```

Checkout [this link](https://www.nextflow.io/docs/latest/sharing.html) for more information how to run workflows directly from github.

Alternatively you can just use `git pull` and run it via `nextflow run main.nf`.

## Build docker (to run with a locally build docker image)

```shell
docker build -t xic-extractor:local . -f docker/Dockerfile
```

Run nextflow then with: `-with-docker xic-extractor:local` 

## Parameters

Description of the parameters:

### raw_spectra

Folder to Thermor (.raw) or Bruker (.d) spectra

### trafoxmls

Folder containing the trafoXML (OpenMS) information to correct the retention time (applied after XIC Extraction).

### extraction_csv:

These are all the required columns for the input csv: 

| identifier                            | mz_start                | mz_end                  | rt_start                                              | rt_end                                                | mz                           | ppm                          | ms_level                                    |
|---------------------------------------|-------------------------|-------------------------|-------------------------------------------------------|-------------------------------------------------------|------------------------------|------------------------------|---------------------------------------------|
| Identifier in the final HDF5 (String) | <Either this or mz/ppm> | <Either this or mz/ppm> | Retention time in Minutes (leave empty for no filter) | Retention time in Minutes (leave empty for no filter) | <Either this or mz-start/end | <Either this or mz-start/end | MS level which to extract (ms1, ms2 or ms3) |

### group_regex

Regex (applied on the filename) to be used for grouping (only relevant for plotting). Here the first matching regex group is used as the group.

### export_visualization

Bool (true or false). If set to true, the extracted XICs are plotting as png and html (plotly). Additionally a grouping can be applied.