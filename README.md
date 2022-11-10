
# getBiomes


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/get-biomes?include_prereleases&label=version)](https://github.com/genomewalker/get-biomes/releases) [![get-biomes](https://github.com/genomewalker/get-biomes/workflows/getBiomes_ci/badge.svg)](https://github.com/genomewalker/get-biomes/actions) [![PyPI](https://img.shields.io/pypi/v/get-biomes)](https://pypi.org/project/get-biomes/)


A simple tool to get biomes related data from MGnify and ENA

# Installation

We recommend having [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First, we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/get-biomes/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install get-biomes
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/genomewalker/get-biomes.git
```

By cloning in a dedicated conda environment

```bash
git clone git@github.com:genomewalker/get-biomes.git
cd get-biomes
conda env create -f environment.yml
conda activate get-biomes
pip install -e .
```


# Usage

getBiomes will take a file with a list of biomes and will retrieve all the samples belonging to those biomes from MGnify and ENA. It will then download the samples and create a directory structure to store the samples.

```
$ getBiomes --help

usage: getBiomes [-h] [--version] {search,download} ...

A simple tool to get biomes related data from MGnify and ENA

positional arguments:
  {search,download}  positional arguments
    search           Search for biomes in MGnify and ENA
    download         Download the biomes gathered by the subcommand search

optional arguments:
  -h, --help         show this help message and exit
  --version          Print program version
```

## Search


```
$ getBiomes search --help

usage: getBiomes search [-h] [--debug] -b BIOMES [--mgnify-filter MG_FILTER] [--ena-filter ENA_FILTER]
                        [--exclude-terms EXCLUDE_TERMS] [-p PREFIX] [-t THREADS] [--combine] [--clean]

optional arguments:
  -h, --help            show this help message and exit
  --debug               Print debug messages

search required arguments:
  -b BIOMES, --biomes BIOMES
                        A txt file containing MGnify biomes. Ex: root:Environmental:Aquatic:Marine

search optional arguments:
  --mgnify-filter MG_FILTER
                        Key-value pairs to filter the MGnify metadata. Valid values are:
                        experiment_type, biome_name, lineage, geo_loc_name, latitude_gte, latitude_lte,
                        longitude_gte, longitude_lte, species, instrument_model, instrument_platform,
                        metadata_key, metadata_value_gte, metadata_value_lte, metadata_value,
                        environment_material, environment_feature, study_accession or include
  --ena-filter ENA_FILTER
                        Key-value pairs to filter the ENA metadata. Valid values are: read_count,
                        instrument_model, instrument_platform, library_layout, library_strategy,
                        library_selection or library_source
  --exclude-terms EXCLUDE_TERMS
                        A comma-separated list of terms to exclude from the metadata
  -p PREFIX, --prefix PREFIX
                        Prefix for the output file
  -t THREADS, --threads THREADS
                        Number of threads to use
  --combine             Combine all output files into one
  --clean               Remove existing output files
```


One would run getBiomes search as:

```bash
getBiomes search -b test-biomes.txt -t 24
```

Where `test-biomes.txt` is a file containing the biomes to retrieve. For example:

```
root:Environmental:Aquatic:Marine
root:Environmental:Aquatic:Freshwater
```

By default, getBiomes will retrieve all the samples from MGnify and ENA that belong to the biomes specified in the input file. However, it is possible to filter the samples retrieved from MGnify using the `--filter` option. For example, to retrieve only the samples from the biomes specified in the input file that have been sequenced using Illumina, that are WGS and with more than 10M reads. In addition, we will remove any entry that contains the words `human` and `16S`:

```bash
getBiomes -b test-biomes.txt -t 24 --mgnify-filter '{"instrument_platform":"illumina","metadata_key":"investigation type","metadata_value":"metagenome"}' --ena-filter '{"library_layout":"PAIRED","library_strategy":"WGS","library_source":"METAGENOMIC","library_selection":"RANDOM", "read_count":10000000}' --exclude-terms human,16S
```

The output file will contain the following columns:

```
accession
sample_accession
sample_name
longitude
latitude
geo_loc_name
studies
biome
sample_desc
environment_biome
environment_feature
environment_material
study_accession
experiment_accession
run_accession
read_count
instrument_model
instrument_platform
library_layout
library_strategy
library_selection
library_source
fastq_ftp
query_biome
```

## Download

```
$ getBiomes download --help

usage: getBiomes download [-h] [--debug] -i INPUT [-o OUTDIR] [--clean] [-t THREADS]

optional arguments:
  -h, --help            show this help message and exit
  --debug               Print debug messages (default: False)

download required arguments:
  -i INPUT, --input INPUT
                        A txt file containing MGnify biomes. Ex: root:Environmental:Aquatic:Marine
                        (default: None)

download optional arguments:
  -o OUTDIR, --outdir OUTDIR
                        The directory to save the fastq files (default: None)
  --clean               Remove existing output files (default: False)
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 1)
```

Once the samples have been retrieved, one can download them using the subcommand `download`. For example:

```bash
getBiomes download -i test-biomes__combined.tsv -o test-biomes -t 24
```

Where `test-biomes__combined.tsv` is the output file from the subcommand `search` and `test-biomes` is the directory where the samples will be downloaded. The output directory contains the file `download_report.tsv` with the status of the downloaded files. One can continue downloading the samples that failed in a previous run. If `--clean` is specified, the output directory will be removed before downloading the samples. 

