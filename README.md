
# getBiomes


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/get-biomes?include_prereleases&label=version)](https://github.com/genomewalker/get-biomes/releases) [![get-biomes](https://github.com/genomewalker/get-biomes/workflows/getBiomes_ci/badge.svg)](https://github.com/genomewalker/get-biomes/actions) [![PyPI](https://img.shields.io/pypi/v/get-biomes)](https://pypi.org/project/get-biomes/) [![Conda](https://img.shields.io/conda/v/genomewalker/get-biomes)](https://anaconda.org/genomewalker/get-biomes)


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

### Using conda

```bash
conda install -c conda-forge -c bioconda -c genomewalker get-biomes
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

getBiomes will take a file with a list of biomes and will retrieve all the samples belonging to those biomes from MGnify and ENA.

```
$ getBiomes --help

usage: getBiomes -b BIOMES [--mgnify-filter MG_FILTER] [--ena-filter ENA_FILTER] [-p PREFIX]
                 [-t THREADS] [--combine] [--clean] [--debug] [--version] [-h]

A simple tool to get biomes related data from MGnify and ENA

required arguments:
  -b BIOMES, --biomes BIOMES
                        A txt file containing MGnify biomes. Ex: root:Environmental:Aquatic:Marine
                        (default: None)

optional arguments:
  --mgnify-filter MG_FILTER
                        Key-value pairs to filter the MGnify metadata. Valid values are:
                        experiment_type, biome_name, lineage, geo_loc_name, latitude_gte, latitude_lte,
                        longitude_gte, longitude_lte, species, instrument_model, instrument_platform,
                        metadata_key, metadata_value_gte, metadata_value_lte, metadata_value,
                        environment_material, environment_feature, study_accession or include (default:
                        None)
  --ena-filter ENA_FILTER
                        Key-value pairs to filter the ENA metadata. Valid values are: read_count,
                        instrument_model, instrument_platform, library_layout, library_strategy,
                        library_selection or library_source (default: None)
  -p PREFIX, --prefix PREFIX
                        Prefix for the output file (default: None)
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 1)
  --combine             Combine all output files into one (default: False)
  --clean               Remove existing output files (default: False)
  --debug               Print debug messages (default: False)
  --version             Print program version
  -h, --help            show this help message and exit
```

One would run getBiomes as:

```bash
getBiomes -b test-biomes.txt -t 24
```

Where `test-biomes.txt` is a file containing the biomes to retrieve. For example:

```
root:Environmental:Aquatic:Marine
root:Environmental:Aquatic:Freshwater
```

By default, getBiomes will retrieve all the samples from MGnify and ENA that belong to the biomes specified in the input file. However, it is possible to filter the samples retrieved from MGnify using the `--filter` option. For example, to retrieve only the samples from the biomes specified in the input file that have been sequenced using Illumina, that are WGS and with more than 10M reads we can do:

```bash
getBiomes -b test-biomes.txt -t 24 --mgnify-filter '{"instrument_platform":"illumina","metadata_key":"investigation type","metadata_value":"metagenome"}' --ena-filter '{"library_layout":"PAIRED","library_strategy":"WGS","library_source":"METAGENOMIC","library_selection":"RANDOM", "read_count":10000000}'
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
