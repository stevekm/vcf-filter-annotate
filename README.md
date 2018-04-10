# vcf-filter-annotate

Pipeline for filtering and annotating variants in .vcf format

# Install

clone this repo

```
git clone https://github.com/stevekm/vcf-filter-annotate.git
cd vcf-filter-annotate
```

# Usage

## Setup

__NOTE:__ local usage requires Docker to be running

If you do not have ANNOVAR reference databases already downloaded, then first run the following to download them:

```
make annovar
```

## Run

__NOTE:__ Profiles for various configurations have been included in the `nextflow.config` file. One of them must be specified to run.

To run locally:

```
make run EP='-profile local'
```

To run on NYU phoenix compute cluster:

```
make run EP='-profile phoenix'
```

# Software

- Java 8 (Nextflow)

- other software dependencies are packaged in Docker containers (specified in `nextflow.config`)
