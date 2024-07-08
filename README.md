# NOTICE: The annovep repository has moved to https://github.com/cbmr-data/annovep

# Custom annotations using VEP

## Installation

### Container image

The pipeline is designed to be run run using podman (preferred) or docker:

1. To build the docker image, simply run `make`.
2. To run the pipeline, use the `annovep` wrapper script in the `bin` folder

The following environmental variables may be set when running `annovep`:

- `ANNOVEP_IMAGE` the docker image to run.
- `ANNOVEP_CACHE` the location of the annotation files, defaulting to `~/annovep`
- `ANNOVEP_RUNNER` the container manager to use. `podman` is preferred if both `podman` and `docker` is available.

### Setup

To download and build custom annotations, run

```bash
annovep setup
```

The location of the annotations may be changed by setting `ANNOVEP_CACHE` (see above)

Note that this process downloads about 350 GB of data in total and uses about 150 GB of space when done.
See below for a list of annotations used.

## Usage

To run the pipeline, simply specify the location of a VCF file and (optionally) the prefix used for output files. If no output prefix is provided, the name of the input file is used.

```bash
annovep pipeline input.vcf.gz output_prefix
```

### Commandline options

A number of command-line options are used to specify the location of data, scripts, and plugins. The following therefore excludes command-line options that are automatically set when running the pipeline through `podman`/`docker`:

```
usage: annovep pipeline [-h]
                        [--annotations FILE]
                        [--enable ANNOTATION] [--disable ANNOTATION]
                        [--output-format {tsv,json,sql,sqlite3}]
                        [--include-json] [--fork N] [--buffer-size N]
                        [--log-level {debug,info,warning,error}]
                        in_file [out_prefix]

positional arguments:
  in_file
  out_prefix

optional arguments:
  -h, --help            show this help message and exit
  --annotations FILE    Optional files containing additional annotations
                        (default: [])
  --enable ANNOTATION   Enable annotations disabled by default (default: {})
  --disable ANNOTATION  Disable annotations enabled by default (default: {})

Output:
  --output-format {tsv,json,sql,sqlite3}
                        Output format for aggregated annotations. Maybe be
                        specified zero or more times. Defaults to TSV if not
                        specified (default: [])
  --include-json        Include JSON data in SQL output, excluding sample
                        specific information (default: False)

VEP options:
  --fork N              Use forking to improve VEP runtime (default: None)
  --buffer-size N       Number of VCF records read by VEP per cycle (default:
                        100000)

Logging:
  --log-level {debug,info,warning,error}
                        Log messages at the specified level. This option
                        applies to the `--log-file` option and to log messages
                        printed to the terminal. (default: info)
```

## Shiny server

A shiny server is provided for browsing and filtering annotations written as a SQLite3 database:

1. Use the pipeline to generate a SQLite3 file:

```bash
annovep pipeline input.vcf.gz output_prefix --output-format sqlite3
```

2. Place the sqlite3 file and the contents of the `shiny` folder in a folder:

```bash
ls
bin  output_prefix.db  run.R  server.R  settings.R  ui.R  www
```

3. And edit `settings.R` according to your preferences:

```bash
cat settings.R
# The following settings can be overwritten using a `settings.R` file in the same
# folder as `server.R`. Setting `settings$password` is required to allow logins.

## Load libraries for from external locations here:
# library(flexo, lib="/path/to/R/")

## Title used on shiny page
# settings$title <- "AnnoVEP"

## Password required to view page
settings$password <- "my-password"

## Filename of sqlite3 database
settings$database <- "output_prefix.db"

## Default chromosome to show
# settings$chrom <- NULL

## Default genes
# settings$genes <- NULL

## Vector of column names to show by default
# settings$columns <- c(...)
```

4. Install [shiny](https://shiny.rstudio.com/), [flexo](https://github.com/coolbutuseless/flexo), and run shiny as descibed in the shiny documentation.

## Annotation sources

The following sources of annotations are used. For information about data sources and annotation fields used, refer to `annotations/*.yaml` and `scripts/setup_vep.sh`.

- [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) release 104 with plugins
  - [Ancestral Allele](https://github.com/Ensembl/VEP_plugins/blob/release/104/AncestralAllele.pm)
  - [ExACpLI](https://github.com/Ensembl/VEP_plugins/blob/postreleasefix/104/ExACpLI.pm)
  - [GERP Conservation Scores](https://github.com/Ensembl/VEP_plugins/blob/release/104/Conservation.pm)
  - [LOFTEE](https://github.com/konradjk/loftee) v[1.0.3](https://github.com/konradjk/loftee/releases/tag/v1.0.3)
- Custom annotation based on
  - [1000 Genomes](https://www.internationalgenome.org/) phased hapoltypes [20201028](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased)
  - [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) [20210821](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/)
  - [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) release [155](https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/)
  - [Ensemble](https://www.ensembl.org/) GTF features release [104](https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/)
  - [GnomAD](https://gnomad.broadinstitute.org/) coverage summary r3.0.1
  - [GnomAD](https://gnomad.broadinstitute.org/) sites r3.0
