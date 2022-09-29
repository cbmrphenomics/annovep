# Custom annotations using VEP

## Installation

The pipeline is designed to be run run using podman (preferred) or docker:

1. To build the docker image, simply run `make`.
2. To run the pipeline, use the `annovep` wrapper script in the `bin` folder

The following environmental variables may be set when running `annovep`:

- `ANNOVEP_IMAGE` the docker image to run.
- `ANNOVEP_CACHE` the location of the annotation files, defaulting to `~/annovep`
- `ANNOVEP_RUNNER` the container manager to use. `podman` is preferred if both `podman` and `docker` is available.

## Setup

To download and build custom annotations, run

```bash
annovep setup
```

The location of the annotations may be changed by setting `ANNOVEP_CACHE` (see above)

Note that this process downloads about 350 GB of data in total and uses about 150 GB of space when done.

## Usage

The follow omits command-line options that are not relevant when running the pipeline through `podman`/`docker`:

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

## Example usage

```bash
annovep pipeline input.vcf.gz output_prefix
```

## Shiny server

A shiny server is provided for browsing and filtering the annotations:

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
