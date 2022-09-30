#!/bin/bash

# Exit on unset variables
set -o nounset
# Exit on unhandled failure in pipes
set -o pipefail
# Have functions inherit ERR traps
set -o errtrace

# Print debug message and terminate script on non-zero return codes
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# [1/2] Most versions of Bash read scripts line by line, causing misbehavior if
# the file changes during runtime. The {} forces Bash to read the entire thing
{
    readonly PROGNAME="setup_cpanm"

    . "$(dirname "$(readlink -f "$0")")/utilities.sh"

    info "Creating 'cpanm' directory in $PWD:"
    info "Execute \"export PERL5LIB=\$PERL5LIB:${PWD}/cpanm/lib/perl5\" to activate it"

    mkdir -p cpanm
    export PERL5LIB=$PERL5LIB:$PWD/cpanm/lib/perl5

    info "Installing JSON"
    cpanm -l $PWD/cpanm JSON

    info "Set::IntervalTree"
    cpanm -l $PWD/cpanm Set::IntervalTree

    info "Creating 'deps' directory in $PWD:"
    mkdir -p "deps"

    {
        cd "deps"

        wget "https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2"
        unp "htslib-1.13.tar.bz2"

        {
            cd "htslib-1.13"

            ./configure
            make -j8

            cd ..
        }

        wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
        tar xzf v335_base.tar.gz

        {
            cd "kent-335_base/src/lib"

            export CFLAGS="-fPIC"
            make -j8

            cd ../../..
        }

        cd ..
    }

    cpanm -l $PWD/cpanm Bio::DB::BigFile
    # << deps/kent-335_base/src/lib

    cpanm -l $PWD/cpanm --configure-args="--htslib $PWD/htslib-1.13" Bio::DB::HTS::Tabix

    # [2/2] Prevent Bash from reading past this point once script is done
    exit $?
}
