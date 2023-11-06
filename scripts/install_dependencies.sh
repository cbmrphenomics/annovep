#!/bin/bash

{

# Exit on unset variables
set -o nounset
# Exit on unhandled failure in pipes
set -o pipefail
# Have functions inherit ERR traps
set -o errtrace

# Print debug message and terminate script on non-zero return codes
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

apt-get update

# Needed for wrapper and merging scripts; python3(-minimal) is required for pip
apt-get install -y \
	curl \
	python3.7 \
	python3.7-dev \
	python3-minimal \
	python3-pip \
	vcftools \
	samtools

# Needed for LOFTEE plugin (exact version unclear, so install both)
apt-get install -y \
	libdbd-sqlite2-perl \
	libdbd-sqlite3-perl

apt-get autoremove -y
apt-get clean
rm -rf /var/lib/apt/lists/*

exit $?
}
