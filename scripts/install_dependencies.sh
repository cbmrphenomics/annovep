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
	apt-get update

	# Needed for wrapper and merging scripts
	apt-get install -y \
		curl \
		cython3 \
		python3 \
		python3-coloredlogs \
		python3-pip \
		python3-pysam

	# Needed to index FASTAs downloaded for VEP
	apt-get install -y \
		samtools

	# Needed for LOFTEE plugin (exact version unclear, so install both)
	apt-get install -y \
		libdbd-sqlite2-perl \
		libdbd-sqlite3-perl

	pip3 install --no-cache-dir \
		liftover

	# Only need to install liftover module
	apt-get autoremove --purge -y \
		cython3

	apt-get autoremove -y
	apt-get clean
	rm -rf /var/lib/apt/lists/*

	# [2/2] Prevent Bash from reading past this point once script is done
	exit $?
}
