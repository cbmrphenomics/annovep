#!/bin/bash

# Exit on unset variables
set -o nounset

# Print debug message and terminate script on non-zero return codes
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# [1/2] Most versions of Bash read scripts line by line, causing misbehavior if
# the file changes during runtime. The {} forces Bash to read the entire thing
{
	function fetch() {
		local filename="${1}"
		local dataset="${2}"

		if [ ! -e "${ANNOVAR_ROOT}/hg38_${filename}.txt" ]; then
			info "Downloading Annovar dataset ${filename} (${dataset})"

			/opt/annovar/annotate_variation.pl \
				-buildver hg38 \
				-downdb \
				-webfrom annovar \
				"$dataset" \
				"${ANNOVAR_CACHE}"
		else
			info "Annovar dataset ${filename} (${dataset}) already downloaded"
		fi
	}

	fetch "gnomad30_genome" "gnomad30_genome"
	fetch "AFR.sites.2015_08" "1000g2015aug"
	fetch "AMR.sites.2015_08" "1000g2015aug"
	fetch "EAS.sites.2015_08" "1000g2015aug"
	fetch "EUR.sites.2015_08" "1000g2015aug"
	fetch "SAS.sites.2015_08" "1000g2015aug"

	# [2/2] Prevent Bash from reading past this point once script is done
	exit $?
}
