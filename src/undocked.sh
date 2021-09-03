#!/bin/bash

# Exit on unset variables
set -o nounset
# Exit on unhandled failure in pipes
set -o pipefail

# Print debug message and terminate script on non-zero return codes
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# [1/2] Most versions of Bash read scripts line by line, causing misbehavior if
# the file changes during runtime. The {} forces Bash to read the entire thing
{
	export ANNOVEP_CACHE="${ANNOVEP_CACHE:-${HOME}/annovep}"
	readonly INSTALL="${ANNOVEP_CACHE}/install"

	# VEP exports
	export VEP_ROOT="${VEP_ROOT:-${INSTALL}/vep/src/ensembl-vep}"
	export VEP_PLUGINS="${VEP_PLUGINS:-${INSTALL}/vep-plugins/Plugins}"

	# Misc exports
	export LIFTOVER_CACHE="${LIFTOVER_CACHE:-${ANNOVEP_CACHE}/liftover}"

	. "${INSTALL}/annovep/pipeline.sh" "$@"

	# [2/2] Prevent Bash from reading past this point once script is done
	exit $?
}
