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
    readonly PROGNAME="setup_vep"
    readonly ANNOVEP_ROOT="$(dirname "$(readlink -f "$0")")"

    . "${ANNOVEP_ROOT}/utilities.sh"

    function download() {
        local -r dst_file="${1}"
        local -r tmp_file="${1}.${RANDOM}.tmp"
        local -r src_url="${2}"

        if [ ! -e "${dst_file}" ]; then
            info "Downloading $(basename ${dst_file})"

            log_command mkdir -p $(dirname "${tmp_file}")
            if ! log_command curl --fail -L "${src_url}" -o "${tmp_file}"; then
                rm -fv "${tmp_file}"
                exit 1
            fi

            log_command mv "${tmp_file}" "${dst_file}"
        else
            info "Already downloaded ${dst_file}"
        fi

    }

    function download_ancestral_fa() {
        local -r dst_tar="${1}"
        local -r dst_fa="${dst_tar/.tar.gz/.fa.gz}"
        local -r src_url="${2}"
        local -r tmp_fa="${dst_fa}.${RANDOM}.tmp"

        download "${dst_tar}" "${src_url}"

        if [ ! -e "${dst_fa}" ]; then
            info "Unpacking ancestral FASTA to ${dst_fa}"

            # The ancestral genotypes are distributed as a TAR file containing (among other
            # things) per-chromosome FASTA files.
            if ! log_command tar zxf "${dst_tar}" --to-stdout --wildcards '*.fa' | bgzip >"${tmp_fa}"; then
                log_command rm -fv "${tmp_fa}"
                exit 1
            fi

            log_command mv "${tmp_fa}" "${dst_fa}"
        else
            info "Already unpacked ancestral FASTA"
        fi

        if [ ! -e "${dst_fa}.fai" ]; then
            info "Indexing ancestral FASTA file"
            log_command samtools faidx "${dst_fa}"
        else
            info "Ancestral FASTA already indexed"
        fi
    }

    function download_vcf() {
        download "${1}" "${2}"
        download "${1}.tbi" "${2}.tbi"
    }

    ####################################################################################
    ## Main VEP cache

    log_command perl ${VEP_ROOT}/INSTALL.pl \
        --AUTO "cf" \
        --SPECIES "homo_sapiens" \
        --ASSEMBLY "GRCh38" \
        --CACHEDIR "${ANNOVEP_CACHE}/cache"

    ####################################################################################
    ## Plugins cache
    readonly plugins_cache="${ANNOVEP_CACHE}/plugins"

    # FASTA used by the AncestralAllele plugin
    readonly ancestral_fa="${plugins_cache}/homo_sapiens_ancestor_GRCh38.tar.gz"
    readonly ancestral_fa_url="ftp://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz"
    download_ancestral_fa "${ancestral_fa}" "${ancestral_fa_url}"

    # GERP scores for the Conservation plugin
    download "${plugins_cache}/gerp_conservation_scores.homo_sapiens.GRCh38.bw" "http://ftp.ensembl.org/pub/current_compara/conservation_scores/90_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

    # FASTA used by the LOFTEE plugin
    download "${plugins_cache}/human_ancestor.fa.gz" "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz"
    download "${plugins_cache}/human_ancestor.fa.gz.fai" "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai"
    download "${plugins_cache}/human_ancestor.fa.gz.gzi" "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi"

    # SQLite database used by the LOFTEE plugin
    if [ ! -e "${plugins_cache}/phylocsf_gerp.sql" ]; then
        download "${plugins_cache}/phylocsf_gerp.sql.gz" "https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz"

        log_command gunzip -k "${plugins_cache}/phylocsf_gerp.sql.gz"
    fi

    download "${plugins_cache}/ExACpLI_values.txt" "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/104/ExACpLI_values.txt"

    ####################################################################################
    ## Custom annotations
    readonly custom_cache="${ANNOVEP_CACHE}/custom"

    # ClinVAR VCF used for --custom annotation
    download_vcf "${custom_cache}/clinvar_20210821.vcf.gz" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20210821.vcf.gz"

    # Custom dbSNP annotation
    readonly dbsnp_raw="${custom_cache}/dbsnp_155_20210513.vcf.gz"
    readonly dbsnp_custom="${custom_cache}/dbsnp_155_20210513_custom.vcf.gz"
    if [ ! -e "${dbsnp_custom}" -o "${dbsnp_raw}" -nt "${dbsnp_custom}" ]; then
        # DBSNP VCF used for --custom annotation (requires processing)
        download_vcf "${custom_cache}/dbsnp_155_20210513.vcf.gz" "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz"

        readonly dbsnp_tmp="${dbsnp_custom}.${RANDOM}"
        log_command python3 "${ANNOVEP_ROOT}/convert_to_custom.py" dbsnp "${dbsnp_raw}" | bgzip >"${dbsnp_tmp}"
        log_command mv -v "${dbsnp_tmp}" "${dbsnp_custom}"
    else
        info "Already created custom dbSNP database"
    fi

    if [ "${dbsnp_custom}" -nt "${dbsnp_custom}.tbi" ]; then
        log_command tabix -fp vcf "${dbsnp_custom}"
    fi

    # Neighbouring genes from Ensemble gene annotation
    readonly ensembl_raw="${custom_cache}/Homo_sapiens.GRCh38.104.gff3.gz"
    readonly neighbours="${custom_cache}/Homo_sapiens.GRCh38.104.neighbours.bed.gz"
    if [ ! -e "${neighbours}" -o "${ensembl_raw}" -nt "${neighbours}" ]; then
        download "${ensembl_raw}" "http://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz"

        readonly neighbours_tmp="${neighbours}.${RANDOM}"
        log_command python3 "${ANNOVEP_ROOT}/convert_to_custom.py" neighbours "${ensembl_raw}" | bgzip >"${neighbours_tmp}"
        log_command mv -v "${neighbours_tmp}" "${neighbours}"
    else
        info "Already created database of neighbouring genes"
    fi

    if [ "${neighbours}" -nt "${neighbours}.tbi" ]; then
        log_command tabix -fp bed "${neighbours}"
    fi

    readonly ensemble_raw=

    # [2/2] Prevent Bash from reading past this point once script is done
    exit $?
}
