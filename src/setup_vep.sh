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
    . "${ANNOVEP_ROOT}/utilities.sh"

    function download() {
        local -r src_url="${1}"
        local -r dst_file="${2}"
        local -r tmp_file="${2}.${RANDOM}.tmp"

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
        local -r src_url="${1}"
        local -r dst_tar="${VEP_CACHE}/$(basename "${1}")"
        local -r dst_fa="${dst_tar/.tar.gz/.fa.gz}"
        local -r tmp_fa="${dst_fa}.${RANDOM}.tmp"

        download "${src_url}" "${dst_tar}"

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
        local -r src_url="${2}"
        local -r dst_file="${VEP_CACHE}/$1"

        download "${src_url}" "${dst_file}"

        if [ ! -e "${dst_file}.tbi" ]; then
            info "Indexing VCF file ${dst_file}"
            log_command tabix -p vcf "${dst_file}"
        else
            info "Already indexed VCF file ${dst_file}"
        fi
    }

    # Download genome cache
    log_command perl ${VEP_ROOT}/INSTALL.pl \
        --AUTO "cf" \
        --SPECIES "homo_sapiens" \
        --ASSEMBLY "GRCh38" \
        --CACHEDIR "${VEP_CACHE}"

    # FASTA used by the AncestralAllele plugin
    download_ancestral_fa "ftp://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz"

    # ClinVAR VCF used for --custom annotation
    download_vcf "clinvar_20210821.vcf.gz" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20210821.vcf.gz"
    # TODO:
    download_vcf "dbsnp_20180418.vcf.gz" "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz"

    # FASTA used by the LOFTEE plugin
    download "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz" "${VEP_CACHE}/human_ancestor.fa.gz"
    download "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai" "${VEP_CACHE}/human_ancestor.fa.gz.fai"
    download "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi" "${VEP_CACHE}/human_ancestor.fa.gz.gzi"

    # SQLite database used by the LOFTEE plugin
    if [ ! -e "${VEP_CACHE}/phylocsf_gerp.sql" ]; then
        download "https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz" "${VEP_CACHE}/phylocsf_gerp.sql.gz"

        log_command gunzip -k "${VEP_CACHE}/phylocsf_gerp.sql.gz"
    fi

    # TODO
    # download "http://ftp.ensembl.org/pub/current_compara/conservation_scores/90_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw" "${VEP_CACHE}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

    # [2/2] Prevent Bash from reading past this point once script is done
    exit $?
}
