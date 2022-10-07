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

    function tabix_index() {
        local -r fmt="${1}"
        local -r filename="${2}"

        if [ "${filename}" -nt "${filename}.tbi" ]; then
            log_command tabix -fp "${fmt}" "${filename}"
        fi
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
    download "${plugins_cache}/gerp_conservation_scores.homo_sapiens.GRCh38.bw" "http://ftp.ensembl.org/pub/release-98/compara/conservation_scores/90_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

    # FASTA used by the LOFTEE plugin
    download "${plugins_cache}/human_ancestor.fa.gz" "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz"
    download "${plugins_cache}/human_ancestor.fa.gz.fai" "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai"
    download "${plugins_cache}/human_ancestor.fa.gz.gzi" "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi"

    # SQLite database used by the LOFTEE plugin
    if [ ! -e "${plugins_cache}/phylocsf_gerp.sql" ]; then
        download "${plugins_cache}/phylocsf_gerp.sql.gz" "https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz"

        log_command gunzip -k "${plugins_cache}/phylocsf_gerp.sql.gz"
        log_command rm -v "${plugins_cache}/phylocsf_gerp.sql.gz"
    fi

    download "${plugins_cache}/ExACpLI_values.txt" "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/104/ExACpLI_values.txt"

    ####################################################################################
    ## Custom annotations
    readonly custom_cache="${ANNOVEP_CACHE}/custom"
    readonly custom_sources="${custom_cache}/sources"

    function custom_dbsnp_annotation() {
        local -r dbsnp_raw="${custom_sources}/dbsnp_155_20210513.vcf.gz"
        local -r dbsnp_custom="${custom_cache}/dbsnp_155_20210513_custom.vcf.gz"
        if [ ! -e "${dbsnp_custom}" -o "${dbsnp_raw}" -nt "${dbsnp_custom}" ]; then
            # DBSNP VCF used for --custom annotation (requires processing)
            download_vcf "${dbsnp_raw}" "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz"

            local -r dbsnp_tmp="${dbsnp_custom}.${RANDOM}.tmp"
            log_command python3 "${ANNOVEP_ROOT}/convert_to_custom.py" dbsnp "${dbsnp_raw}" | bgzip >"${dbsnp_tmp}"
            log_command mv -v "${dbsnp_tmp}" "${dbsnp_custom}"
            log_command rm -v "${dbsnp_raw}" "${dbsnp_raw}.tbi"
        else
            info "Already created custom dbSNP database"
        fi

        tabix_index vcf "${dbsnp_custom}"
    }

    function neighbours_from_ensemble() {
        local -r ensembl_raw="${custom_sources}/Homo_sapiens.GRCh38.104.gff3.gz"
        local -r neighbours="${custom_cache}/Homo_sapiens.GRCh38.104.neighbours.bed.gz"
        if [ ! -e "${neighbours}" -o "${ensembl_raw}" -nt "${neighbours}" ]; then
            download "${ensembl_raw}" "http://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz"

            local -r neighbours_tmp="${neighbours}.${RANDOM}.tmp"
            log_command python3 "${ANNOVEP_ROOT}/convert_to_custom.py" neighbours "${ensembl_raw}" | bgzip >"${neighbours_tmp}"
            log_command mv -v "${neighbours_tmp}" "${neighbours}"
            log_command rm -v "${ensembl_raw}"
        else
            info "Already created database of neighbouring genes"
        fi

        tabix_index bed "${neighbours}"
    }

    function gnomad_coverage() {
        local -r gnomad_cov="${custom_cache}/gnomAD.genomes.r3.0.1.coverage.vcf.gz"
        local -r gnomad_cov_raw="${custom_sources}/gnomAD.genomes.r3.0.1.coverage.summary.tsv.bgz"

        if [ ! -e "${gnomad_cov}" -o "${gnomad_cov_raw}" -nt "${gnomad_cov}" ]; then
            download "${gnomad_cov_raw}" "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"

            local -r gnomad_cov_tmp="${gnomad_cov}.${RANDOM}.tmp"
            log_command python3 "${ANNOVEP_ROOT}/convert_to_custom.py" gnomad:cov "${gnomad_cov_raw}" | bgzip >"${gnomad_cov_tmp}"
            log_command mv -v "${gnomad_cov_tmp}" "${gnomad_cov}"
            log_command rm -v "${gnomad_cov_raw}"
        else
            info "Already created custom GnomAD coverage annotations"
        fi

        tabix_index vcf "${gnomad_cov}"
    }

    function gnomad_sites() {
        local -r vcf_final="${custom_cache}/gnomAD.genomes.r3.0.0.sites.vcf.gz"
        if [ ! -e "${vcf_final}" ]; then
            for chrom in $(seq 1 22) X Y; do
                local vcf_part="${custom_sources}/gnomad.genomes.r3.0.sites.chr${chrom}.custom.vcf.gz"

                if [ ! -e "${vcf_part}" ]; then
                    local vcf_part_raw="${custom_sources}/gnomad.genomes.r3.0.sites.chr${chrom}.vcf.gz"
                    download "${vcf_part_raw}" "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.chr${chrom}.vcf.bgz"

                    local vcf_part_tmp="${vcf_part}.${RANDOM}.tmp"
                    log_command python3 "${ANNOVEP_ROOT}/convert_to_custom.py" gnomad:sites "${vcf_part_raw}" | bgzip >"${vcf_part_tmp}"
                    log_command mv -v "${vcf_part_tmp}" "${vcf_part}"
                    log_command rm -v "${vcf_part_raw}"
                else
                    info "Already created custom GnomAD site annotations for chromosome ${chrom}"
                fi
            done

            local -r vcf_tmp="${vcf_final}.${RANDOM}.tmp"
            local -r vcf_parts=$(printf "${custom_sources}/gnomad.genomes.r3.0.sites.chr%s.custom.vcf.gz " $(seq 1 22) X Y)
            log_command vcf-concat ${vcf_parts} | bgzip >"${vcf_tmp}"
            tabix_index vcf "${vcf_tmp}"
            log_command mv -v "${vcf_tmp}" "${vcf_final}"
            log_command mv -v "${vcf_tmp}.tbi" "${vcf_final}.tbi"

            # log_command rm -v ${vcf_parts}
        else
            info "Already created custom GnomAD site annotations"
        fi

        tabix_index vcf "${vcf_final}"
    }

    function population_stats_1kg() {
        local -r vcf_final="${custom_cache}/1000Genomes_20200805.vcf.gz"
        if [ ! -e "${vcf_final}" ]; then
            local -r root_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

            for chrom in $(seq 1 22) X; do
                if [ "${chrom}" = "X" ]; then
                    local vcf_part_raw_name="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.filtered.eagle2-phased.v2.vcf.gz"
                else
                    local vcf_part_raw_name="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
                fi

                local vcf_part_raw="${custom_sources}/${vcf_part_raw_name}"
                local vcf_part="${custom_sources}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.custom.vcf.gz"

                if [ ! -e "${vcf_part}" ]; then
                    download "${vcf_part_raw}" "${root_url}/${vcf_part_raw_name}"

                    local vcf_part_tmp="${vcf_part}.${RANDOM}.tmp"
                    log_command python3 "${ANNOVEP_ROOT}/convert_to_custom.py" 1k_genomes "${vcf_part_raw}" | bgzip >"${vcf_part_tmp}"
                    log_command mv -v "${vcf_part_tmp}" "${vcf_part}"
                    log_command rm -v "${vcf_part_raw}"
                else
                    info "Already created custom 1000 Genomes annotations for chromosome ${chrom}"
                fi
            done

            local -r vcf_tmp="${vcf_final}.${RANDOM}.tmp"
            local -r vcf_parts=$(printf "${custom_sources}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr%s.custom.vcf.gz " $(seq 1 22) X)
            log_command vcf-concat ${vcf_parts} | bgzip >"${vcf_tmp}"
            tabix_index vcf "${vcf_tmp}"
            log_command mv -v "${vcf_tmp}" "${vcf_final}"
            log_command mv -v "${vcf_tmp}.tbi" "${vcf_final}.tbi"
            # log_command rm -v "${vcf_parts}"
        else
            info "Already created custom 1K Genomes annotations"
        fi

        tabix_index vcf "${vcf_final}"
    }

    # ClinVAR VCF used for --custom annotation
    download_vcf "${custom_cache}/clinvar_20210821.vcf.gz" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20210821.vcf.gz"

    # Download dbSNP and generate custom annotation
    custom_dbsnp_annotation

    # Determine neighbouring genes from Ensemble gene annotation
    neighbours_from_ensemble

    # Summarize GnomAD coverage and sites statistics
    gnomad_coverage
    gnomad_sites

    # 1000 genomes population statistics
    population_stats_1kg

    exit $? # [2/2] Prevent Bash from reading past this point once script is done
}
