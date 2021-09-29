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
    readonly ANNOVEP_ROOT="$(dirname "$(readlink -f "$0")")"

    . "${ANNOVEP_ROOT}/utilities.sh"

    MISSING_FILES=0
    function require_files() {
        local -r desc=$1
        shift

        for filename in "${@}"; do
            if [ -e "${filename}" ]; then
                info "[✓] ${desc} file found at ${filename}"
            else
                error "[☓] ${desc} file not found at ${filename}"
                MISSING_FILES=1
            fi
        done
    }

    if [ $# -lt 2 ]; then
        error "Wrong number of arguments (${#}); usage is"
        error "  annovep pipeline <input.vcf.gz> <output_prefix> [<vep args, ...>]"
        exit 1
    fi

    readonly input_vcf="${1}"
    readonly output_prefix="${2}"
    readonly output_vep_json="${2}.vep.json.gz"
    readonly output_vep_html="${2}.vep.html"
    readonly output_tsv="${2}.tsv"
    shift 2

    require_files "Input VCF file" "${input_vcf}"

    ####################################################################################
    ## Paths VEP plugin annotations

    readonly PLUGINS_CACHE="${ANNOVEP_CACHE}/plugins"

    # http://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele
    readonly VEP_ANCESTRAL="${PLUGINS_CACHE}/homo_sapiens_ancestor_GRCh38.fa.gz"
    require_files "Ancestral FASTA" "${VEP_ANCESTRAL}" "${VEP_ANCESTRAL}.fai" "${VEP_ANCESTRAL}.gzi"

    # https://github.com/Ensembl/VEP_plugins/blob/release/104/Conservation.pm
    readonly VEP_CONSERVATION="${PLUGINS_CACHE}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
    require_files "GERP Scores" "${VEP_CONSERVATION}"

    # https://github.com/konradjk/loftee
    readonly VEP_LOFTEE_FA="${PLUGINS_CACHE}/human_ancestor.fa.gz"
    require_files "loftee ancestral FASTA" "${VEP_LOFTEE_FA}" "${VEP_LOFTEE_FA}.fai" "${VEP_LOFTEE_FA}.gzi"

    readonly VEP_LOFTEE_SQL="${PLUGINS_CACHE}/phylocsf_gerp.sql"
    require_files "loftee conservation database" "${VEP_LOFTEE_SQL}"

    readonly VEP_EXACPLI="${PLUGINS_CACHE}/ExACpLI_values.txt"
    require_files "ExACpLI values" "${VEP_EXACPLI}"

    ####################################################################################
    ## Paths and fields for VEP custom annotations

    readonly CUSTOM_CACHE="${ANNOVEP_CACHE}/custom"

    # http://m.ensembl.org/info/docs/tools/vep/script/vep_custom.html#custom_example
    readonly VEP_CLINVAR="${CUSTOM_CACHE}/clinvar_20210821.vcf.gz"
    readonly VEP_CLINVAR_FIELDS="ALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG"
    require_files "ClinVar" "${VEP_CLINVAR}" "${VEP_CLINVAR}.tbi"

    readonly VEP_1K_GENOMES="${CUSTOM_CACHE}/1000Genomes_20200805.vcf.gz"
    readonly VEP_1K_GENOMES_FIELDS="$(printf "AF_%s_unrel," AMR AFR EAS EUR SAS)"
    require_files "1000 genomes (custom)" "${VEP_1K_GENOMES}" "${VEP_1K_GENOMES}.tbi"

    # Custom made gnomAD VCFs containing coverage statistics
    #   $ wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz
    #   $ python3 convert_to_custom gnomad:coverage gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz | bgzip > gnomAD.genomes.r3.0.1.coverage.vcf.gz
    readonly VEP_GNOMAD_COVERAGE="${CUSTOM_CACHE}/gnomAD.genomes.r3.0.1.coverage.vcf.gz"
    readonly VEP_GNOMAD_COVERAGE_FIELDS="gnomAD_mean,gnomAD_median,gnomAD_over_15,gnomAD_over_50"
    require_files "gnomAD coverage (custom)" "${VEP_GNOMAD_COVERAGE}" "${VEP_GNOMAD_COVERAGE}.tbi"

    # Custom made gnomAD VCF containing allele frequencies
    #   $ wget $(printf "https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.chr%s.vcf.gz " $(seq 1 22) X Y)
    #   $ python3 convert_to_custom gnomad:sites gnomad.genomes.r3.0.sites.chr*.vcf.gz | bgzip > gnomad.genomes.r3.0.0.sites.vcf.gz
    readonly VEP_GNOMAD_SITES="${CUSTOM_CACHE}/gnomAD.genomes.r3.0.0.sites.vcf.gz"
    readonly VEP_GNOMAD_SITES_FIELDS="AF,AF_afr,AF_ami,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_sas"
    require_files "gnomAD sites (custom)" "${VEP_GNOMAD_SITES}" "${VEP_GNOMAD_SITES}.tbi"

    # Custom made dbSNP VCF containing aggregated information
    readonly VEP_DBSNP="${CUSTOM_CACHE}/dbsnp_155_20210513_custom.vcf.gz"
    readonly VEP_DBSNP_FIELDS="ids,alts,functions"
    require_files "dbSNP summaries (custom)" "${VEP_DBSNP}" "${VEP_DBSNP}.tbi"

    # Custom made dbSNP VCF containing aggregated information
    readonly VEP_NEIGHBOURS="${CUSTOM_CACHE}/Homo_sapiens.GRCh38.104.neighbours.bed.gz"
    require_files "database of neighbouring genes" "${VEP_NEIGHBOURS}" "${VEP_NEIGHBOURS}.tbi"

    # if require_files failed one or more times
    if [ ${MISSING_FILES} -ne 0 ]; then
        exit 1
    fi

    ####################################################################################

    if [ "${output_vep_json}" -nt "${input_vcf}" ]; then
        info "VEP has already been run on input VCF"
    else
        info "Running VEP on input VCF:"

        # The following are required to the meet expectations of the combine script:
        #   --dont_skip: Ensure that variants on unknown sequences are still written
        #   --allow_non_variant: Ensure that non-variant sites are still written

        log_command perl "${VEP_ROOT}/vep" \
            --verbose \
            --offline \
            --cache \
            --symbol \
            --format "vcf" \
            --json \
            --force_overwrite \
            --compress_output "gzip" \
            --dont_skip \
            --allow_non_variant \
            --dir_cache "${ANNOVEP_CACHE}/cache" \
            --dir_plugins "${VEP_PLUGINS}" \
            --assembly "GRCh38" \
            --plugin "AncestralAllele,${VEP_ANCESTRAL}" \
            --plugin "Conservation,${VEP_CONSERVATION}" \
            --plugin "ExACpLI,${VEP_EXACPLI}" \
            --plugin "LoF,loftee_path:${VEP_PLUGINS},human_ancestor_fa:${VEP_LOFTEE_FA},conservation_file:${VEP_LOFTEE_SQL}" \
            --custom "${VEP_1K_GENOMES},1KGenomes,vcf,exact,0,${VEP_1K_GENOMES_FIELDS}" \
            --custom "${VEP_CLINVAR},ClinVar,vcf,exact,0,${VEP_CLINVAR_FIELDS}" \
            --custom "${VEP_DBSNP},dbSNP,vcf,exact,0,${VEP_DBSNP_FIELDS}" \
            --custom "${VEP_GNOMAD_COVERAGE},gnomAD_coverage,vcf,overlap,0,${VEP_GNOMAD_COVERAGE_FIELDS}" \
            --custom "${VEP_GNOMAD_SITES},gnomAD_sites,vcf,exact,0,${VEP_GNOMAD_SITES_FIELDS}" \
            --custom "${VEP_NEIGHBOURS},neighbours,bed,overlap,0" \
            --polyphen b \
            --input_file "${input_vcf}" \
            --output_file "${output_vep_json}" \
            --stats_file "${output_vep_html}" \
            "${@}"
    fi

    info "Aggregating results ..."
    log_command python3 "${ANNOVEP_ROOT}/vep_json_to_tsv.py" \
        --liftover-cache "${ANNOVEP_CACHE}/liftover" \
        --vep-output "${output_vep_json}" \
        "${input_vcf}" \
        "${output_tsv}"

    # [2/2] Prevent Bash from reading past this point once script is done
    exit $?
}
