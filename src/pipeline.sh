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

    function require_file() {
        if [ -e "${2}" ]; then
            info "[✓] ${1} file found at ${2}"
        else
            error "[☓] ${1} file not found at ${2}"
            exit 1
        fi
    }

    if [ $# -lt 2 -o $# -gt 3 ]; then
        error "Wrong number of arguments (${#}); usage is"
        error "  annovep pipeline <input.vcf.gz> <output_prefix> [<threads>]"
        exit 1
    fi

    readonly input_vcf="${1}"
    readonly output_prefix="${2}"
    readonly output_av_input="${2}.hg38_multianno.input"
    readonly output_av_table="${2}.hg38_multianno.txt"
    readonly output_vep_json="${2}.vep.json.gz"
    readonly output_vep_html="${2}.vep.html"
    readonly output_tsv="${2}.tsv"
    readonly threads=${3:-1}

    require_file "Input VCF file" "${input_vcf}"

    # We manually convert the VCF to the Annovar input format. While table_annovar.pl
    # can do this, enabling VCF input also enables VCF output that we don't need.
    if [ "${output_av_input}" -nt "${input_vcf}" ]; then
        info "Input VCF already converted to Annovar input format"
    else
        info "Converting input VCF to Annovar input format:"

        log_command perl "${ANNOVAR_ROOT}/convert2annovar.pl" \
            --allsample \
            --withfreq \
            --format "vcf4" \
            --outfile "${output_av_input}" \
            "${input_vcf}"
    fi

    if [ "${output_av_table}" -nt "${output_av_input}" ]; then
        info "Annovar has already been run on '${output_av_input}'"
    else
        info "Running annovar on '${output_av_input}':"
        log_command perl "${ANNOVAR_ROOT}/table_annovar.pl" \
            --buildver "hg38" \
            --nastring "." \
            --dot2underline \
            --polish \
            --remove \
            --protocol "gnomad30_genome,AFR.sites.2015_08,AMR.sites.2015_08,EAS.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08" \
            --operation "f,f,f,f,f,f" \
            --outfile "${output_prefix}" \
            "${output_av_input}" \
            "${ANNOVAR_CACHE}" \
            --thread "${threads}"
    fi

    ## Plugin annotations

    # http://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele
    readonly VEP_ANCESTRAL="${VEP_CACHE}/homo_sapiens_ancestor_GRCh38.fa.gz"
    require_file "Ancestral FASTA" "${VEP_ANCESTRAL}"

    # https://github.com/Ensembl/VEP_plugins/blob/release/104/Conservation.pm
    readonly VEP_CONSERVATION="${VEP_CACHE}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
    require_file "GERP Scores" "${VEP_CONSERVATION}"

    # https://github.com/konradjk/loftee
    readonly VEP_LOFTEE_FA="${VEP_CACHE}/human_ancestor.fa.gz"
    require_file "loftee ancestral FASTA" "${VEP_LOFTEE_FA}"
    require_file "loftee ancestral FASTA (FAI)" "${VEP_LOFTEE_FA}.fai"
    require_file "loftee ancestral FASTA (GZI)" "${VEP_LOFTEE_FA}.gzi"

    readonly VEP_LOFTEE_SQL="${VEP_CACHE}/phylocsf_gerp.sql"
    require_file "loftee conservation database" "${VEP_LOFTEE_SQL}"

    ## Custom annotations

    # http://m.ensembl.org/info/docs/tools/vep/script/vep_custom.html#custom_example
    readonly VEP_CLINVAR="${VEP_CACHE}/clinvar_20210821.vcf.gz"
    readonly VEP_CLINVAR_FIELDS="ALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG"
    require_file "ClinVar" "${VEP_CLINVAR}"

    readonly VEP_1K_GENOMES="${VEP_CACHE}/1000Genomes_20200805.vcf.gz"
    readonly VEP_1K_GENOMES_FIELDS="$(printf "AF_%s_unrel," AMR AFR EAS EUR SAS)"
    require_file "1000 genomes (custom)" "${VEP_1K_GENOMES}"

    # Custom made gnomAD VCFs containing coverage statistics
    readonly VEP_GNOMAD_COVERAGE="${VEP_CACHE}/gnomAD_coverage.vcf.gz"
    readonly VEP_GNOMAD_COVERAGE_FIELDS="gnomAD_mean,gnomAD_median,gnomAD_over_15,gnomAD_over_50"
    require_file "gnomAD coverage" "${VEP_GNOMAD_COVERAGE}"

    # Custom made gnomAD VCFs containing SNP filters
    readonly VEP_GNOMAD_SITES="${VEP_CACHE}/gnomad.genomes.r3.0.sites.vcf.gz"
    readonly VEP_GNOMAD_SITES_FIELDS="AF,AF_afr,AF_ami,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_sas"
    require_file "gnomAD sites (custom)" "${VEP_GNOMAD_SITES}"

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
            --compress_output "gzip" \
            --dont_skip \
            --allow_non_variant \
            --dir_cache "${VEP_CACHE}" \
            --dir_plugins "${VEP_PLUGINS}" \
            --plugin "AncestralAllele,${VEP_ANCESTRAL}" \
            --plugin "Conservation,${VEP_CONSERVATION}" \
            --plugin "LoF,loftee_path:${VEP_PLUGINS},human_ancestor_fa:${VEP_LOFTEE_FA},conservation_file:${VEP_LOFTEE_SQL}" \
            --custom "${VEP_1K_GENOMES},1KGenomes,vcf,exact,0,${VEP_1K_GENOMES_FIELDS}" \
            --custom "${VEP_CLINVAR},ClinVar,vcf,exact,0,${VEP_CLINVAR_FIELDS}" \
            --custom "${VEP_GNOMAD_COVERAGE},gnomAD_coverage,vcf,overlap,0,${VEP_GNOMAD_COVERAGE_FIELDS}" \
            --custom "${VEP_GNOMAD_SITES},gnomAD_sites,vcf,exact,0,${VEP_GNOMAD_SITES_FIELDS}" \
            --polyphen b \
            --input_file "${input_vcf}" \
            --output_file "${output_vep_json}" \
            --stats_file "${output_vep_html}" \
            --fork "${threads}"
    fi

    info "Aggregating results ..."
    log_command python3 "${ANNOVEP_ROOT}/combine_results.py" \
        --liftover-cache "${LIFTOVER_CACHE}" \
        --vep-output "${output_vep_json}" \
        "${input_vcf}" \
        "${output_tsv}"

    # [2/2] Prevent Bash from reading past this point once script is done
    exit $?
}
