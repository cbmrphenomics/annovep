#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import collections
import functools
import gzip
import json
import logging
import re
import sys
from pathlib import Path
from typing import Dict

import coloredlogs
import liftover


########################################################################################


def _build_consequence_ranks():
    # https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    consequences = [
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        # See below for special handling of variants in NMD transcripts
        # "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "regulatory_region_variant",
        "feature_truncation",
        # See below; is ranked below custom NMD terms
        # "intergenic_variant",
    ]

    # Consequences in NMD transcripts are given the lowest priority
    consequences.extend(f"NMD_{_name}" for _name in tuple(consequences))
    consequences.append("NMD_transcript_variant")
    # Intergenic is the absolutely most insignificant term
    consequences.append("intergenic_variant")

    return {consequence: rank for rank, consequence in enumerate(consequences)}


def _build_columns():
    onek_af = "Frequency of existing variant in 1000 Genomes combined {} population"
    gnomad_af = "Frequency of existing variant in gnomAD genomes {} population"

    return {
        "Chr": "Chromosome/Contig recorded in input VCF",
        "Pos": "Position recorded in input VCF",
        "ID": "ID recorded in input VCF",
        "Ref": "Reference allele recorded in input VCF",
        "Alt": "The single ALT allele described by this row",
        "Alts": "The full ALT string from the input VCF",
        "Filters": "Filters recorded in input VCF",
        "DP": "Sum of read depth for this position",
        "Freq": "Frequency of alternative allele in samples",
        "GT_00": "Number of samples with ref/ref genotype",
        "GT_01": "Number of samples with ref/alt genotype",
        "GT_11": "Number of samples with alt/alt genotype",
        "GT_NA": "Number of samples with missing genotypes",
        "GT_other": "Number of samples with other genotypes",
        "Info": "INFO string from input VCF",
        "Liftover_hg19": "Corresponding coordinates in hg19 genome. May be . if "
        "position is invalid or missing from the genome.",
        "Liftover_hg38": "Corresponding coordinates in hg38 genome. May be . if "
        "position is invalid or missing from the genome.",
        "VEP_allele": "The pos:ref:alt corresponding to VEP output",
        "Ancestral_allele": "",
        "Genes_overlapping": "Genes overlapping allele",
        "Genes_upstream": "Neighbouring genes upstream of allele",
        "Genes_downstream": "Neighbouring genes downstream of allele",
        "Func_n_most_significant": "Number of consequences ranked as most"
        "significant in terms of impact",
        "Func_most_significant": "The most significant functional consequence",
        "Func_least_significant": "The last significant functional consequence for "
        "the same gene as the most significant consequence",
        "Func_most_significant_canonical": "The most significant functional "
        "consequence for canonical transcripts only",
        "Func_gene_id": "Gene with the most significant consequence",
        "Func_transcript_id": "Transcript with the most significant consequence",
        "Func_gene_symbol": "Gene symbol (e.g. HGNC)",
        "Func_gene_symbol_source": "Source of gene symbol",
        "Func_cdna_position": "Relative position of base pair in cDNA sequence",
        "Func_cds_position": "Relative position of base pair in coding sequence",
        "Func_protein_position": "Relative position of amino acid in protein",
        "Func_amino_acids": "Reference and variant amino acids",
        "Func_codons": "Reference and variant codon sequence",
        "Func_impact": "Subjective impact classification of consequence type",
        "Func_strand": "Strand of the feature (1/-1)",
        "Func_polyphen": "PolyPhen prediction",
        "Func_polyphen_score": "PolyPhen score",
        "Func_conservation_score": "The conservation score for this site",
        "Func_lof": "Loss-of-function annotation (HC/LC = High/Low Confidence)",
        "Func_lof_filter": "Reason for LoF not being HC",
        "Func_lof_flags": "Possible warning flags for LoF",
        "Func_lof_info": "Info used for LoF annotation",
        "Func_ExACpLI": "Probabililty of a gene being loss-of-function intolerant",
        "dbSNP_ids": "dbSNP ids for alleles alleles matching this pos:ref/alt",
        "dbSNP_alts": "dbSNP allele strings records matching pos:ref/*",
        "dbSNP_functions": "GO terms recorded in dbSNP",
        "ClinVar_ID": "The ClinVar Allele ID",
        "ClinVar_disease": "ClinVar's preferred disease name",
        "ClinVar_significance": "Clinical significance for this single variant",
        "gnomAD_mean": "gnomAD genomes mean coverage for this site",
        "gnomAD_median": "gnomAD genomes median coverage for this site",
        "gnomAD_over_15": "gnomAD genomes fraction with coverage over 15x",
        "gnomAD_over_50": "gnomAD genomes fraction with coverage over 50x",
        "gnomAD_filter": "gnomAD genomes FILTER",
        "gnomAD_AFR_AF": gnomad_af.format("African/American"),
        "gnomAD_AMI_AF": gnomad_af.format("Amish"),
        "gnomAD_AMR_AF": gnomad_af.format("American"),
        "gnomAD_ASJ_AF": gnomad_af.format("Ashkenazi Jewish"),
        "gnomAD_EAS_AF": gnomad_af.format("East Asian"),
        "gnomAD_FIN_AF": gnomad_af.format("Finnish"),
        "gnomAD_NFE_AF": gnomad_af.format("Non-Finnish European"),
        "gnomAD_OTH_AF": gnomad_af.format("other combined"),
        "gnomAD_SAS_AF": gnomad_af.format("South Asian"),
        "1KG_AFR_AF": onek_af.format("African"),
        "1KG_AMR_AF": onek_af.format("American"),
        "1KG_EAS_AF": onek_af.format("East Asian"),
        "1KG_EUR_AF": onek_af.format("European"),
        "1KG_SAS_AF": onek_af.format("South Asian"),
    }


########################################################################################


def abort(line, *args, **kwargs):
    if args or kwargs:
        line = line.format(*args, **kwargs)

    print(line, file=sys.stderr)
    sys.exit(1)


def parse_vcf(line):
    fields = line.rstrip("\r\n").split("\t")
    chr, pos, id, ref, alt, qual, filters, info, *fmt_and_samples = fields

    samples = []
    if fmt_and_samples:
        fmt_keys = fmt_and_samples[0].split(":")
        for sample in fmt_and_samples[1:]:
            samples.append(dict(zip(fmt_keys, sample.split(":"))))

    return {
        "Chr": chr,
        "Pos": int(pos),
        "ID": None if id == "." else id,
        "Ref": ref,
        # . is treated as a actual value, rather than an empty list. This is done so
        # that (limited) information can be retrieved for non-specific variants.
        "Alts": alt.split(","),
        "Qual": float(qual) if qual != "." else None,
        "Filters": [] if filters == "." else filters.split(";"),
        "Info": info,
        "Samples": samples,
    }


@functools.lru_cache()
def parse_vcf_genotypes(genotypes, _re=re.compile(r"[|/]")):
    if genotypes in (None, "./.", ".|."):
        return (None, None)

    result = tuple(int(value) for value in _re.split(genotypes))
    if len(result) != 2:
        raise ValueError(genotypes)

    return result


########################################################################################
## Annotations


class Annotator:
    def __init__(self, liftover_cache=None) -> None:
        self._consequence_ranks = _build_consequence_ranks()
        self._lifter = liftover.get_lifter("hg38", "hg19", liftover_cache)

    def annotate(self, vep):
        row = parse_vcf(vep["input"])
        samples = row.pop("Samples")

        row["Ref"] = self._validate_sequence(row["Ref"], "ACGTN*")
        row["DP"] = self._calculate_depth(samples)

        genotype_counts = self._count_genotypes(samples)
        frequencies = self._calculate_allele_freqs(genotype_counts)

        # Construct the cleaned up alleles / positions used by VEP
        vep_alleles = self._construct_vep_alleles(row)

        for allele_idx, allele in enumerate(row["Alts"], start=1):
            allele = self._validate_sequence(allele, "ACGTN*.")

            copy = dict(row)
            copy["Alt"] = allele
            copy["Freq"] = frequencies.get(allele_idx, ".")

            gt_00 = genotype_counts.get((0, 0), 0)
            gt_01 = genotype_counts.get((0, allele_idx), 0)
            gt_10 = genotype_counts.get((allele_idx, 0), 0)
            gt_11 = genotype_counts.get((allele_idx, allele_idx), 0)
            gt_na = genotype_counts.get((None, None), 0)

            copy["GT_00"] = gt_00
            copy["GT_01"] = gt_01 + gt_10
            copy["GT_11"] = gt_11
            copy["GT_NA"] = gt_na
            copy["GT_other"] = (
                sum(genotype_counts.values(), 0)
                - gt_00
                - gt_01
                - gt_10
                - gt_10
                - gt_11
                - gt_na
            )

            # Cleaned up coordinates/sequences used by VEP
            vep_allele = vep_alleles[allele]
            copy[":vep:"] = vep_allele

            # The position and sequences that VEP reports for this allele
            copy["VEP_allele"] = "{start}:{ref}:{alt}".format(**vep_allele)

            # Add functional annotation
            consequence = self._get_allele_consequence(vep, vep_allele["alt"])
            self._add_gene_info(consequence, copy)
            self._add_ancestral_allele(consequence, copy)

            copy["Func_conservation_score"] = consequence.get("conservation", ".")
            copy["Func_polyphen"] = consequence.get("polyphen_prediction", ".")
            copy["Func_polyphen_score"] = consequence.get("polyphen_score", ".")
            copy["Func_ExACpLI"] = consequence.get("exacpli", ".")

            # add custom annotation
            self._add_custom_annotation(vep, copy)
            self._add_neighbouring_genes(vep, copy)
            self._add_liftover_annotations(vep, copy)

            yield copy

    def _validate_sequence(self, sequence, whitelist):
        sequence = sequence.upper()

        # Don't bother supporting old/weird VCFs
        invalid_characters = set(sequence).difference(whitelist)
        if invalid_characters:
            raise ValueError(sequence)

        return sequence

    def _calculate_depth(self, samples):
        depths = []
        for sample in samples:
            depth = sample.get("DP", ".")
            if depth != ".":
                depths.append(int(depth))

        return sum(depths) if depths else "."

    def _count_genotypes(self, samples):
        counts = collections.defaultdict(int)
        for sample in samples:
            counts[parse_vcf_genotypes(sample.get("GT"))] += 1

        return dict(counts)

    def _calculate_allele_freqs(self, counts):
        allele_counts = collections.defaultdict(int)
        for alleles, count in counts.items():
            if alleles != (None, None):
                for allele in alleles:
                    allele_counts[allele] += count

        frequencies = {}
        total = sum(allele_counts.values())
        for key, value in allele_counts.items():
            frequencies[key] = "{:.4g}".format(value / total)

        return frequencies

    def _construct_vep_alleles(self, row):
        start = row["Pos"]
        ref = row["Ref"]
        alts = row["Alts"]
        if alts == ["."]:
            # VEP has limited support for including records for non-variants
            return {".": {"start": start, "ref": ref, "alt": ref, "alleles": ref}}

        if any(len(ref) != len(allele) for allele in alts):
            # find out if all the alts start with the same base, ignoring "*"
            any_non_star = False
            first_bases = {ref[0]}
            for alt in alts:
                if not alt.startswith("*"):
                    any_non_star = True
                    first_bases.add(alt[0])

            if any_non_star and len(first_bases) == 1:
                start += 1
                ref = ref[1:] or "-"
                alts = list(alts)
                for idx, alt in enumerate(alts):
                    if alt.startswith("*"):
                        alts[idx] = alt
                    else:
                        alts[idx] = alt[1:] or "-"

        # String generated by VEP summarizing ref/alleles in a JSON record
        allele_string = "{}/{}".format(ref, "/".join(alts))

        return {
            vcf_alt: {
                "start": start,
                "ref": ref,
                "alt": vep_alt,
                "alleles": allele_string,
            }
            for vcf_alt, vep_alt in zip(row["Alts"], alts)
        }

    def _get_allele_consequence(self, vep, allele):
        # The JSON record contains transcript, integenic, or no consequences
        transcript_consequences = vep.get("transcript_consequences", ())
        intergenic_consequences = vep.get("intergenic_consequences", ())
        assert not (transcript_consequences and intergenic_consequences), vep

        consequences = []
        canonical_consequences = []
        for consequence in transcript_consequences or intergenic_consequences:
            if consequence["variant_allele"] == allele:
                consequence_terms = consequence["consequence_terms"]
                if "NMD_transcript_variant" in consequence_terms:
                    # terms in Nonsence Mediated Decay variants are of no significance
                    consequence_terms = [
                        term if term == "NMD_transcript_variant" else f"NMD_{term}"
                        for term in consequence_terms
                    ]

                # Gene ID will be missing for intergenetic consequences
                gene = consequence.get("gene_id", ".")
                is_canonical = consequence.get("canonical")

                for term in consequence_terms:
                    entry = (self._consequence_ranks[term], term, gene, consequence)

                    consequences.append(entry)
                    if is_canonical:
                        canonical_consequences.append(entry)

        if not consequences:
            if allele == vep["allele_string"] or allele.startswith("*"):
                # No consequences for non-variable sites
                return {}

            raise ValueError((allele, vep))

        consequences.sort(key=lambda it: it[0])
        canonical_consequences.sort(key=lambda it: it[0])

        # One of the most significant consequences is picked "randomly"
        _, most_significant, gene_id, consequence = consequences[0]
        _, most_significant_canonical, _, _ = consequences[0]

        consequence["most_significant_canonical"] = most_significant_canonical
        consequence["most_significant"] = most_significant
        consequence["n_most_significant"] = 0
        for _, term, _, _ in consequences:
            if term != most_significant:
                break

            consequence["n_most_significant"] += 1

        for _, term, gene_id_, _ in reversed(consequences):
            if gene_id_ == gene_id:
                consequence["least_significant"] = term
                break

        # Effects on NMD transcritps are considered to be the lowest possible
        if most_significant.startswith("NMD_"):
            consequence["impact"] = "MODIFIER"

        return consequence

    def _add_gene_info(self, consequence, dst):
        for key in (
            "amino_acids",
            "codons",
            "cdna_end",
            "cdna_start",
            "gene_id",
            "transcript_id",
            "gene_symbol",
            "gene_symbol_source",
            "impact",
            "strand",
            "least_significant",
            "most_significant",
            "n_most_significant",
            "most_significant_canonical",
            "lof",
            "lof_filter",
            "lof_flags",
            "lof_info",
        ):
            dst[f"Func_{key}"] = consequence.get(key, ".")

        for key in (
            "cdna",
            "cds",
            "protein",
        ):
            dst[f"Func_{key}_position"] = self._format_coordinates(consequence, key)

    def _format_coordinates(self, consequence, key):
        start = consequence.get(f"{key}_start")
        end = consequence.get(f"{key}_end")

        if start is None:
            if end is None:
                return "."

            return f"?-{end}"
        elif end is None:
            return f"{start}-?"
        elif start == end:
            return start

        return f"{start}-{end}"

    def _add_ancestral_allele(self, consequence, dst):
        allele = consequence.get("aa")
        # The explicty check for falsey values is used to catch both missing values and
        # as a workaround for bug where "aa" is -nan (converted to None in _read_record)
        if not allele:
            allele = "." * len(dst["Ref"])

        dst["Ancestral_allele"] = allele

    def _do_add_custom_annotation(self, src, dst, name, fields, default="."):
        data = {}

        alt = dst[":vep:"]["alt"]
        # annotation = {"fields": ..., "name": "chr:start-end"}
        annotations = src.get("custom_annotations", {}).get(name, [])
        for annotation in annotations:
            if annotation.get("allele", alt) == alt:
                # data = {"FILTER": ".", field1: value1, ...}
                data = annotation.get("fields", {})
                break

        if not isinstance(fields, dict):
            fields = dict(zip(fields, fields))

        for src_key, dst_key in fields.items():
            dst[dst_key] = data.get(src_key, default)

    def _add_custom_annotation(self, src, dst):
        self._do_add_custom_annotation(
            src=src,
            dst=dst,
            name="1KGenomes",
            fields={
                "AF_AFR_unrel": "1KG_AFR_AF",
                "AF_AMR_unrel": "1KG_AMR_AF",
                "AF_EAS_unrel": "1KG_EAS_AF",
                "AF_EUR_unrel": "1KG_EUR_AF",
                "AF_SAS_unrel": "1KG_SAS_AF",
            },
        )

        self._do_add_custom_annotation(
            src=src,
            dst=dst,
            name="dbSNP",
            fields={
                "ids": "dbSNP_ids",
                "alts": "dbSNP_alts",
                "functions": "dbSNP_functions",
            },
        )

        self._do_add_custom_annotation(
            src=src,
            dst=dst,
            name="gnomAD_coverage",
            fields=(
                "gnomAD_mean",
                "gnomAD_median",
                "gnomAD_over_15",
                "gnomAD_over_50",
            ),
        )

        self._do_add_custom_annotation(
            src=src,
            dst=dst,
            name="gnomAD_sites",
            fields={
                "FILTER": "gnomAD_filter",
                "AF": "gnomAD_ALL_AF",
                "AF_afr": "gnomAD_AFR_AF",
                "AF_ami": "gnomAD_AMI_AF",
                "AF_amr": "gnomAD_AMR_AF",
                "AF_asj": "gnomAD_ASJ_AF",
                "AF_eas": "gnomAD_EAS_AF",
                "AF_fin": "gnomAD_FIN_AF",
                "AF_nfe": "gnomAD_NFE_AF",
                "AF_oth": "gnomAD_OTH_AF",
                "AF_sas": "gnomAD_SAS_AF",
            },
        )

        self._do_add_custom_annotation(
            src=src,
            dst=dst,
            name="ClinVar",
            fields={
                "ALLELEID": "ClinVar_ID",
                "CLNDN": "ClinVar_disease",
                "CLNSIG": "ClinVar_significance",
            },
        )

    def _add_neighbouring_genes(self, src, dst, nnearest=3):
        # Start coordinate of VEP allele
        astart = src["start"]
        # End coordinate of the allele. This is shared between all ALTs
        aend = src["end"]

        # Alleles may cover multiple bases, so we may have multiple instances of each
        # category, some of which may also be overlapping with the allele
        neighbours = {"downstream": set(), "upstream": set(), "overlap": set()}

        # annotation = {"fields": ..., "name": "chr:start-end"}
        annotations = src.get("custom_annotations", {}).get("neighbours", [])
        for annotation in annotations:
            values = annotation["name"].split(";")

            # The first value (the category) is skipped
            for gene in values[1:]:
                nstart_end, name = gene.split(":")
                nstart, nend = nstart_end.split("-")
                nstart = int(nstart)
                nend = int(nend)

                if aend < nstart:
                    neighbours["downstream"].add((nstart - aend, name))
                elif astart > nend:
                    neighbours["upstream"].add((astart - nend, name))
                else:
                    neighbours["overlap"].add(name)

        def _to_list(values):
            values = sorted(values)[:nnearest]
            if not values:
                return "."

            return ";".join(f"{distance}:{name}" for distance, name in values)

        dst["Genes_overlapping"] = ";".join(sorted(neighbours["overlap"])) or "."
        dst["Genes_upstream"] = _to_list(neighbours["upstream"])
        dst["Genes_downstream"] = _to_list(neighbours["downstream"])

    def _add_liftover_annotations(self, vep, row):
        src_chrom = row["Chr"]
        src_pos = row["Pos"]
        src_liftover = "{}:{}+".format(src_chrom, src_pos)

        # Returns list of overlapping liftover coordinates, an empty list if the
        # position does not exist in the target genome, or KeyError if unknown.
        try:
            coordinates = self._lifter.query(src_chrom, int(src_pos))
        except KeyError:
            coordinates = []
            # Src coordinates were not usable, so we also write them as missing
            src_liftover = "."

        positions = []
        for (chrom, pos, strand) in coordinates:
            positions.append("{}:{}{}".format(chrom, pos, strand))

        dst_liftover = ";".join(positions) or "."

        row[f"Liftover_hg19"] = dst_liftover
        row[f"Liftover_hg38"] = src_liftover


class Output:
    def __init__(self, keys, out_prefix, extension):
        self.keys = keys
        self._handle = sys.stdout
        if out_prefix is not None:
            self._handle = open(f"{out_prefix}{extension}", "wt")

    def finalize(self):
        self._handle.close()

    def process_row(self, data):
        raise NotImplementedError()

    def _print(self, line, *args, **kwargs):
        if args or kwargs:
            line = line.format(*args)

        print(line, file=self._handle, **kwargs)


class JSONOutput(Output):
    def __init__(self, keys, out_prefix):
        super().__init__(keys, out_prefix, ".json")

    def process_row(self, data):
        json.dump({key: data[key] for key in self.keys}, self._handle)
        self._handle.write("\n")


class TSVOutput(Output):
    def __init__(self, keys, out_prefix):
        super().__init__(keys, out_prefix, ".tsv")

        self._print("#{}", "\t".join(self.keys))

        if out_prefix is not None:
            with open(f"{out_prefix}.tsv.columns", "wt") as handle:
                print("Name\tDescription", file=handle)

                for name, description in self.keys.items():
                    print(name, description, sep="\t", file=handle)

    def process_row(self, data):
        row = []
        for key in self.keys:
            value = data[key]
            if isinstance(value, (tuple, list)):
                if value:
                    value = ",".join(map(str, value))
                else:
                    value = "."
            elif value is None:
                value = "."
            else:
                value = str(value)

            row.append(value)

        self._print("\t".join(row))


class SQLOutput(Output):
    def __init__(self, keys, out_prefix):
        super().__init__(keys, out_prefix, ".sql")

        self._row = 0
        for table in ("Annotations",):
            self._print("DROP TABLE IF EXISTS [{}];", table)

        self._print("CREATE TABLE [Annotations] (")
        self._print("    pid INTEGER PRIMARY KEY ASC", end="")
        for key in self.keys:
            self._print(",\n    [{}]", key, end="")
        self._print("\n);")

    def finalize(self):
        self._handle.close()

    def process_row(self, data):
        self._row += 1
        values = [str(self._row)]
        for key in self.keys:
            value = data[key]
            if isinstance(value, (int, float)):
                value = repr(value)
            else:
                value = "'{}'".format(str(value).replace("'", "''"))

            values.append(value)

        self._print("INSERT INTO [Annotations] VALUES ({});", ", ".join(values))


OUTPUT_FORMATS = {
    "json": JSONOutput,
    "tsv": TSVOutput,
    "sql": SQLOutput,
}


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument(
        "in_json",
        type=Path,
    )

    parser.add_argument(
        "out_prefix",
        nargs="?",
        type=Path,
    )

    parser.add_argument(
        "--output-format",
        action="append",
        type=str.lower,
        choices=OUTPUT_FORMATS.keys(),
        help="Output format for aggregated annotations. Maybe be specified zero or "
        "more times. Defaults to TSV if not specified",
    )

    parser.add_argument(
        "--liftover-cache",
        type=Path,
    )

    parser.add_argument(
        "--log-level",
        default="info",
        choices=("debug", "info", "warning", "error"),
        type=str.lower,
        help="Log messages at the specified level. This option applies to the "
        "`--log-file` option and to log messages printed to the terminal.",
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    coloredlogs.install(
        level=args.log_level,
        fmt="%(asctime)s %(levelname)s %(message)s",
    )

    log = logging.getLogger("main")
    log.info("reading variants from '%s'", args.in_json)

    with gzip.open(args.in_json, "rb") as handle:
        header = _build_columns()
        annotator = Annotator(liftover_cache=args.liftover_cache)

        output_formats = set(args.output_format) if args.output_format else ["tsv"]
        if args.out_prefix is None and len(output_formats) > 1:
            log.error("[out_prefix] must be set when writing more than one format")
            return 1

        writers: Dict[str, Output] = {}
        for key in output_formats:
            cls = OUTPUT_FORMATS[key]
            writers[key] = cls(keys=header, out_prefix=args.out_prefix)

        nan_re = re.compile(rb":(-)?NaN\b", flags=re.I)

        for line in handle:
            # Workaround for non-standard JSON output observed in some records, where
            # an expected string value was -nan. Python accepts "NaN", but null seems
            # more reasonable for downstream compatibility
            line = nan_re.sub(b":null", line)

            record = json.loads(line)

            for row in annotator.annotate(record):
                for writer in writers.values():
                    writer.process_row(row)

        for writer in writers.values():
            writer.finalize()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
