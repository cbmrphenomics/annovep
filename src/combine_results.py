#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import collections
import gzip
import json
import logging
import re
import sys
import textwrap
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import coloredlogs
import liftover
import pysam


# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
VEP_CONSEQUENCES = [
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
    "NMD_transcript_variant",
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
    "intergenic_variant",
]

VEP_CONSEQUENCE_RANKS = {
    consequence: rank for rank, consequence in enumerate(VEP_CONSEQUENCES)
}


def abort(line, *args, **kwargs):
    if args or kwargs:
        line = line.format(*args, **kwargs)

    print(line, file=sys.stderr)
    sys.exit(1)


class Annotator:
    def annotate(self, vcf, row):
        raise NotImplementedError()

    def keys(self):
        raise NotImplementedError()


class AnnotateBasicsInfo(Annotator):
    def __init__(self, keepindelref=False) -> None:
        self._keepindelref = keepindelref

    def annotate(self, vcf, row):
        row["Chr"] = vcf.contig
        row["Pos"] = vcf.pos
        row["ID"] = vcf.id or "."
        row["Ref"] = self._validate_sequence(vcf.ref, "ACGTN*")
        row["Filters"] = ";".join(vcf.filter or ".")
        row["DP"] = self._calculate_depth(vcf)

        # There is no easy way to get just the INFO field a str
        vcf_as_text = str(vcf).split("\t", 8)
        # For VCF files with no genotypes, the INFO field may contain trailing newlines
        row["Info"] = vcf_as_text[7].rstrip()

        genotype_counts = self._calculate_genotype_counts(vcf)
        frequencies = self._calculate_allele_freqs(genotype_counts)

        # Workaround for non-variants; additional logic is found in the VEP annotator
        alts = ["."] if vcf.alts is None else vcf.alts

        # Construct the cleaned up alleles / positions used by VEP
        vep_alleles = self._construct_vep_alleles(vcf)

        for allele in alts:
            allele = self._validate_sequence(allele, "ACGTN*.")

            copy = dict(row)
            copy["Alt"] = allele
            copy["Freq"] = frequencies.get(allele, ".")

            gt_00 = genotype_counts.get((vcf.ref, vcf.ref), 0)
            gt_01 = genotype_counts.get((vcf.ref, allele), 0)
            gt_10 = genotype_counts.get((allele, vcf.ref), 0)
            gt_11 = genotype_counts.get((allele, allele), 0)
            gt_na = genotype_counts.get((None, None), 0)

            copy["GT_00"] = gt_00
            copy["GT_01"] = gt_01 + gt_10
            copy["GT_11"] = gt_11
            copy["GT_NA"] = gt_na
            copy["GT_Other"] = (
                sum(genotype_counts.values(), 0)
                - gt_00
                - gt_01
                - gt_10
                - gt_10
                - gt_11
                - gt_na
            )

            # Cleaned up coordinates/sequences used by VEP
            copy[":vep:"] = vep_alleles[allele]

            yield copy

    def keys(self):
        return {
            "Chr": "Chromosome/Contig recorded in input VCF",
            "Pos": "Position recorded in input VCF",
            "ID": "ID recorded in input VCF",
            "Ref": "Reference allele recorded in input VCF",
            "Alt": "Single alt. allele recorded in input VCF",
            "Filters": "Filters recorded in input VCF",
            "DP": "Sum of read depth for this position",
            "Freq": "Frequency of alternative allele in samples",
            "GT_00": "Number of samples with ref/ref genotype",
            "GT_01": "Number of samples with ref/alt genotype",
            "GT_11": "Number of samples with alt/alt genotype",
            "GT_NA": "Number of samples with missing genotypes",
            "GT_Other": "Number of samples with other genotypes",
            "Info": "INFO string from input VCF",
        }

    def _validate_sequence(self, sequence, whitelist):
        sequence = sequence.upper()

        # Don't bother supporting old/weird VCFs
        invalid_characters = set(sequence).difference(whitelist)
        if invalid_characters:
            raise ValueError(sequence)

        return sequence

    def _calculate_depth(self, vcf):
        depths = []
        for sample in vcf.samples.values():
            depth = sample["DP"]
            if depth is not None:
                depths.append(depth)

        return sum(depths) if depths else "."

    def _calculate_genotype_counts(self, vcf):
        counts = collections.defaultdict(int)
        for sample in vcf.samples.values():
            counts[sample.alleles] += 1

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

    def _construct_vep_alleles(self, vcf):
        start = vcf.pos
        ref = vcf.ref
        if vcf.alts is None:
            # VEP has limited support for including records for non-variants
            return {".": {"start": start, "ref": ref, "alt": ref, "alleles": ref}}

        alts = list(vcf.alts)
        if any(len(vcf.ref) != len(allele) for allele in alts):
            # find out if all the alts start with the same base, ignoring "*"
            first_bases = {ref[0]}
            for alt in alts:
                if not alt.startswith("*"):
                    first_bases.add(alt[0])

            if len(first_bases) == 1:
                start += 1
                ref = ref[1:] or "-"
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
            for vcf_alt, vep_alt in zip(vcf.alts, alts)
        }


class AnnotateLiftOver(Annotator):
    def __init__(self, source: str, cache: Path) -> None:
        if source == "hg19":
            destination = "hg38"
        elif source == "hg38":
            destination = "hg19"
        else:
            raise ValueError(source)

        self._src_db = source.lower()
        self._dst_db = destination.lower()
        self._lifter = liftover.get_lifter(source, destination, cache)

    def annotate(self, vcf, row):
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

        row[f"Liftover_{self._src_db}"] = src_liftover
        row[f"Liftover_{self._dst_db}"] = dst_liftover

        return [row]

    def keys(self):
        return {
            "Liftover_hg19": "Corresponding coordinates in hg19 genome. May be . if "
            "position is invalid or missing from the genome.",
            "Liftover_hg38": "Corresponding coordinates in hg38 genome. May be . if "
            "position is invalid or missing from the genome.",
        }


class AnnotateVEP(Annotator):
    CACHE_SIZE = 1000

    def __init__(self, filepath: Path) -> None:
        self._log = logging.getLogger("vep")
        self._log.info("reading VEP annotations from '%s'", filepath)

        self._handle = gzip.open(filepath, "rt")

        self._cached_for = None
        self._cached_record = {}

    def annotate(self, vcf, row):
        # https://m.ensembl.org/info/docs/tools/vep/vep_formats.html#json
        vep = self._read_record(vcf, row)

        # The position and sequences that VEP reports for this allele
        row["VEPAllele"] = "{start}:{ref}:{alt}".format(**row[":vep:"])

        # Add functional annotation
        consequence = self._get_allele_consequence(vep, row[":vep:"]["alt"])
        self._add_gene_info(consequence, row)
        self._add_ancestral_allele(consequence, row)

        row["Func_conservation_score"] = consequence.get("conservation", ".")
        row["Func_polyphen"] = consequence.get("polyphen_prediction", ".")
        row["Func_polyphen_score"] = consequence.get("polyphen_score", ".")

        # add custom annotation
        self._add_1k_genomes_annotation(vep, row)
        self._add_dbsnp_annotation(vep, row)
        self._add_gnomad_annotation(vep, row)
        self._add_clinvar_annotation(vep, row)

        return [row]

    def keys(self):
        onek_af = "Frequency of existing variant in 1000 Genomes combined {} population"
        gnomad_af = "Frequency of existing variant in gnomAD genomes {} population"

        return {
            "VEPAllele": "The pos:ref:alt corresponding to VEP output",
            "AncestralAllele": "",
            "Func_n_most_significant": "Number of consequences ranked as most"
            "significant in terms of impact",
            "Func_most_significant": "The most significant functional consequence",
            "Func_least_significant": "The last significant functional consequence for "
            "the same gene as the most significant consequence",
            "Func_gene_id": "Gene with the most significant consequence",
            "Func_transcript_id": "Transcript with the most significant consequence",
            "Func_gene_symbol": "Gene symbol (e.g. HGNC)",
            "Func_gene_symbol_source": "Source of gene symbol",
            "Func_all_gene_symbols": "Gene symbols for all genes affected by variant",
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
            "dbSNP_ids": "dbSNP ids for alleles alleles matching this pos:ref/alt",
            "dbSNP_alts": "dbSNP allele strings records matching pos:ref/*",
            "dbSNP_functions": "GO terms recorded in dbSNP",
            "ClinVar_ID": "The ClinVar Allele ID",
            "ClinVar_Disease": "ClinVar's preferred disease name",
            "ClinVar_Significance": "Clinical significance for this single variant",
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

    def _read_record(self, vcf, row, _regex=re.compile(r":(-)?NaN\b", flags=re.I)):
        # There is a single VEP associated with each VCF record (provided that the
        # --allow_non_variant option used. Thus we only need to read a new record when
        # a new VCF record is passed.
        if self._cached_for is not vcf:
            line = self._handle.readline()
            if not line:
                raise RuntimeError("vep annotations are truncated")

            # Workaround for non-standard JSON output observed in some records, where
            # an expected string value was -nan. Python accepts "NaN", but null seems
            # more reasonable
            line = _regex.sub(":null", line)

            self._cached_record = json.loads(line)
            self._cached_for = vcf

            expected_pos = (row["Chr"], row[":vep:"]["start"], row[":vep:"]["alleles"])
            observed_pos = (
                self._cached_record["seq_region_name"],
                self._cached_record["start"],
                self._cached_record["allele_string"],
            )

            # VEP output is expected to be in the same order as the input
            if expected_pos != observed_pos:
                raise ValueError(f"VEP: {observed_pos} != {expected_pos}")

        return self._cached_record

    def _get_allele_consequence(self, vep, allele):
        # The JSON record contains transcript, integenic, or no consequences
        transcript_consequences = vep.get("transcript_consequences", ())
        intergenic_consequences = vep.get("intergenic_consequences", ())
        assert not (transcript_consequences and intergenic_consequences), vep

        consequences = []
        gene_symbols = set()
        for consequence in transcript_consequences or intergenic_consequences:
            if consequence["variant_allele"] == allele:
                for term in consequence["consequence_terms"]:
                    rank = VEP_CONSEQUENCE_RANKS[term]
                    # Gene ID will be missing for intergenetic consequences
                    gene = consequence.get("gene_id", ".")

                    consequences.append((rank, term, gene, consequence))

                    gene_symbol = consequence.get("gene_symbol")
                    if gene_symbol is not None:
                        gene_symbols.add(gene_symbol)

        if not consequences:
            if allele == vep["allele_string"] or allele.startswith("*"):
                # No consequences for non-variable sites
                return {}

            raise ValueError((allele, vep))

        consequences.sort(key=lambda it: it[0])
        # One of the most significant consequences is picked "randomly"
        _, most_significant, gene_id, consequence = consequences[0]

        consequence["all_gene_symbols"] = ";".join(gene_symbols) or "."

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
            "all_gene_symbols",
            "impact",
            "strand",
            "least_significant",
            "most_significant",
            "n_most_significant",
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

        dst["AncestralAllele"] = allele

    def _add_custom_annotation(self, src, dst, name, fields, default="."):
        data = {}

        alt = dst[":vep:"]["alt"]
        # annotation = {"fields": ..., "name": "chr:start-end"}
        annotations = src.get("custom_annotations", {}).get(name, {})
        for annotation in annotations:
            if annotation.get("allele", alt) == alt:
                # data = {"FILTER": ".", field1: value1, ...}
                data = annotation.get("fields", {})
                break

        if not isinstance(fields, dict):
            fields = dict(zip(fields, fields))

        for src_key, dst_key in fields.items():

            dst[dst_key] = data.get(src_key, default)

    def _add_1k_genomes_annotation(self, src, dst):
        self._add_custom_annotation(
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

    def _add_dbsnp_annotation(self, src, dst):
        self._add_custom_annotation(
            src=src,
            dst=dst,
            name="dbSNP",
            fields={
                "ids": "dbSNP_ids",
                "alts": "dbSNP_alts",
                "functions": "dbSNP_functions",
            },
        )

    def _add_gnomad_annotation(self, src, dst):
        self._add_custom_annotation(
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

        self._add_custom_annotation(
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

    def _add_clinvar_annotation(self, src, dst):
        self._add_custom_annotation(
            src=src,
            dst=dst,
            name="ClinVar",
            fields={
                "ALLELEID": "ClinVar_ID",
                "CLNDN": "ClinVar_Disease",
                "CLNSIG": "ClinVar_Significance",
            },
        )


class Output:
    def __init__(self, keys):
        self.keys = keys


class JSONOutput(Output):
    def print_header(self, out):
        pass

    def print_row(self, out, data):
        json.dump({key: data[key] for key in self.keys}, out)
        out.write("\n")


class TSVOutput(Output):
    def print_header(self, out):
        for name, description in self.keys.items():
            label = f"# {name}: "
            lines = textwrap.wrap(description, width=88 - len(label))

            for line in lines or [""]:
                print(f"{label}{line}", file=out)
                label = "# {}  ".format(" " * len(name))

        print("\t".join(map(str, self.keys)), file=out)

    def print_row(self, out, data):
        print("\t".join(map(str, (data[key] for key in self.keys))), file=out)


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument(
        "in_vcf",
        type=Path,
    )

    parser.add_argument(
        "out_file",
        nargs="?",
        type=Path,
    )

    parser.add_argument(
        "--output-format",
        type=str.lower,
        default="tsv",
        choices=("json", "tsv"),
        help="Output format for aggregated annotations",
    )

    parser.add_argument(
        "--database",
        default="hg38",
        choices=("hg19", "hg38"),
    )

    parser.add_argument(
        "--vep-output",
        required=True,
        type=Path,
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
    log.info("reading variants from '%s'", args.in_vcf)

    with pysam.VariantFile(str(args.in_vcf)) as handle:
        # order of operations
        annotations = [
            AnnotateBasicsInfo(),
            AnnotateLiftOver(
                source=args.database,
                cache=args.liftover_cache,
            ),
            AnnotateVEP(args.vep_output),
        ]

        header = {}
        for annotator in annotations:
            header.update(annotator.keys())

        out_handle = sys.stdout
        if args.out_file:
            out_handle = open(args.out_file, "wt")

        if args.output_format == "json":
            writer = JSONOutput(header)
        else:
            writer = TSVOutput(header)

        try:
            writer.print_header(out_handle)
            for read in handle.fetch():
                rows = [{}]
                for annotator in annotations:
                    annotated_rows = []
                    for row in rows:
                        annotated_rows.extend(annotator.annotate(read, row))

                    rows = annotated_rows

                for row in rows:
                    writer.print_row(out_handle, row)
        except BrokenPipeError:
            # Gracefully handle use of head, tail, less, etc.
            return 0
        finally:
            out_handle.close()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
