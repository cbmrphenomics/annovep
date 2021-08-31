#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import collections
import gzip
import json
import logging
import re
import sys
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
        row["Ref"] = self._validate_sequence(vcf.ref, "ACGTN*")
        row["DP"] = self._calculate_depth(vcf)

        # There is no easy way to get just the INFO field a str
        vcf_as_text = str(vcf).split("\t", 8)
        # For VCF files with no genotypes, the INFO field may contain trailing newlines
        row["Info"] = vcf_as_text[7].rstrip()

        genotype_counts = self._calculate_genotype_counts(vcf)
        frequencies = self._calculate_allele_freqs(genotype_counts)

        # Construct the cleaned up alleles / positions used by Annovar/VEP
        vep_alleles = self._construct_vep_alleles(vcf)
        annovar_alleles = self._construct_annovar_alleles(vcf)

        for allele in vcf.alts:
            allele = self._validate_sequence(allele, "ACGTN*")

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

            # Cleaned up coordinates/sequences used by annovar and vep
            copy[":annovar:"] = annovar_alleles[allele]
            copy[":vep:"] = vep_alleles[allele]

            yield copy

    def keys(self):
        return (
            "Chr",
            "Pos",
            "Ref",
            "Alt",
            "DP",
            "Freq",
            "GT_00",
            "GT_01",
            "GT_11",
            "GT_NA",
            "GT_Other",
            "Info",
        )

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

    def _construct_annovar_alleles(self, vcf):
        alleles = {}
        for vcf_alt in vcf.alts:
            start, end, ref, alt = self._construct_annovar_variant(
                vcf.pos, vcf.ref, vcf_alt
            )

            alleles[vcf_alt] = {"start": start, "end": end, "ref": ref, "alt": alt}

        return alleles

    def _construct_annovar_variant(self, start, ref, alt):
        # Annovar doens't understand "*"
        if alt == "*":
            alt = "0"

        end = start + len(ref) - 1
        if len(ref) == 1 and len(alt) == 1:  # SNV
            return (start, end, ref, alt)
        elif self._keepindelref:
            return (start, end, ref, alt)
        elif ref.startswith(alt):  # Deletion
            start = start + len(alt)
            (ref, alt) = (ref[len(alt) :], "-")
        elif alt.startswith(ref):  # Insertion
            start = end
            (ref, alt) = ("-", alt[len(ref) :])
        elif ref[:-1] == alt[:-1]:  # Block substitution
            start = end
            (ref, alt) = (ref[-1], alt[-1])

        # 20150324: further adjust when only part of alt and ref matches
        return self._trim_ref_and_alt(start, end, ref, alt)

    def _trim_ref_and_alt(self, start, end, ref, alt):
        while ref and alt and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            end -= 1

        while ref and alt and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            start += 1

        if not ref:
            ref = "-"
            # now it is an insertion so the start should decrease by 1 (20150925)
            start -= 1
        elif not alt:
            alt = "-"

        return (start, end, ref, alt)

    def _construct_vep_alleles(self, vcf):
        start = vcf.pos
        ref = vcf.ref
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

        allele_string = "/".join([ref] + alts)

        return {
            vcf_alt: {
                "start": start,
                "ref": ref,
                "alt": vep_alt,
                "alleles": allele_string,
            }
            for vcf_alt, vep_alt in zip(vcf.alts, alts)
        }


class AnnotateAnnovar(Annotator):
    CACHE_SIZE = 1000

    def __init__(self, filepath: Path) -> None:
        self._log = logging.getLogger("annovar")
        self._log.info("reading Annovar annotations from '%s'", filepath)

        self._handle = filepath.open("rt")
        self._header = self._handle.readline().rstrip().split("\t")

        self._mapping = {
            "gnomAD_ALL_AF": "AF",
            "gnomAD_AFR_AF": "AF_afr",
            "gnomAD_AMI_AF": "AF_ami",
            "gnomAD_AMR_AF": "AF_amr",
            "gnomAD_ASJ_AF": "AF_asj",
            "gnomAD_EAS_AF": "AF_eas",
            "gnomAD_FIN_AF": "AF_fin",
            "gnomAD_NFE_AF": "AF_nfe",
            "gnomAD_OTH_AF": "AF_oth",
            "gnomAD_SAS_AF": "AF_sas",
            "1KG_AFR_AF": "AFR.sites.2015_08",
            "1KG_AMR_AF": "AMR.sites.2015_08",
            "1KG_EAS_AF": "EAS.sites.2015_08",
            "1KG_EUR_AF": "EUR.sites.2015_08",
            "1KG_SAS_AF": "SAS.sites.2015_08",
        }

    def annotate(self, vcf, row):
        data = self._read_record(row)

        for dst, src in self._mapping.items():
            row[dst] = data[src]

        yield row

    def keys(self):
        return self._mapping.keys()

    def _read_record(self, row):
        annovar = row[":annovar:"]
        expected = (
            row["Chr"],
            annovar["start"],
            annovar["end"],
            annovar["ref"],
            annovar["alt"],
        )

        # Annovar will occassionally generate multiple lines for the same allele
        line = self._handle.readline().rstrip()
        row = dict(zip(self._header, line.split("\t")))
        observed = (
            row["Chr"],
            int(row["Start"]),
            int(row["End"]),
            row["Ref"],
            row["Alt"],
        )

        if expected != observed:
            raise ValueError(f"Annovar: {observed} != {expected}")

        return row


class AnnotateLiftOver(Annotator):
    def __init__(self, source: str, cache: Path) -> None:
        if source == "hg19":
            destination = "hd38"
        elif source == "hg38":
            destination = "hg19"
        else:
            raise ValueError(source)

        self._src_db = source.upper()
        self._dst_db = destination.upper()
        self._lifter = liftover.get_lifter(source, destination, cache)

    def annotate(self, vcf, row):
        src_chrom = row["Chr"]
        src_pos = row["Pos"]

        row[f"Chr_{self._src_db}"] = src_chrom
        row[f"Pos_{self._src_db}"] = src_pos

        # Returns list of overlapping liftover coordinates. May be an empty list if
        # the position does not exist in the target genome.
        # FIXME: Are multiple positions possible? And if so, how best to handle?
        coordinates = self._lifter.query(src_chrom, int(src_pos))
        if len(coordinates) == 1:
            ((chrom, pos, _strand),) = coordinates
        else:
            chrom = pos = "."

        row[f"Chr_{self._dst_db}"] = chrom
        row[f"Pos_{self._dst_db}"] = pos

        return [row]

    def keys(self):
        return ("Chr_HG19", "Pos_HG19", "Chr_HG38", "Pos_HG38")


class AnnotateQC(Annotator):
    def __init__(self) -> None:
        self._log = logging.getLogger("qc")
        self._log.warning("qc pass not implemented")

    def annotate(self, vcf, row):
        row["QC"] = "."

        return [row]

    def keys(self):
        return ["QC"]


class AnnotateVEP(Annotator):
    CACHE_SIZE = 1000

    def __init__(self, filepath: Path) -> None:
        self._log = logging.getLogger("vep")
        self._log.info("reading VEP annotations from '%s'", filepath)

        self._handle = gzip.open(filepath, "rt")

        self._cached_for = None
        self._cached_record = {}

    def annotate(self, vcf, row):
        vep = self._read_record(vcf, row)

        consequence = self._get_allele_consequence(vep, row[":vep:"]["alt"])

        self._add_gene_names(consequence, row)
        self._add_ancestral_allele(consequence, row)
        self._add_polyphen_prediction(consequence, row)
        self._add_gnomad_annotation(vep, row)
        self._add_clinvar_annotation(vep, row)

        return [row]

    def keys(self):
        return (
            "AncestralAllele",
            "Gene_id",
            "Gene_symbol",
            "Gene_symbol_source",
            "Polyphen_score",
            "Polyphen_prediction",
            "ClinVar_ID",
            "ClinVar_Disease",
            "ClinVar_Significance",
            # gnomAD are last since annovar gnomADs are appended afterwards
            "gnomAD_mean",
            "gnomAD_median",
            "gnomAD_over_15",
            "gnomAD_over_50",
            "gnomAD_filters",
        )

    def _read_record(self, vcf, row):
        # There is a single VEP associated with each VCF record (provided that the
        # --allow_non_variant option used. Thus we only need to read a new record when
        # a new VCF record is passed.
        if self._cached_for is not vcf:
            line = self._handle.readline()
            # Workaround for non-standard JSON output observed in some records; Python
            # accepts "NaN", but null seems more reasonable
            line = re.sub(r":(-)?nan\b", r":null", line, flags=re.I)

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

        for consequence in transcript_consequences or intergenic_consequences:
            if consequence["variant_allele"] == allele:
                return consequence

        # FIXME: Are starts the only missing items? And is this the best way to deal
        #        with them?
        assert allele.startswith("*"), vep

        return {}

    def _add_gene_names(self, consequence, dst):
        dst["Gene_id"] = consequence.get("gene_id")
        dst["Gene_symbol"] = consequence.get("gene_symbol")
        dst["Gene_symbol_source"] = consequence.get("gene_symbol_source")

    def _add_ancestral_allele(self, consequence, dst):
        allele = consequence.get("aa")
        # The explicty check for falsey values is used to catch both missing values and
        # as a workaround for bug where "aa" is -nan (converted to None in _read_record)
        if not allele:
            allele = "." * len(dst["Ref"])

        dst["AncestralAllele"] = allele

    def _add_polyphen_prediction(self, consequence, dst):
        dst["Polyphen_score"] = consequence.get("polyphen_score", ".")
        dst["Polyphen_prediction"] = consequence.get("polyphen_prediction", ".")

    def _add_custom_annotation(self, src, dst, name, fields, default="."):
        data = {}

        # annotation = {"fields": ..., "name": "chr:start-end"}
        annotations = src.get("custom_annotations", {}).get(name, {})
        for annotation in annotations:
            if annotation.get("allele", dst["Alt"]) == dst["Alt"]:
                # data = {"FILTER": ".", field1: value1, ...}
                data = annotation.get("fields", {})
                break

        if not isinstance(fields, dict):
            fields = dict(zip(fields, fields))

        for src_key, dst_key in fields.items():
            dst[dst_key] = data.get(src_key, default)

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
            name="gnomAD_filters",
            fields=("gnomAD_filters",),
        )

    def _add_clinvar_annotation(self, src, dst):
        # Available annotation (need to be added in pipeline script):
        #   AF_ESP        float  allele frequencies from GO-ESP
        #   AF_EXAC       float  allele frequencies from ExAC
        #   AF_TGP        float  allele frequencies from TGP
        #   ALLELEID      int    the ClinVar Allele ID
        #   CLNDN         str    ClinVar's preferred disease name for the concept
        #                        specified by disease identifiers in CLNDISDB
        #   CLNDNINCL     str    For included Variant : ClinVar's preferred disease name
        #                        for the concept specified by disease identifiers in
        #                        CLNDISDB
        #   CLNDISDB      str    Tag-value pairs of disease database name and
        #                        identifier, e.g. OMIM:NNNNNN
        #   CLNDISDBINCL  str    For included Variant: Tag-value pairs of disease
        #                        database name and identifier, e.g. OMIM:NNNNNN
        #   CLNHGVS       str    Top-level (primary assembly, alt, or patch) HGVS
        #                        expression.
        #   CLNREVSTAT    str    ClinVar review status for the Variation ID
        #   CLNSIG        str    Clinical significance for this single variant
        #   CLNSIGCONF    str    Conflicting clinical significance for this single
        #                        variant
        #   CLNSIGINCL    str    Clinical significance for a haplotype or genotype that
        #                        includes this variant. Reported as pairs of
        #                        VariationID:clinical significance.
        #   CLNVC         str    Variant type
        #   CLNVCSO       str    Sequence Ontology id for variant type
        #   CLNVI         str    the variant's clinical sources reported as tag-value
        #                        pairs of database and variant identifier
        #   DBVARID       str    nsv accessions from dbVar for the variant
        #   GENEINFO      str    Gene(s) for the variant reported as gene symbol:gene
        #                        id. The gene symbol and id are delimited by a colon (:)
        #                        and each pair is delimited by a vertical bar (|)

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


class AnnotateSplitByAllele(Annotator):
    def annotate(self, vcf, row):
        for data in row.pop("alleles").values():
            data.update(row)

            yield data

    def keys(self):
        return ()


def setup_annotators(args, log, vcf) -> List[Annotator]:
    # order of operations
    annotations: Dict[str, Optional[Annotator]] = {
        "basic": AnnotateBasicsInfo(),
        "liftover": None,
        "vep": None,
        # "split": AnnotateSplitByAllele(),
        "annovar": None,
        "qc": AnnotateQC(),
    }

    if args.liftover_cache:
        annotations["liftover"] = AnnotateLiftOver(
            source=args.database,
            cache=args.liftover_cache,
        )

    if args.annovar_output:
        annotations["annovar"] = AnnotateAnnovar(args.annovar_output)

    if args.vep_output:
        annotations["vep"] = AnnotateVEP(args.vep_output)

    selection: List[Annotator] = []
    for key, value in annotations.items():
        if value is None:
            log.warning("%s annotations disabled", key)
        else:
            selection.append(value)

    return selection


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
        "out_tsv",
        nargs="?",
        type=Path,
    )

    parser.add_argument(
        "--database",
        default="hg38",
        choices=("hg19", "hg38"),
    )

    parser.add_argument(
        "--annovar-output",
        type=Path,
    )

    parser.add_argument(
        "--vep-output",
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
        annotations = setup_annotators(args, log, handle)

        header = []
        for annotator in annotations:
            header.extend(annotator.keys())

        out_handle = sys.stdout
        if args.out_tsv:
            out_handle = open(args.out_tsv, "wt")

        try:
            print("\t".join(map(str, header)), file=out_handle)
            for read in handle.fetch():
                rows = [{}]
                for annotator in annotations:
                    annotated_rows = []
                    for row in rows:
                        annotated_rows.extend(annotator.annotate(read, row))

                    rows = annotated_rows

                for row in rows:
                    print(
                        "\t".join(map(str, (row[key] for key in header))),
                        file=out_handle,
                    )
        except BrokenPipeError:
            # Gracefully handle use of head, tail, less, etc.
            return 0
        finally:
            out_handle.close()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
