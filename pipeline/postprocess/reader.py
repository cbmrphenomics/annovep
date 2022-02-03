import json
import logging
import re

from utils import open_rb


class VEPReader:
    def __init__(self, filename):
        self._log = logging.getLogger(__name__)
        self._handle = open_rb(filename)
        self.metadata = self._read_metadata()
        self._first_record = None

    def _read_metadata(self):
        metadata = {}
        for line in self._handle:
            record = self._read_record(line)
            record_id = ";".join(record["ID"])

            if record_id == "AnnoVEP:Samples":
                metadata["samples"] = record["Info"]
            elif record_id.startswith("AnnoVEP:"):
                self._log.warning("unexpected metadata %r", record_id)
            else:
                self._first_record = record
                break

        return metadata

    def _read_record(self, line, nan_re=re.compile(rb":(-)?NaN\b", flags=re.I)):
        # Workaround for non-standard JSON output observed in some records, where
        # an expected string value was -nan. Python accepts "NaN", but null seems
        # more reasonable for downstream compatibility
        line = nan_re.sub(b":null", line)
        data = json.loads(line)

        vcf_record = data["input"]
        fields = vcf_record.rstrip("\r\n").split("\t")
        chr, pos, id, ref, alt, qual, filters, info, *fmt_and_samples = fields
        chr = decode_contig_name(chr)

        data["input"] = "\t".join((chr, pos, id, ref, alt, qual, filters, info))

        samples = []
        if fmt_and_samples:
            fmt_keys = fmt_and_samples[0].split(":")
            for sample in fmt_and_samples[1:]:
                samples.append(dict(zip(fmt_keys, sample.split(":"))))

        return {
            "Chr": chr,
            "Pos": int(pos),
            "ID": [] if id == "." else id.split(";"),
            "Ref": ref,
            # . is treated as a actual value, rather than an empty list. This is done so
            # that (limited) information can be retrieved for non-specific variants.
            "Alts": alt.split(","),
            "Quality": None if qual == "." else float(qual),
            "Filters": [] if filters == "." else filters.split(";"),
            "Info": [] if info == "." else info.split(";"),
            "Samples": samples,
            "VEP": data,
        }

    def __iter__(self):
        chrom = None
        count = 0

        if self._first_record is not None:
            yield self._first_record
            chrom = self._first_record["Chr"]
            count = 1

            self._first_record = None

        for line in self._handle:
            record = self._read_record(line)
            if record["Chr"] != chrom:
                if chrom is not None:
                    self._log.info("processed %i records on %r", count, chrom)
                chrom = record["Chr"]
                count = 1

            yield record

        if chrom is not None:
            self._log.info("processed %i records on %r", count, chrom)


def decode_contig_name(name):
    """Decode contig name encoded by `preprocess_vcf.py`"""
    if name.startswith("annovep_"):
        return bytes.fromhex(name[8:]).decode("utf-8")

    return name
