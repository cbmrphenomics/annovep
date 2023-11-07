from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Iterator

from typing_extensions import TypedDict
from utils import open_rb


class VEPRecord(TypedDict):
    Chr: str
    Pos: int
    ID: list[str]
    Ref: str
    Alts: list[str]
    Quality: None | float
    Filters: list[str]
    Info: list[str]
    Samples: list[dict[str, str]]
    VEP: object


class VEPReader:
    def __init__(self, filename: Path):
        self._first_record = None

        self._log = logging.getLogger(__name__)
        self._handle = open_rb(filename)
        self.metadata = self._read_metadata()
        self.timestamp = filename.stat().st_mtime

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

    def _read_record(
        self,
        line: bytes,
        nan_re: re.Pattern[bytes] = re.compile(rb":(-)?NaN\b", flags=re.I),
    ) -> VEPRecord:
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

        samples: list[dict[str, str]] = []
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

    def __iter__(self) -> Iterator[VEPRecord]:
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
            else:
                count += 1

            yield record

        if chrom is not None:
            self._log.info("processed %i records on %r", count, chrom)


def decode_contig_name(name: str) -> str:
    """Decode contig name encoded by `preprocess_vcf.py`"""
    if name.startswith("annovep_"):
        return bytes.fromhex(name[8:]).decode("utf-8")

    return name
