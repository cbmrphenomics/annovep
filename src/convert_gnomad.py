#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import bz2
import gzip
import io
import logging
import sys
from os import fspath
from pathlib import Path
from typing import IO, Union, cast

import coloredlogs

TEMPLATE = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n"


class Counter:
    def __init__(self, log):
        self._log = log
        self._chrom = None
        self._count = 0

    def __call__(self, chrom):
        if chrom != self._chrom:
            self.print()
            self._chrom = chrom
            self._count = 1
        else:
            self._count += 1

    def print(self):
        if self._chrom is not None:
            self._log.info("Processed %i sites on %s", self._count, self._chrom)

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, *args, **kwargs):
        self.print()


def filters_to_vcf(counter, in_handle, out_handle):
    for line in in_handle:
        site_info, filters = line.split()
        chrom, start, ref, alt = site_info.split(":")

        out_handle.write(
            TEMPLATE.format(
                chrom=chrom,
                pos=start,
                id=".",
                ref=ref,
                alt=alt,
                qual=".",
                filter=".",
                info=f"gnomAD_filters={filters}",
            )
        )

        counter(chrom)


def coverage_to_vcf(counter, in_handle, out_handle):
    header = in_handle.readline().rstrip().split("\t")

    for line in in_handle:
        fields = line.rstrip().split("\t")
        row = dict(zip(header, fields))

        info = [
            "gnomAD_mean={}".format(row["mean"]),
            "gnomAD_median={}".format(row["median_approx"]),
            "gnomAD_over_15={}".format(row["over_15"]),
            "gnomAD_over_50={}".format(row["over_50"]),
        ]

        out_handle.write(
            TEMPLATE.format(
                chrom=row["#chr"],
                pos=row["pos"],
                id=".",
                ref=".",
                alt=".",
                qual=".",
                filter=".",
                info=";".join(info),
            )
        )

        counter(row["#chr"])


def open_ro(filename: Union[str, Path]) -> IO[str]:
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(fspath(filename), "rb")
    try:
        header = handle.read(2)
        handle.seek(0)

        if header == b"\x1f\x8b":
            handle = cast(IO[bytes], gzip.GzipFile(mode="rb", fileobj=handle))
        elif header == b"BZ":
            handle = bz2.BZ2File(handle, "rb")

        return io.TextIOWrapper(handle)
    except:
        handle.close()
        raise


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("format", choices=("coverage", "filters"), type=str.lower)
    parser.add_argument("in_txt", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    coloredlogs.install(
        level="INFO",
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(levelname)s %(message)s",
        format="%(asctime)s %(levelname)s %(message)s",
    )

    log = logging.getLogger("__main__")
    with open_ro(args.in_txt) as in_handle:
        with Counter(log) as counter:
            if args.format == "filters":
                log.info("Converting gnomAD filters to VCF")
                return filters_to_vcf(counter, in_handle, sys.stdout)
            elif args.format == "coverage":
                log.info("Converting gnomAD coverage to VCF")
                return coverage_to_vcf(counter, in_handle, sys.stdout)
            else:
                raise NotImplementedError(args.format)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
