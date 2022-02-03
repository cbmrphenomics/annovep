#!/usr/bin/env python3
# -*- coding: utf8 -*-
import logging
from typing import Dict

from . import output
from .annotations import Annotator
from .reader import VEPReader


def main(args, annotations):
    if not any(fmt.startswith("sql") for fmt in args.output_format):
        args.include_json = False

    log = logging.getLogger("convert_vep")
    log.info("reading VEP annotations from '%s'", args.in_file)

    output_formats = set(args.output_format)
    if not output_formats:
        output_formats = ["tsv"]

    if args.out_prefix is None and len(output_formats) > 1:
        log.error("[out_prefix] must be set when writing more than one format")
        return 1

    vep_reader = VEPReader(args.in_file)

    annotations = Annotator(
        annotations=annotations,
        metadata=vep_reader.metadata,
        liftover_cache=args.data_liftover,
    )

    header = output.columns(annotations.fields)
    writers: Dict[str, output.Output] = {}
    for key in output_formats:
        cls = output.FORMATS[key]
        writers[key] = cls(keys=header, out_prefix=args.out_prefix)

    try:
        for record in vep_reader:
            if args.include_json:
                for writer in writers.values():
                    writer.process_json(record)

            for row in annotations.annotate(record):
                for writer in writers.values():
                    writer.process_row(row)

        for writer in writers.values():
            writer.finalize()
    except BrokenPipeError:
        pass

    return 0
