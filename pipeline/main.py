#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import sys
from pathlib import Path

import coloredlogs

from annotation import load_annotations
from pipeline import main as pipeline_main
from postprocess import main as postprocess_main
from preprocess import main as preprocess_main


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def add_common_args(parser):
    parser.add_argument("--root", type=Path, default=Path("~/annovep").expanduser())
    parser.add_argument(
        "--log-level",
        default="info",
        choices=("debug", "info", "warning", "error"),
        type=str.lower,
        help="Log messages at the specified level. This option applies to the "
        "`--log-file` option and to log messages printed to the terminal.",
    )


def parse_args(argv):
    mainparser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    mainparser.set_defaults(
        main=None,
        root=None,
        annotations=[],
        data_cache=None,
        data_custom=None,
        data_plugins=None,
        data_liftover=None,
        install=None,
        install_annovep=None,
        install_plugins=None,
    )

    subparsers = mainparser.add_subparsers()

    # Pipeline
    parser = subparsers.add_parser("run")
    parser.set_defaults(main=pipeline_main)

    parser.add_argument("in_vcf", type=Path)
    parser.add_argument("out_prefix", type=Path)

    parser.add_argument("--annotations", default=[], action="append", type=Path)

    group = parser.add_argument_group("Data locations")
    parser.add_argument("--data-cache", type=Path)
    parser.add_argument("--data-custom", type=Path)
    parser.add_argument("--data-plugins", type=Path)
    parser.add_argument("--data-liftover", type=Path)

    group = parser.add_argument_group("Installation locations")
    group.add_argument("--install", type=Path)
    group.add_argument("--install-plugins", type=Path)
    group.add_argument("--install-annovep", type=Path)

    group = parser.add_argument_group("VEP options")
    group.add_argument("--fork", type=int, default=1)

    add_common_args(parser)

    # Pre-processing
    parser = subparsers.add_parser("pre-process")
    parser.set_defaults(main=preprocess_main)

    parser.add_argument("in_vcf", type=Path)

    add_common_args(parser)

    # Post-processing
    parser = subparsers.add_parser("post-process")
    parser.set_defaults(main=postprocess_main)

    parser.add_argument("in_json", type=Path)
    parser.add_argument("out_prefix", type=Path)

    parser.add_argument("--annotations", default=[], action="append", type=Path)
    parser.add_argument("--data-liftover", type=Path)

    parser.add_argument(
        "--output-format",
        action="append",
        type=str.lower,
        choices=("tsv", "json", "sql"),
        help="Output format for aggregated annotations. Maybe be specified zero or "
        "more times. Defaults to TSV if not specified",
    )

    parser.add_argument(
        "--include-json",
        action="store_true",
        help="Include JSON data in SQL output, excluding sample specific information",
    )

    add_common_args(parser)

    return mainparser


def main(argv):
    coloredlogs.install(
        level="INFO",
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
        # Workaround for coloredlogs disabling colors in docker containers
        isatty=sys.stderr.isatty(),
    )

    parser = parse_args(argv)
    args = parser.parse_args(argv)

    if args.data_cache is None:
        args.data_cache = args.root / "cache"
    if args.data_custom is None:
        args.data_custom = args.root / "custom"
    if args.data_plugins is None:
        args.data_plugins = args.root / "plugins"
    if args.data_liftover is None:
        args.data_liftover = args.root / "liftover"

    if args.install is None:
        args.install = args.root / "install"
    if args.install_plugins is None:
        args.install_plugins = args.install / "vep-plugins" / "Plugins"

    variables = {
        # Data folders
        "data-cache": args.data_cache,
        "data-custom": args.data_custom,
        "data-plugins": args.data_plugins,
        "data-liftover": args.data_liftover,
        # Installation folders
        "install": args.install,
        "install-plugins": args.install_plugins,
        "install-annovep": args.install_annovep,
    }

    annotations = tuple(load_annotations(args.annotations, variables))

    if args.main is None:
        parser.print_usage()
        return 1

    return args.main(args, annotations)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
