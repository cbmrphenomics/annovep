#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import logging
import sys
from pathlib import Path

import coloredlogs

from annotation import AnnotationError, load_annotations
from pipeline import main as pipeline_main
from postprocess import main as postprocess_main
from preprocess import main as preprocess_main


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)

    # Pipeline
    parser.add_argument("in_file", type=Path)
    parser.add_argument("out_prefix", type=Path, nargs="?")

    parser.add_argument(
        "--annotations",
        metavar="FILE",
        type=Path,
        default=[],
        action="append",
        help="Optional files containing additional annotations",
    )

    parser.add_argument(
        "--do",
        type=str.lower,
        choices=("run", "pre-process", "post-process"),
        default="run",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--root",
        metavar="DIR",
        type=Path,
        default=Path("~/annovep").expanduser(),
        help="The root location of the AnnoVEP install",
    )

    group = parser.add_argument_group("Data locations")
    parser.add_argument("--data-cache", metavar="DIR", type=Path)
    parser.add_argument("--data-custom", metavar="DIR", type=Path)
    parser.add_argument("--data-plugins", metavar="DIR", type=Path)
    parser.add_argument("--data-liftover", metavar="DIR", type=Path)

    group = parser.add_argument_group("Installation locations")
    group.add_argument("--install", metavar="DIR", type=Path)
    group.add_argument("--install-plugins", metavar="DIR", type=Path)
    group.add_argument("--install-annovep", metavar="DIR", type=Path)

    group = parser.add_argument_group("Output")
    group.add_argument(
        "--output-format",
        default=[],
        action="append",
        type=str.lower,
        choices=("tsv", "json", "sql", "sqlite3"),
        help="Output format for aggregated annotations. Maybe be specified zero or "
        "more times. Defaults to TSV if not specified",
    )

    group.add_argument(
        "--include-json",
        action="store_true",
        help="Include JSON data in SQL output, excluding sample specific information",
    )

    group = parser.add_argument_group("VEP options")
    group.add_argument("--fork", metavar="N", type=int)
    group.add_argument("--buffer_size", metavar="N", default=100_000, type=int)

    group = parser.add_argument_group("Logging")
    group.add_argument(
        "--log-level",
        default="info",
        choices=("debug", "info", "warning", "error"),
        type=str.lower,
        help="Log messages at the specified level. This option applies to the "
        "`--log-file` option and to log messages printed to the terminal.",
    )

    return parser


def main(argv):
    parser = parse_args(argv)
    args = parser.parse_args(argv)

    coloredlogs.install(
        level=args.log_level,
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
        # Workaround for coloredlogs disabling colors in docker containers
        isatty=sys.stderr.isatty(),
    )

    if args.out_prefix is None:
        args.out_prefix = args.in_file

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

    log = logging.getLogger("annovep")
    try:
        annotations = load_annotations(log, args.annotations, variables)
    except AnnotationError as error:
        log.error("error while loading annotations: %s", error)
        return 1

    if args.do == "run":
        return pipeline_main(args, annotations)
    elif args.do == "pre-process":
        return preprocess_main(args, annotations)
    elif args.do == "post-process":
        return postprocess_main(args, annotations)

    raise NotImplementedError(args.do)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
