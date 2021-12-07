#!/usr/bin/env python3
# -*- coding: utf8 -*-
import logging
import os
import subprocess
import sys

from utils import cmd_to_str, join_procs, update_required


def exec(log, command, stdin=subprocess.DEVNULL, stdout=None):
    log.info("Running %s", cmd_to_str(command))

    return subprocess.Popen(
        command,
        stdin=stdin,
        stdout=stdout,
        close_fds=True,
    )


def exec_self(log, command, stdin=subprocess.DEVNULL, stdout=None):
    import main

    return exec(
        log=log,
        command=[sys.executable, main.__file__] + command,
        stdin=stdin,
        stdout=stdout,
    )


def run_vep(args, log, annotations):
    command = [
        "vep",
        "--verbose",
        "--offline",
        "--cache",
        "--symbol",
        # Marks canonical transcripts
        "--canonical",
        "--format",
        "vcf",
        "--json",
        "--force_overwrite",
        "--compress_output",
        "gzip",
        # Ensure that variants on unknown sequences are still written
        "--dont_skip",
        # Ensure that non-variant sites are still written
        "--allow_non_variant",
        "--dir_cache",
        args.data_cache,
        "--dir_plugins",
        args.install_plugins,
        "--assembly",
        "GRCh38",
        "--polyphen",
        "b",
        "--output_file",
        args.out_vep_json,
        "--stats_file",
        args.out_vep_html,
    ]

    for annotation in annotations:
        command.extend(annotation.params)

    preproc = exec_self(
        log=log,
        command=[
            "pre-process",
            args.in_vcf,
            "--log-level",
            args.log_leve,
        ],
        stdout=subprocess.PIPE,
    )

    vepproc = exec(log=log, command=command, stdin=preproc.stdout)

    return join_procs(log, [preproc, vepproc])


def run_post_proc(args, log):
    command = [
        "post-process",
        args.out_vep_json,
        args.out_prefix,
        "--log-level",
        args.log_leve,
        "--data-liftover",
        args.data_liftover,
    ]

    for annotation in args.annotations:
        command += ["--annotations", annotation]

    proc = exec_self(log=log, command=command)

    return join_procs(log, [proc])


def main(args, annotations):
    log = logging.getLogger("annovep")

    any_errors = False
    for annotation in annotations:
        log.info("Checking files for annotation %s", annotation.name)
        for filename in annotation.files:
            if not os.path.exists(filename):
                log.error("Required %s file %r not found", annotation.name, filename)
                any_errors = True

    if any_errors:
        return 1

    args.out_vep_json = f"{args.out_prefix}.vep.json.gz"
    args.out_vep_html = f"{args.out_prefix}.vep.html"

    if update_required(
        output=args.out_vep_json,
        inputs=[args.in_vcf] + args.annotations,
    ):
        log.info("Running VEP")
        if not run_vep(args, log, annotations):
            return 1
    else:
        log.info("VEP annotations already up to date")

    log.info("Running post-processing")
    if not run_post_proc(args, log):
        return 1

    return 0
