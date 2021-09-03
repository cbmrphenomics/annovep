#!/usr/bin/env python3
# -*- coding: utf8 -*-
import logging
import os
import pwd
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Union

import coloredlogs

OK = "[✓]"
ERR = "[☓]"

MARKER = "not_mounted"

DATA_ROOT = Path("/data")
USER_ROOT = DATA_ROOT / "user"
CACHE_ROOT = DATA_ROOT / "cache"

ANNOVEP_ROOT = Path("/opt/annovep")

VEP_ROOT = Path("/opt/vep/src/ensembl-vep")
VEP_PLUGINS = Path("/opt/vep-plugins/Plugins/")


COMMANDS: Dict[str, List[Union[Path, str]]] = {
    "bash": ["/bin/bash"],
    # Conversion of formats/tables
    "convert:gnomadflags": ["python3", ANNOVEP_ROOT / "convert_gnomad_flags.py"],
    # VEP setup and direct execution
    "vep": ["perl", VEP_ROOT / "vep"],
    "vep:setup": ["bash", ANNOVEP_ROOT / "setup_vep.sh"],
    # Main pipeline
    "pipeline": ["bash", ANNOVEP_ROOT / "pipeline.sh"],
}


def quote_command(args):
    return [shlex.quote(str(value)) for value in args]


def user_exists(uid):
    try:
        pwd.getpwuid(uid)
    except KeyError:
        return False

    return True


def user_add(uid, name):
    return not subprocess.call(["useradd", "-u", str(uid), name])


def check_container(log):
    log.info("Checking container ..")

    folders = [
        ["cache", CACHE_ROOT],
        ["user data", USER_ROOT],
    ]

    for name, root in folders:
        if not root.is_dir():
            log.error(" %s annovep %s folder does not exist: %s", ERR, name, root)
            return False
        log.info("  %s annovep %s folder exists ..", OK, name)

        if (root / MARKER).exists():
            log.error(" %s annovep %s folder not mounted in docker!", ERR, name)
            return False
        log.info("  %s annovep %s folder mounted ..", OK, name)

    return True


def check_databases(log):
    log.info("VEP cache accessible at '%s'", CACHE_ROOT)

    return True


def change_user(log):
    target_uid = os.environ.get("ANNOVEP_USER", "1000")
    try:
        target_uid = int(target_uid)
    except Exception:
        log.info("invalid ANNOVEP_USER: %r", target_uid)
        return False

    if target_uid < 1000:
        # FIXME: Better solution needed
        log.error("cannot change user to %i", target_uid)
        return False

    target_uid = int(target_uid)
    if not user_exists(target_uid):
        log.info("adding dummy user")
        user_add(target_uid, "dummy")

    log.info("changing to gid %s", target_uid)
    os.setgid(target_uid)
    log.info("changing to uid %s", target_uid)
    os.setuid(target_uid)

    return True


def main(argv):
    coloredlogs.install(
        level="INFO",
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(levelname)s %(message)s",
        format="%(asctime)s %(levelname)s %(message)s",
    )

    log = logging.getLogger("annovep")
    if not argv:
        log.error("No command specified. Available commands are ")
        for nth, key in enumerate(COMMANDS, start=1):
            log.error("  %i. %s", nth, key)

        return 1

    command, *argv = argv
    commandline = COMMANDS.get(command.lower())
    if commandline is None:
        log.error("unknown command %r", command)
        return 1

    if not change_user(log):
        return 1

    if not check_container(log):
        if command != "bash":
            return 1

    if not check_databases(log):
        if command != "bash":
            return 1

    # Set home to user data folder
    env = dict(os.environ)

    env["HOME"] = "/home/dummy"

    env["ANNOVEP_CACHE"] = str(CACHE_ROOT)

    env["VEP_ROOT"] = str(VEP_ROOT)
    env["VEP_PLUGINS"] = str(VEP_PLUGINS)

    log.info("Running %s", " ".join(quote_command(commandline + argv)))

    return subprocess.call(
        commandline + argv,
        stdin=sys.stdin,
        cwd=USER_ROOT,
        env=env,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
