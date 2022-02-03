import bz2
import gzip
import io
import os
import shlex
import signal
import subprocess
import time
from os import fspath


def cmd_to_str(command, max_length=float("inf")):
    text = " ".join(shlex.quote(str(value)) for value in command)
    if len(text) > max(0, max_length - 3):
        return text[: max_length - 3] + "..."

    return text


def update_required(output, inputs):
    if not os.path.exists(output):
        return True

    mtime = os.stat(output).st_mtime
    for dependency in inputs:
        if os.stat(dependency).st_mtime > mtime:
            return True

    return False


def open_rb(filename):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(fspath(filename), "rb")
    try:
        header = handle.read(2)
        handle.seek(0)

        if header == b"\x1f\x8b":
            handle = gzip.GzipFile(mode="rb", fileobj=handle)
        elif header == b"BZ":
            handle = bz2.BZ2File(handle, "rb")

        return handle
    except:
        handle.close()
        raise


def open_ro(filename):
    return io.TextIOWrapper(open_rb(filename))


def join_procs(log, procs):
    start_time = time.time()
    sleep_time = 0.05
    commands = list(enumerate(procs))
    return_codes = [None] * len(commands)

    assert all(hasattr(cmd, "args") for (_, cmd) in commands)

    while commands and not any(return_codes):
        try:
            # Wait for arbitrary command
            commands[0][1].wait(sleep_time if len(commands) > 1 else None)
        except subprocess.TimeoutExpired:
            sleep_time = min(1, sleep_time * 2)

        for (index, command) in list(commands):
            if command.poll() is not None:
                return_code = command.wait()
                return_codes[index] = return_code
                commands.remove((index, command))
                sleep_time = 0.05

                if return_code < 0:
                    return_code = signal.Signals(-return_code).name

                log.info("Command finished: %s", cmd_to_str(command.args, 80))
                log.info("  Runtime:        %.1fs", time.time() - start_time)
                log.info("  Return-code:    %s", return_code)

    if any(return_codes):
        for index, command in commands:
            log.warning("Terminating command: %s", cmd_to_str(command.args, 80))
            command.terminate()
            return_codes[index] = command.wait()

        log.error("Errors occured during processing!")

    return not any(return_codes)
