#!/usr/bin/env python3
# -*- coding: utf8 -*-

try:
    from annovep._version import VERSION
except ModuleNotFoundError:
    VERSION = "Undefined"

__all__ = [
    "VERSION",
]
