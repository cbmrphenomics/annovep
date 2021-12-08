from typing import Dict, Optional, NamedTuple

import ruamel.yaml


class AnnotationError(Exception):
    pass


class Field(NamedTuple):
    name: str
    type: str
    help: str
    split_by: Optional[str]


def _str_list(value):
    """Require that a value is List[str]"""
    if not isinstance(value, list):
        raise AnnotationError("not a list")
    elif any(not isinstance(it, str) for it in value):
        raise AnnotationError("non-string values found")

    return value


def _str_dict(value):
    """Require that a value is Dict[str, str]"""
    if not isinstance(value, dict):
        raise AnnotationError("not a dict")
    elif not all(isinstance(it, str) for it in value):
        raise AnnotationError("non-string keys found")
    elif not all(isinstance(it, str) for it in value.values()):
        raise AnnotationError("non-string values found")

    return value


def _truthy(type):
    """Wraper around _str_list/_str_dict to also require a truthy (non-empty) value"""

    def _is_truthy(value):
        value = type(value)
        if not value:
            raise AnnotationError("no value")

        return value

    return _is_truthy


def _parse(layout, name, data):
    if not isinstance(name, str):
        raise AnnotationError(f"Name {name!r} is not a string")
    elif not isinstance(data, dict):
        raise AnnotationError(f"{name!r} settings are not a dict")

    output = {}
    for key, info in layout.items():
        value = data.pop(key, None)
        if value is None:
            if "default" not in info:
                raise AnnotationError(f"no {key} specified for {name}")

            # No validation of default value to allow e.g. None
            output[key] = info["default"]
            continue

        vtype = info["type"]
        if type(vtype) is type and not isinstance(value, vtype):
            raise AnnotationError(f"{key} for {name} is not a {vtype.__name__}")

        try:
            output[key] = vtype(value)
        except AnnotationError as error:
            raise AnnotationError(f"invalid {key} for {name}: {error}")

    if data:
        raise AnnotationError(f"unexpected settings in {name!r}: {data!r}")

    return output


def combine_variables(variables, user_variables):
    """Combine user Variables with built-in variables"""
    user_variables = dict(user_variables)
    user_variables.update(variables)

    return {
        key: str(value).format(**user_variables)
        for key, value in user_variables.items()
    }


def apply_variables(variables, values):
    for value in values:
        yield value.format(**variables)


def _parse_fields(data, name) -> Dict[str, Field]:
    default_type = data.pop("FieldType")
    data = data.pop("Fields")

    if data is None:
        return {}
    elif not isinstance(data, dict):
        raise AnnotationError(f"Fields for plugin {name!r} are not a dict")

    layout = {
        "Name": {"type": str, "default": None},
        "Help": {"type": str, "default": ""},
        "FieldType": {"type": str, "default": default_type},
        "Split-by": {"type": str, "default": None},
    }

    for key, value in data.items():
        if not isinstance(key, str):
            raise AnnotationError(f"Fields for plugin {name!r} has non-str key {key!r}")
        elif value is None:
            value = {"Name": None}

        value = _parse(layout, f"Field {key} for {name}", value)
        field_type = value["FieldType"]
        if field_type not in ("str", "int", "float"):
            raise AnnotationError(f"Bad type {field_type!r} in {key!r} for {name!r}")

        data[key] = Field(
            name=value["Name"],
            help=value["Help"],
            type=value["FieldType"],
            split_by=value["Split-by"],
        )

    return data


_PLUGIN_LAYOUT = {
    "Rank": {"type": int, "default": 0},
    "Files": {"type": _truthy(_str_list)},
    "Parameters": {"type": _truthy(_str_list)},
    "Variables": {"type": _str_dict, "default": {}},
    "FieldType": {"type": str, "default": "str"},
    "Fields": {"type": lambda it: it, "default": {}},
}


class Plugin:
    def __init__(self, name, data, variables):
        data = _parse(_PLUGIN_LAYOUT, name, data)

        variables = combine_variables(variables, data.pop("Variables"))

        self.rank = data.pop("Rank")
        self.name = name
        self.files = apply_variables(variables, data.pop("Files"))
        self._params = apply_variables(variables, data.pop("Parameters"))
        self.fields = _parse_fields(data, name)

        assert not data, data

    @property
    def params(self):
        params = [self.name]
        params.extend(str(value) for value in self._params)

        return ["--plugin", ",".join(params)]


_CUSTOM_LAYOUT = {
    "Rank": {"type": int, "default": 0},
    "File": {"type": str},
    "Variables": {"type": _str_dict, "default": {}},
    "FieldType": {"type": str, "default": "str"},
    "Fields": {"type": lambda it: it, "default": {}},
}


class Custom:
    def __init__(self, name, data, variables, type):
        if type not in ("vcf", "bed"):
            raise AnnotationError(f"invalid custom annotation type {type!r}")

        data = _parse(_CUSTOM_LAYOUT, name, data)
        variables = combine_variables(variables, data.pop("Variables"))

        self.rank = data.pop("Rank")
        self.name = name
        self._type = type
        self._file = data.pop("File").format(**variables)
        self.fields = _parse_fields(data, name)

        assert not data, data

    @property
    def files(self):
        return [self._file, f"{self._file}.tbi"]

    @property
    def params(self):
        params = [self._file, self.name, self._type, "overlap", "0"]
        for name in self.fields:
            if not (name.startswith(":") and name.endswith(":")):
                params.append(name)

        return ["--custom", ",".join(params)]


def load_annotations(filepaths, variables=None):
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)

    for filepath in filepaths:
        with filepath.open("rt") as handle:
            data = yaml.load(handle)

        for idx, (name, settings) in enumerate(data.items()):
            if not isinstance(settings, dict):
                raise AnnotationError(f"{name} is not a dict")

            type = settings.pop("Type")
            if type == "Plugin":
                value = Plugin(name, settings, variables)
            elif type == "VCF":
                value = Custom(name, settings, variables, type="vcf")
            elif type == "BED":
                value = Custom(name, settings, variables, type="bed")
            else:
                raise AnnotationError(f"Unknown annotation type {type!r} for {name!r}")

            # Ensure that ordering is stable in relation to annotation files
            value.rank = (value.rank, idx)

            yield value
