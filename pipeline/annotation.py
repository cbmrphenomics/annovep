import ruamel.yaml


def _pop_str_list(data, name, key):
    values = data.pop(key, [])
    if not values:
        raise AnnotationError(f"No {key} for plugin {name!r}")
    elif not isinstance(values, list):
        raise AnnotationError(f"{key} for plugin {name!r} are not a list")

    return tuple(values)


def _replace_vars(values, variables):
    result = []
    for value in values:
        result.append(value.format(**variables))

    return tuple(result)


class AnnotationError(Exception):
    pass


class Plugin:
    def __init__(self, name, data, variables):
        if not isinstance(name, str):
            raise AnnotationError(f"Plugin name {name!r} is not a string")
        elif data is None:
            raise AnnotationError(f"No settings for plugin {name!r}")
        elif not isinstance(data, dict):
            raise AnnotationError(f"Plugin {name!r} settings are not a dict")

        user_variables = data.pop("Variables", {})
        if not isinstance(user_variables, dict):
            raise AnnotationError(f"Variables for plugin {name!r} not a dict")
        variables = dict(variables)
        for key, value in user_variables.items():
            variables[key] = value.format(**variables)

        self.name = name
        self.files = _replace_vars(_pop_str_list(data, name, "Files"), variables)
        self._params = _replace_vars(_pop_str_list(data, name, "Parameters"), variables)

        if data:
            raise AnnotationError(f"Unexpected settings in plugin {name!r}: {data!r}")

    @property
    def params(self):
        params = [self.name]
        params.extend(str(value) for value in self._params)

        return ["--plugin", ",".join(params)]


class Custom:
    def __init__(self, name, data, variables, type):
        if not isinstance(name, str):
            raise AnnotationError(f"Plugin name {name!r} is not a string")
        elif data is None:
            raise AnnotationError(f"No settings for plugin {name!r}")
        elif not isinstance(data, dict):
            raise AnnotationError(f"Plugin {name!r} settings are not a dict")
        elif type not in ("vcf", "bed"):
            raise AnnotationError(f"invalid custom annotation type {type!r}")

        file_ = data.pop("File", None)
        if not file_:
            raise AnnotationError(f"No File specified for custom {name!r}")
        elif not isinstance(file_, str):
            raise AnnotationError(f"{name!r} File is not a string: {self._file!r}")

        if type == "vcf":
            self.fields = _pop_str_list(data, name, "Fields")
        else:
            self.fields = ()

        self.name = name
        self._type = type
        self._file = file_.format(**variables)

        if data:
            raise AnnotationError(f"Unexpected settings in plugin {name!r}: {data!r}")

    @property
    def files(self):
        return [self._file, f"{self._file}.tbi"]

    @property
    def params(self):
        params = [self._file, self.name, self._type, "overlap", "0"]
        params.extend(self.fields)

        return ["--custom", ",".join(params)]


def parse_annotation(cls, data, variables, **kwargs):
    if not isinstance(data, dict):
        raise AnnotationError("Plugins is not a dict")

    for name, settings in data.items():
        yield cls(name, settings, variables, **kwargs)


def load_annotations(filepath, variables):
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)

    with filepath.open("rt") as handle:
        data = yaml.load(handle)

    for key, annotations in data.items():
        if key == "Plugins":
            yield from parse_annotation(Plugin, annotations, variables)
        elif key == "VCFs":
            yield from parse_annotation(Custom, annotations, variables, type="vcf")
        elif key == "BEDs":
            yield from parse_annotation(Custom, annotations, variables, type="bed")
        else:
            raise AnnotationError("Unexpected annotation type {key!r}")
