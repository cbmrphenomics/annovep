import collections
import json
import sqlite3
import sys
import zlib

from . import consequences


class IntegerCol(str):
    """Indicates that a column contains integer values."""


class FloatCol(str):
    """Indicates that a column contains floating point values."""


COL_TYPES = {
    "str": lambda it: it,
    "int": IntegerCol,
    "float": FloatCol,
}


def columns(annotations):
    columns = collections.OrderedDict()
    for field in annotations.values():
        wrapper = COL_TYPES[field.type]
        columns[field.name] = wrapper(field.help)

    return columns


# Columns that contain consequence terms (see `consequences.ranks()`)
CONSEQUENCE_COLUMNS = (
    "Func_most_significant",
    "Func_least_significant",
    "Func_most_significant_canonical",
)


class Output:
    def __init__(self, keys, out_prefix, extension):
        self.keys = keys
        self._handle = sys.stdout
        if out_prefix is not None:
            self._handle = open(f"{out_prefix}{extension}", "wt")

    def finalize(self):
        if self._handle is not sys.stdout:
            self._handle.close()

    def process_json(self, data):
        pass

    def process_row(self, data):
        raise NotImplementedError()

    def _print(self, line="", *args, **kwargs):
        if args:
            line = line.format(*args)

        print(line, file=self._handle, **kwargs)


class JSONOutput(Output):
    def __init__(self, keys, out_prefix):
        super().__init__(keys, out_prefix, ".json")

    def process_row(self, data):
        json.dump({key: data[key] for key in self.keys}, self._handle)
        self._handle.write("\n")


class TSVOutput(Output):
    def __init__(self, keys, out_prefix):
        super().__init__(keys, out_prefix, ".tsv")

        self._print("#{}", "\t".join(self.keys))

        if out_prefix is not None:
            with open(f"{out_prefix}.tsv.columns", "wt") as handle:
                print("Name\tDescription", file=handle)

                for name, description in self.keys.items():
                    print(name, description, sep="\t", file=handle)

    def process_row(self, data):
        row = [self._to_string(data[key]) for key in self.keys]

        self._print("\t".join(row))

    @staticmethod
    def _to_string(value):
        if isinstance(value, (tuple, list)):
            return ";".join(map(str, value or "."))
        elif value is None:
            return "."

        return str(value)


class SQLOutput(Output):
    # Columns that are renamed for the DB/shiny interface
    COLUMN_MAPPING = {
        "Chr": "Hg38_chr",
        "Pos": "Hg38_pos",
    }

    def __init__(self, keys, out_prefix):
        super().__init__(keys, out_prefix, ".sql")

        self._consequence_ranks = self._build_consequence_ranks()
        self._contigs = {
            "hg19": collections.defaultdict(int),
            "hg38": collections.defaultdict(int),
        }
        self._genes = {}
        self._n_overlap = 0
        self._n_row = 0
        self._n_json = 0

        self._print("PRAGMA TEMP_STORE=MEMORY;")
        self._print("PRAGMA JOURNAL_MODE=OFF;")
        self._print("PRAGMA SYNCHRONOUS=OFF;")
        self._print("PRAGMA LOCKING_MODE=EXCLUSIVE;")

        self._print("BEGIN;")
        self._print()
        self._print_descriptions()
        self._print()
        self._print_consequence_terms()
        self._print()
        self._print_gene_tables()
        self._print()

        for table in ("Annotations",):
            self._print("DROP TABLE IF EXISTS [{}];", table)
        self._print()

        query = [
            "CREATE TABLE [Annotations] (\n",
            "    [pk] INTEGER PRIMARY KEY ASC",
        ]

        for key, description in self.keys.items():
            datatype = "TEXT"
            if key in CONSEQUENCE_COLUMNS:
                key = f"{key}_id"
                datatype = "INTEGER REFERENCES [Consequenes]([pk])"

            # Rename columns for SQL output only
            key = self.COLUMN_MAPPING.get(key, key)

            if isinstance(description, IntegerCol):
                datatype = "INTEGER"
            elif isinstance(description, FloatCol):
                datatype = "REAL"

            query.append(f",\n    [{key}] {datatype}")
        query.append("\n);")
        self._print("".join(query))

        self._print()
        self._print_json_table()
        self._print()
        self._print("END;")
        self._print()
        self._print("BEGIN;")

    def finalize(self):
        self._print("END;")
        self._print()
        self._print("BEGIN;")

        for idx, (gene, info) in enumerate(sorted(self._genes.items())):
            self._print(
                "INSERT INTO [Genes] VALUES ({}, {}, {}, {}, {}, {}, {}, {});",
                idx,
                self._to_string(gene),
                self._to_string(info["Chr"]),
                self._to_string(info["MinPos"]),
                self._to_string(info["MaxPos"]),
                self._to_string(info["Variants"]),
                self._to_string(info["Most_significant"]),
                self._to_string(info["Most_significant_canonical"]),
            )

        self._print()
        self._print_contig_names()
        self._print()
        self._print("CREATE INDEX IPositions_hg38 ON Annotations (Hg38_chr, Hg38_pos);")
        self._print("CREATE INDEX IPositions_hg19 ON Annotations (Hg19_chr, Hg19_pos);")
        self._print("CREATE INDEX IPositions_json ON JSON (Hg38_chr, Hg38_pos);")
        self._print("END;")
        self._print()

        # Collect information for the query planner
        self._print("ANALYZE;")

        self._handle.close()

    def process_json(self, data):
        data = dict(data)

        # Remove any sample specific information and leave only summary information
        vcf_fields = data.pop("input").split("\t", 8)[:8]
        data["input"] = "\t".join(vcf_fields).rstrip("\r\n")

        # Keys are sorted to improve compression ratio
        blob = zlib.compress(json.dumps(data, sort_keys=True).encode("utf-8")).hex()

        self._print(
            "INSERT INTO [JSON] VALUES ({}, {}, {}, {});",
            self._n_json,
            self._to_string(vcf_fields[0]),
            self._to_string(int(vcf_fields[1])),
            "X'{}'".format(blob),
        )

        self._n_json += 1

    def process_row(self, data):
        self._n_row += 1
        self._contigs["hg38"][data["Chr"]] += 1
        self._contigs["hg19"][data["Hg19_chr"]] += 1

        data = dict(data)

        # VEP consequence terms
        for key in CONSEQUENCE_COLUMNS:
            value = data.get(key)
            if value is not None:
                value = self._consequence_ranks[value]

            data[key] = value

        values = [str(self._n_row)]
        values.extend(self._to_string(data[key]) for key in self.keys)

        self._print("INSERT INTO [Annotations] VALUES ({});", ", ".join(values))

        for gene in data.get("Genes_overlapping", ()):
            gene_info = self._genes.get(gene)
            if gene_info is None:
                self._genes[gene] = {
                    "Chr": data["Chr"],
                    "MinPos": data["Pos"],
                    "MaxPos": data["Pos"],
                    "Variants": 1,
                    "Most_significant": data["Func_most_significant"],
                    "Most_significant_canonical": data[
                        "Func_most_significant_canonical"
                    ],
                }
            elif gene_info["Chr"] != data["Chr"]:
                raise ValueError(f"gene {gene!r} found on multiple contigs")
            else:
                gene_info["MinPos"] = min(data["Pos"], gene_info["MinPos"])
                gene_info["MaxPos"] = max(data["Pos"], gene_info["MaxPos"])
                gene_info["Variants"] += 1

                for key in ("Most_significant", "Most_significant_canonical"):
                    gene_info[key] = self._worst_consequence(
                        gene_info[key], data[f"Func_{key.lower()}"]
                    )

    @staticmethod
    def _worst_consequence(consequence_a, consequence_b):
        if consequence_a is None:
            return consequence_b
        elif consequence_b is None:
            return consequence_a

        return max(consequence_a, consequence_b)

    def _print_descriptions(self):
        self._print("DROP TABLE IF EXISTS [Columns];")
        self._print(
            """
            CREATE TABLE [Columns] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Name] TEXT,
              [Table] TEXT,
              [Column] TEXT,
              [Description] TEXT
            );
            """
        )

        for pk, (key, description) in enumerate(self.keys.items()):
            # Rename columns for SQL output only
            key = self.COLUMN_MAPPING.get(key, key)

            table = "Annotations"
            column = key

            if key in CONSEQUENCE_COLUMNS:
                table = "Consequences"
                column = f"{key}_id"

            self._print(
                "INSERT INTO [Columns] VALUES ({}, {}, {}, {}, {});",
                pk,
                self._to_string(key),
                self._to_string(table),
                self._to_string(column),
                self._to_string(description),
            )

    def _print_consequence_terms(self):
        self._print("DROP TABLE IF EXISTS [Consequences];")
        self._print(
            """
            CREATE TABLE [Consequences] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Name] TEXT
            );
            """
        )

        for name, pk in self._consequence_ranks.items():
            self._print(
                "INSERT INTO [Consequences] VALUES ({}, {});",
                pk,
                self._to_string(name),
            )

    def _print_gene_tables(self):
        self._print("DROP TABLE IF EXISTS [Genes];")
        self._print(
            """
            CREATE TABLE [Genes] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Name] TEXT,
              [Hg38_chr] TEXT,
              [Hg38_start] INTEGER,
              [Hg38_end] INTEGER,
              [Variants] INTEGER,
              [Most_significant] INTEGER REFERENCES [Consequenes]([pk]),
              [Most_significant_canonical] INTEGER REFERENCES [Consequenes]([pk])
            );
            """
        )

    def _print_json_table(self):
        self._print("DROP TABLE IF EXISTS [JSON];")
        self._print(
            """
            CREATE TABLE [JSON] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Hg38_chr] TEXT,
              [Hg38_pos] INTEGER,
              [Data] BINARY
            );
            """
        )
        self._print()

    def _print_contig_names(self):
        self._print("DROP TABLE IF EXISTS [Contigs];")
        self._print(
            """
            CREATE TABLE [Contigs] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Build] TEXT,
              [Name] TEXT,
              [Variants] INTEGER
            );
            """
        )

        contigs = []
        overlap = self._contigs["hg19"].keys() & self._contigs["hg38"].keys()

        # Collect hg38 contigs. These should be in the proper order
        for name, variants in self._contigs["hg38"].items():
            contigs.append(("hg38", name, variants))
            if name in overlap:
                contigs.append(("hg19", name, variants))

        # Collect hg19 only contigs; unmapped variants are ignored
        for name in sorted(self._contigs["hg19"].keys() - overlap - set([None])):
            contigs.append(("hg19", name, self._contigs["hg19"][name]))

        for primary_key, (build, name, variants) in enumerate(contigs):
            self._print(
                "INSERT INTO [Contigs] VALUES ({}, {}, {}, {});",
                primary_key,
                self._to_string(build),
                self._to_string(name),
                self._to_string(variants),
            )

    @staticmethod
    def _build_consequence_ranks():
        """Returns consequences with a human friendly ranking: bad > insignificant."""
        human_friendly_ranks = collections.OrderedDict()
        for rank, name in enumerate(reversed(list(consequences.ranks()))):
            human_friendly_ranks[name] = rank

        return human_friendly_ranks

    @staticmethod
    def _to_string(value):
        if isinstance(value, (int, float)):
            return repr(value)
        elif isinstance(value, (tuple, list)):
            if not value:
                return "NULL"

            value = ";".join(map(str, value))
        elif value is None:
            return "NULL"
        else:
            value = str(value)

        return "'{}'".format(value.replace("'", "''"))


class SQLite3Output(SQLOutput):
    def __init__(self, keys, out_prefix):
        self._conn = sqlite3.connect(f"{out_prefix}.sqlite3")
        self._curs = self._conn.cursor()

        super().__init__(keys, None)

    def finalize(self):
        super().finalize()

        self._conn.commit()
        self._conn.close()

    def _print(self, line="", *args, **kwargs):
        if args:
            line = line.format(*args)

        query = line.strip()
        if query:
            self._curs.execute(query)


FORMATS = {
    "json": JSONOutput,
    "tsv": TSVOutput,
    "sql": SQLOutput,
    "sqlite3": SQLite3Output,
}