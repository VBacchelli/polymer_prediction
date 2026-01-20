from schema import JSON_SCHEMA
from etl import convert_row

class ExtendedWrapper:
    COLUMN_MAP = {
        "SMILES": "smiles",
        "Tg": ("properties", "Tg"),
        "Density": ("properties", "density"),
        "Tc": ("properties", "Tc")
    }

    SOURCE = "extendedPolymers"

    def convert_row(self, row):
        return convert_row(
            row=row,
            column_map=self.COLUMN_MAP,
            base_schema=JSON_SCHEMA,
            extra_fields={"source_dataset": self.SOURCE}
        )
