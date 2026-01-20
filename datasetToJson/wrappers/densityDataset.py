from schema import JSON_SCHEMA
from etl import convert_row

class DensityWrapper:
    COLUMN_MAP = {
        "SMILES": "smiles",
        "Density": ("properties", "density")
    }

    SOURCE = "polymerDensity"

    def convert_row(self, row):
        return convert_row(
            row=row,
            column_map=self.COLUMN_MAP,
            base_schema=JSON_SCHEMA,
            extra_fields={"source_dataset": self.SOURCE}
        )
