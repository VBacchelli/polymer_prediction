from schema import JSON_SCHEMA
from etl import convert_row

class BiceranoWrapper:
    COLUMN_MAP = {
        "Polymer": "polymer_name",
        "SMILES": "smiles",
        "BigSMILES": "bigsmiles",
        "Tg (K) exp": ("properties", "Tg"),
    }

    SOURCE = "BiceranoTg"

    def convert_row(self, row):
        return convert_row(
            row=row,
            column_map=self.COLUMN_MAP,
            base_schema=JSON_SCHEMA,
            extra_fields={"source_dataset": self.SOURCE}
        )
