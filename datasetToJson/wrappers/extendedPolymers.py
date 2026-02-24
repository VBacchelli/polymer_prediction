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
        row = dict(row)

        # Conversione Tg da °C a K (se presente e non nullo)
        tg_value = row.get("Tg")
        if tg_value is not None:
            try:
                row["Tg"] = float(tg_value) + 273.15
            except (ValueError, TypeError):
                row["Tg"] = None  

        return convert_row(
            row=row,
            column_map=self.COLUMN_MAP,
            base_schema=JSON_SCHEMA,
            extra_fields={"source_dataset": self.SOURCE}
        )
