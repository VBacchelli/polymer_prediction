from schema import JSON_SCHEMA
from etl import convert_row

class PolyInfoWrapper:
    SOURCE = "polyInfo_kaggle"

    COLUMN_MAP = {
        "SMILES": "smiles",
        # proprietà
        "Tg": ("properties", "Tg")
    }

    def convert_row(self, row):
        doc = convert_row(
            row=row,
            column_map=self.COLUMN_MAP,
            base_schema=JSON_SCHEMA,
            extra_fields={"source_dataset": self.SOURCE}
        )

        # conversione Celsius → Kelvin
        tg_c = doc.get("properties", {}).get("Tg")

        if tg_c is not None:
            try:
                doc["properties"]["Tg"] = float(tg_c) + 273.15
            except (ValueError, TypeError):
                return None  # scarta riga

        return doc
