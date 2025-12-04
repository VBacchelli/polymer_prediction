from schema import JSON_SCHEMA
import copy

class BiceranoWrapper:

    COLUMN_MAP = {
        "Polymer": "polymer_name",
        "SMILES": "smiles",
        "BigSMILES": "bigsmiles",
        "Tg (K) exp": ("properties", "Tg"),
    }

    def convert_row(self, row, row_index):
        doc = copy.deepcopy(JSON_SCHEMA)

        for col in row.index:
            if col in self.COLUMN_MAP:
                target = self.COLUMN_MAP[col]
                
                if isinstance(target, tuple):
                    # È un campo annidato, es. Tg dentro a properties
                    doc[target[0]][target[1]] = row[col]
                else:
                    doc[target] = row[col]

        # aggiungi campi generali
        doc["source_dataset"] = "Bicerano_bigsmiles"

        return doc
