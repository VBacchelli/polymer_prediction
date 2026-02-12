from schema import JSON_SCHEMA
from etl import convert_row

class BiceranoHfWrapper:
    SOURCE = "biceranoPolymers"

    COLUMN_MAP = {
        "Poymer name": "polymer_name",
        "SMILES (Atoms Ce and Th are placeholders for head and tail information, respectively)": "smiles",
        "BigSMILES": "bigsmiles",
        # proprietà
        "Experiment Tg (K)": ("properties", "Tg"),
        "Experiment density at 300K (g/cc)": ("properties", "density"),
        "Calculated glass CTE (10-6/K)": ("properties", "glass_cte"),
        "Calculated rubber CTE (10-6/K)": ("properties", "rubber_cte"),
    }

    def convert_row(self, row):
        return convert_row(
            row=row,
            column_map=self.COLUMN_MAP,
            base_schema=JSON_SCHEMA,
            extra_fields={"source_dataset": self.SOURCE}
        )
