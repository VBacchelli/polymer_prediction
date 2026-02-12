import os
import zipfile
import requests
import pandas as pd
import copy
import re

# ============================================================
# UTILS
# ============================================================
def normalize_star(smiles: str) -> str:
    smiles = re.sub(r"\(\s*\*\s*\)", "[*]", smiles)
    smiles = smiles.replace("[Ce]", "[*]").replace("[Th]", "[*]")
    smiles = re.sub(r"(?<!\[)\*(?!\])", "[*]", smiles)

    return smiles


# ============================================================
# LOADERS
# ============================================================

def download_file(url: str) -> str:
    """
    Scarica un file da URL e restituisce il path locale.
    """
    print("Downloading dataset...")
    r = requests.get(url, stream=True)
    r.raise_for_status()

    # prova a ottenere il nome file dall'header
    cd = r.headers.get("Content-Disposition")
    if cd and "filename=" in cd:
        filename = cd.split("filename=")[1].strip('"')
    else:
        filename = url.split("/")[-1]

    # sanitizza nome file
    filename = filename.replace('"', "").replace(";", "").strip()

    with open(filename, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"Saved file as {filename}")
    return filename


def identify_file_type(filepath: str) -> str:
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()

    if ext == ".zip":
        return "zip"
    if ext == ".csv":
        return "csv"

    raise ValueError(f"Unsupported file format: {ext}")


def extract_zip(zip_path: str, extract_to: str = "extracted") -> str:
    """
    Estrae uno ZIP e restituisce la cartella di output.
    """
    print("Extracting ZIP...")
    os.makedirs(extract_to, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(extract_to)

    return extract_to


def load_first_csv(folder: str) -> pd.DataFrame:
    """
    Cerca e carica il primo CSV trovato in una cartella.
    """
    for f in os.listdir(folder):
        if f.lower().endswith(".csv"):
            path = os.path.join(folder, f)
            print(f"Found CSV: {path}")
            return pd.read_csv(path, encoding="latin1")

    raise FileNotFoundError("No CSV found in extracted folder.")


def load_dataset(source: str) -> pd.DataFrame:
    """
    source può essere:
    - path locale a CSV
    - path locale a ZIP
    - URL
    """
    if source.startswith("http://") or source.startswith("https://"):
        source = download_file(source)

    filetype = identify_file_type(source)

    if filetype == "csv":
        print("Detected CSV file.")
        return pd.read_csv(source, encoding="latin1")

    if filetype == "zip":
        print("Detected ZIP file.")
        extract_dir = extract_zip(source)
        return load_first_csv(extract_dir)

    raise ValueError(f"Cannot load dataset from {source}")


# ============================================================
# CONVERSION
# ============================================================

def convert_row(
    row: pd.Series,
    column_map: dict,
    base_schema: dict,
    extra_fields: dict | None = None
) -> dict:
    """
    Converte una riga del DataFrame in un documento JSON normalizzato.
    """
    doc = copy.deepcopy(base_schema)

    for col, target in column_map.items():
        if col not in row:
            continue

        value = row[col]

        # normalizza asterischi in smiles
        if target == "smiles":
            value = normalize_star(value)

        # salta NaN
        if pd.isna(value):
            continue

        if isinstance(target, tuple):
            parent, child = target
            doc[parent][child] = value
        else:
            doc[target] = value

    if extra_fields:
        doc.update(extra_fields)

    return doc


# ============================================================
# PIPELINE
# ============================================================

def dataframe_to_documents(df: pd.DataFrame, wrapper):
    """
    Generatore di documenti JSON a partire da un DataFrame e un wrapper.
    """
    for _, row in df.iterrows():
        yield wrapper.convert_row(row)


def dataframe_to_json_files(
    df: pd.DataFrame,
    wrapper,
    output_folder: str = "json_output"
):
    """
    Converte un DataFrame in file JSON su disco.
    """
    os.makedirs(output_folder, exist_ok=True)

    for i, doc in enumerate(dataframe_to_documents(df, wrapper)):
        path = os.path.join(output_folder, f"polymer_{i}.json")
        with open(path, "w") as f:
            import json
            json.dump(doc, f, indent=4)

    print(f"Generated {i + 1} JSON files in {output_folder}")
