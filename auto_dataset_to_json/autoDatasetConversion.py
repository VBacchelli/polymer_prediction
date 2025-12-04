import os
import zipfile
import requests
import pandas as pd
import json

# ------ DOWNLOAD DA URL -----
def download_file(url):
    print("Downloading dataset...")
    r = requests.get(url, stream=True)
    r.raise_for_status()

    # Prova a leggere il nome del file dall'header
    cd = r.headers.get("Content-Disposition")
    if cd and "filename=" in cd:
        filename = cd.split("filename=")[1].strip('"')
    else:
        # fallback: usa il nome alla fine dell’URL
        filename = url.split("/")[-1]

    with open(filename, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"Saved file as {filename}")
    return filename

def identify_file_type(filepath):
    _, ext = os.path.splitext(filepath)

    ext = ext.lower()

    if ext == ".zip":
        return "zip"
    if ext == ".csv":
        return "csv"

    raise ValueError(f"Unsupported file format: {ext}")

# ------ ESTRAZIONE ZIP ------------
def extract_zip(zip_path, extract_to="extracted"):
    print("Extracting...")
    os.makedirs(extract_to, exist_ok=True)
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
    print("Extracted to", extract_to)
    return extract_to

# ------ LETTURA DEL PRIMO CSV TROVATO NELLA CARTELLA IN INPUT ------
def load_first_csv(folder):
    for f in os.listdir(folder):
        if f.endswith(".csv"):
            csv_path = os.path.join(folder, f)
            print("Found CSV:", csv_path)
            return pd.read_csv(csv_path, encoding="latin1")
    raise FileNotFoundError("No CSV found in extracted ZIP.")

# ------- GESTIONE DEL FILE IN BASE AL TIPO (ZIP vs CSV vs ALTRO) -------
def load_dataset_from_file(filepath):
    filetype = identify_file_type(filepath)

    if filetype == "csv":
        print("Detected CSV file.")
        return pd.read_csv(filepath, encoding="latin1")

    elif filetype == "zip":
        print("Detected ZIP file. Extracting...")
        extract_dir = "extracted"
        extract_zip(filepath, extract_dir)
        return load_first_csv(extract_dir)

    else:
        raise ValueError(f"Unknown file type for {filepath}")

# ----- CONVERSIONE IN JSON -------
def create_json_documents(df, output_folder="json_output"):
    os.makedirs(output_folder, exist_ok=True)

    for i, row in df.iterrows():
        doc = {
            "polymer_name": row.get("Polymer"),
            "smiles": row.get("SMILES"),
            "bigsmiles": row.get("BigSMILES"),
            "properties": {
                "Tg": row.get("Tg (K) exp"),
                "unit": "K",
            },
            "source_dataset": "Bicerano"
        }

        out_path = os.path.join(output_folder, f"polymer_{i}.json")
        with open(out_path, "w") as f:
            json.dump(doc, f, indent=4)

    print(f"Generated {len(df)} JSON files in {output_folder}")

if __name__ == "__main__":
    url = "https://springernature.figshare.com/ndownloader/files/42507037"

    downloaded_file = download_file(url)
    df = load_dataset_from_file(downloaded_file)

    print(df.head())
    create_json_documents(df)
