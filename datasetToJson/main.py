import os
from dotenv import load_dotenv

from etl import load_dataset, dataframe_to_documents
from wrappers.densityDataset import DensityWrapper
from mongoRepository import MongoRepository

load_dotenv()

MONGO_URI = os.getenv("MONGO_URI")

#URL = "https://huggingface.co/datasets/AdrianM0/bicerano_polymers/resolve/main/HT_MD_polymer_properties.csv?download=true"
PATH="./extracted/1_polymer_density_dataset.csv"

# ETL
df = load_dataset(PATH)
wrapper = DensityWrapper()
docs = dataframe_to_documents(df, wrapper)

# Mongo
repo = MongoRepository(
    uri=MONGO_URI,
    db_name="PolymerPrediction",
    collection_name=wrapper.SOURCE
)

n = repo.insert_documents(docs, reset=True)
print(f"Inserted {n} documents")

repo.close()
