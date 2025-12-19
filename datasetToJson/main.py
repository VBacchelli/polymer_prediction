import os
from dotenv import load_dotenv

from etl import load_dataset, dataframe_to_documents
from wrappers.biceranoTg import BiceranoWrapper
from mongoRepository import MongoRepository

load_dotenv()

MONGO_URI = os.getenv("MONGO_URI")

URL = "https://springernature.figshare.com/ndownloader/files/42507037"

# ETL
df = load_dataset(URL)
wrapper = BiceranoWrapper()
docs = dataframe_to_documents(df, wrapper)

# Mongo
repo = MongoRepository(
    uri=MONGO_URI,
    db_name="PolymerPrediction",
    collection_name="bicerano_tg"
)

n = repo.insert_documents(docs, reset=True)
print(f"Inserted {n} documents")

repo.close()
