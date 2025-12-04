import os
import json
from pymongo import MongoClient

def upload_json_folder_to_mongo(folder, mongo_uri, db_name, collection_name):
    # Connessione al database
    client = MongoClient(mongo_uri)
    collection = client[db_name][collection_name]

    # Carica i JSON da cartella
    docs = []
    for filename in os.listdir(folder):
        if filename.endswith(".json"):
            path = os.path.join(folder, filename)
            with open(path, "r") as f:
                doc = json.load(f)
                docs.append(doc)

    # Inserimento
    if docs:
        result = collection.insert_many(docs)
        print(f"Inserted {len(result.inserted_ids)} documents into MongoDB.")
    else:
        print("No JSON documents found to upload.")
    client.close()

if __name__ == "__main__":
    FOLDER = "json_output"
    MONGO_URI = ""
    DB_NAME = "polymer_db"
    COLLECTION_NAME = "bicerano_tg"

    upload_json_folder_to_mongo(FOLDER, MONGO_URI, DB_NAME, COLLECTION_NAME)
