from pymongo import MongoClient
from typing import Iterable
import os


class MongoRepository:
    """
    Repository MongoDB per il caricamento di documenti JSON.
    """

    def __init__(
        self,
        uri: str,
        db_name: str,
        collection_name: str
    ):
        self.uri = uri
        self.db_name = db_name
        self.collection_name = collection_name

        self.client = MongoClient(self.uri)
        self.collection = self.client[self.db_name][self.collection_name]

    def insert_documents(
        self,
        documents: Iterable[dict],
        reset: bool = False
    ) -> int:
        """
        Inserisce una sequenza di documenti in MongoDB.

        Parameters
        ----------
        documents : Iterable[dict]
            Documenti da inserire (lista o generatore)
        reset : bool
            Se True, svuota la collezione prima dell'inserimento

        Returns
        -------
        int
            Numero di documenti inseriti
        """

        if reset:
            self.collection.delete_many({})

        docs = list(documents)

        if not docs:
            return 0

        result = self.collection.insert_many(docs)
        return len(result.inserted_ids)

    def count(self) -> int:
        """
        Restituisce il numero di documenti nella collezione.
        """
        return self.collection.count_documents({})

    def close(self):
        """
        Chiude la connessione MongoDB.
        """
        self.client.close()
