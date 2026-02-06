# Polymer Dataset ETL

This folder contains the ETL (Extract, Transform, Load) pipeline used to convert
heterogeneous polymer datasets into a unified JSON format and store them in MongoDB.

Each source dataset is processed independently and mapped to a common schema,
so that data coming from different sources can be queried and used consistently.

## Overview

The ETL pipeline works as follows:

1. Load a dataset (CSV or ZIP, local path or URL)
2. Convert each row into a normalized JSON document
3. Store the resulting documents in MongoDB

A shared schema defines the target structure of each document, while a dedicated
wrapper is implemented for each dataset to handle column mapping and dataset-specific logic.

## Schema

A common base schema is used to represent polymer data in a uniform way
(e.g. molecular representations and properties).

This schema is defined once in the file [`schema.py`](./datasetToJson/schema.py) and reused across all datasets.

## Dataset wrappers

Each dataset has its own wrapper class, responsible for:
- mapping dataset-specific columns to the common schema
- handling missing or renamed fields
- adding metadata such as the source dataset name

This makes it easy to add new datasets without modifying the core ETL logic. The wrappers used can be found in the [wrappers](/datasetToJson/wrappers) folder.

## MongoDB

The converted documents are stored in MongoDB, with one collection per dataset.
A [repository layer](./datasetToJson/mongoRepository.py) is used to handle insertion and collection management.

The MongoDB connection is configured via environment variables.

## Run

A typical ETL run can be started with:

```bash
python main.py
```
