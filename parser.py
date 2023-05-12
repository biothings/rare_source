import os
from copy import deepcopy
from typing import Dict, Union
import pandas as pd


def read_rare_diseases(filepath) -> pd.DataFrame:
    separator = ","
    na_value = ""
    quotechar = '"'  # double quote
    column_info = [
        # Each element is a tuple of (column_index, column_name, data_type)
        (0, "Rare_Disease_Name", "string"),  # column 0 (the disease's name)
        (1, "Disease_Aliases", "string"),  # column 1 (a string separated by double slashes)
        (2, "Associated_Genes", "string"),  # column 2 (a string of gene symbols, separated by semicolons)
        (3, "Disease_Annotations", "string"),  # column 3 (a URL)
        (4, "OMIM", "string"),  # column 4 (OMIM ID of the disease, a six-digit number)
        (5, "Orphanet", "string"),  # column 5 (Orphanet ID, a.k.a. ORPHAcode of the disease, numerical)
        (6, "UMLS", "string"),  # column 6 (UMLS CUI of the disease)
        (7, "Mesh", "string"),  # column 7 (MeSH ID of the disease)
        (8, "ICD10CM", "string"),  # column 8 (ICD-10-CM ID of the disease, a seven-character, alphanumeric code)
        # (9, "Gene_Descriptions", "string"),  # column 9 (a string of gene descriptions, separated by semicolons; see column 2)
        # (10, "Links", "string"),  # column 10 (a URL)
    ]
    column_indices = [e[0] for e in column_info]
    column_names = [e[1] for e in column_info]
    column_dtypes = {e[1]: e[2] for e in column_info}

    # Use header=0 to ignore the original header
    data_frame = pd.read_csv(filepath, sep=separator, header=0, names=column_names, usecols=column_indices,
                             dtype=column_dtypes, na_values=[na_value], quotechar=quotechar)

    gene_seperator = ";"
    assert data_frame["Associated_Genes"].str.endswith(gene_seperator).all()
    data_frame["Associated_Genes"] = data_frame["Associated_Genes"].str[:-1]  # remove the tailing semicolon
    data_frame["Associated_Genes"] = data_frame["Associated_Genes"].str.split(gene_seperator)

    # Each disease annotation is a URL in "https://raresource.nih.gov/literature/disease/<GARD_ID>" format
    # GARD is proposed by Genetic and Rare Diseases Information Center of NIH
    # See https://rarediseases.info.nih.gov/about for more details on GARD
    url_path = "https://raresource.nih.gov/literature/disease/"
    assert data_frame["Disease_Annotations"].str.startswith(url_path).all()
    assert not data_frame["Disease_Annotations"].str[len(url_path):].duplicated().any()  # ensure the the GARD IDs are unique
    data_frame["GARD_ID"] = data_frame["Disease_Annotations"].str[len(url_path):]  # remove the heading url path, get GARD IDs

    # Do not split Disease_Aliases here; leave it to get_rare_disease_register()
    # alias_separator = "//"
    # data_frame["Disease_Aliases"] = data_frame["Disease_Aliases"].str.split(alias_separator)

    return data_frame


def construct_rare_disease_entity(rare_disease_row: pd.Series) -> Dict:
    entity = dict()

    # These 3 columns contains No <NA> values
    entity["name"] = rare_disease_row["Rare_Disease_Name"]
    entity["annotation_url"] = rare_disease_row["Disease_Annotations"]
    entity["gard"] = rare_disease_row["GARD_ID"]

    # The following columns may contain <NA> values, so pd.isna() checks are necessary
    if not pd.isna(rare_disease_row["Disease_Aliases"]):
        aliases = rare_disease_row["Disease_Aliases"].split("//")
        if aliases[-1] == "":
            aliases = aliases[0:-1]  # remove the last empty alias, because Disease_Aliases may end with double slashes
        entity["alias"] = aliases

    if not pd.isna(rare_disease_row["OMIM"]):
        entity["omim"] = rare_disease_row["OMIM"]

    if not pd.isna(rare_disease_row["Orphanet"]):
        entity["orphanet"] = rare_disease_row["Orphanet"]

    if not pd.isna(rare_disease_row["UMLS"]):
        entity["umls"] = rare_disease_row["UMLS"]

    if not pd.isna(rare_disease_row["Mesh"]):
        entity["mesh"] = rare_disease_row["Mesh"]

    if not pd.isna(rare_disease_row["ICD10CM"]):
        entity["icd10cm"] = rare_disease_row["ICD10CM"]

    return entity


def get_rare_disease_register(rare_disease_df: pd.DataFrame) -> Dict:
    register = {}  # <GARD_ID : Rare_Disease_Entity>

    for _, row in rare_disease_df.iterrows():
        # These 3 fields should not contain NaN
        rare_disease = construct_rare_disease_entity(row)
        register[row["GARD_ID"]] = rare_disease

    return register


def get_gene_disease_mapping(rare_disease_df: pd.DataFrame) -> Dict:
    mapping = {}  # <Gene_Symbol : GARD_ID_List>

    for _, row in rare_disease_df.iterrows():
        for gene_symbol in row["Associated_Genes"]:
            if gene_symbol not in mapping:
                mapping[gene_symbol] = []

            mapping[gene_symbol].append(row["GARD_ID"])

    return mapping


def read_genes(filepath) -> pd.DataFrame:
    separator = ","
    na_value = ""
    quotechar = '"'  # double quote
    column_info = [
        # Each element is a tuple of (column_index, column_name, data_type)
        (0, "Gene_Information", "string"),  # column 0 (the gene's symbol)
        (1, "Description", "string"),  # column 1
        # (2, "Associated_Rare_Diseases", "string"),  # column 2 (a string of <disease_name?gard_id> sequences, separated by semicolons)
        (3, "Gene_Annotations", "string"),  # column 3 (a URL)
        (4, "Gene_IDs", "string"),  # column 4 (Entrez ID, integer)
        (5, "Ensembl_Gene_ID", "string"),  # column 5 (Ensembl ID, alphanumeric)
        (6, "HGNC_ID", "string"),  # column 6 (HGNC ID, numeric)
        # (7, "Links", "string"),  # column 7 (a URL)
    ]
    column_indices = [e[0] for e in column_info]
    column_names = [e[1] for e in column_info]
    column_dtypes = {e[1]: e[2] for e in column_info}

    # Use header=0 to ignore the original header
    data_frame = pd.read_csv(filepath, sep=separator, header=0, names=column_names, usecols=column_indices,
                             dtype=column_dtypes, na_values=[na_value], quotechar=quotechar)

    return data_frame


def construct_gene_entity(gene_row: pd.Series) -> Dict:
    entity = dict()

    # These 4 columns contains No <NA> values
    entity["symbol"] = gene_row["Gene_Information"]
    entity["description"] = gene_row["Description"]
    entity["annotation_url"] = gene_row["Gene_Annotations"]
    entity["entrezgene"] = gene_row["Gene_IDs"]

    # The following columns may contain <NA> values, so pd.isna() checks are necessary
    if not pd.isna(gene_row["Ensembl_Gene_ID"]):
        entity["ensemblgene"] = gene_row["Ensembl_Gene_ID"]

    if not pd.isna(gene_row["HGNC_ID"]):
        entity["hgnc"] = gene_row["HGNC_ID"]

    return entity


def construct_doc(gene_entity: Dict, gene_disease_mapping: Dict, rare_disease_register: Dict) -> Union[Dict, None]:
    gard_ids = gene_disease_mapping.get(gene_entity["symbol"], None)
    if gard_ids is None:
        return None

    disease_entities = [rare_disease_register[gid] for gid in gard_ids if gid in rare_disease_register]
    if not disease_entities:
        return None

    disease_entities = deepcopy(disease_entities)
    for disease_entity in disease_entities:
        disease_entity["cooccurrence_url"] = f"https://raresource.nih.gov/literature/cooccurrence/{gene_entity['symbol']}/{disease_entity['gard']}"

    doc = {
        "_id": gene_entity["entrezgene"],
        **gene_entity,
        "raresource": {
            "disease": disease_entities
        }
    }

    return doc


def load_data(data_folder):
    rare_disease_filename = "RARe-SOURCE-Browse-Rare-Diseases-03-21-2023.csv"
    rare_disease_filepath = os.path.join(data_folder, rare_disease_filename)
    gene_filename = "RARe-SOURCE-Browse-Genes-03-21-2023.csv"
    gene_filepath = os.path.join(data_folder, gene_filename)

    rare_disease_df = read_rare_diseases(rare_disease_filepath)
    rare_disease_register = get_rare_disease_register(rare_disease_df)
    gene_disease_mapping = get_gene_disease_mapping(rare_disease_df)

    gene_df = read_genes(gene_filepath)
    for _, row in gene_df.iterrows():
        gene_entity = construct_gene_entity(row)
        doc = construct_doc(gene_entity=gene_entity, gene_disease_mapping=gene_disease_mapping, rare_disease_register=rare_disease_register)
        if doc is not None:
            yield doc
