def rare_source_mapping(cls):
    mapping = {
        "symbol": {
            "normalizer": "keyword_lowercase_normalizer",
            "type": "keyword"
        },
        "description": {
            "type": "text"
        },
        "annotation_url": {
            "type": "text",
            "index": False
        },
        "entrezgene": {
            "type": "keyword"
        },
        "ensemblgene": {
            "type": "keyword"
        },
        "hgnc": {
            "type": "keyword"
        },
        "raresource": {
            "properties": {
                "disease": {
                    "properties": {
                        "name": {
                            "type": "text"
                        },
                        "alias": {
                            "type": "text"
                        },
                        "annotation_url": {
                            "type": "text",
                            "index": False
                        },
                        "cooccurrence_url": {
                            "type": "text",
                            "index": False
                        },
                        "gard": {
                            "type": "keyword"
                        },
                        "omim": {
                            "type": "keyword"
                        },
                        "orphanet": {
                            "type": "keyword"
                        },
                        "umls": {
                            "type": "keyword"
                        },
                        "mesh": {
                            "type": "keyword"
                        },
                        "icd10cm": {
                            "type": "keyword"
                        }
                    }
                }
            }
        }
    }

    return mapping
