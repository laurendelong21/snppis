# -*- coding: utf-8 -*-


import re

import pandas as pd

from .api import batched_query
from .constants import DATABASES
from .lookup_predictions import load_pathway_to_snps
from .report_impactful import check_impacted

__all__ = [
    'get_pathway_to_patient_to_score',
]

DBSNP_RE = re.compile(r'^rs\d+$')


def get_pathway_to_patient_to_score(df: pd.DataFrame):
    """

    :param df: A dataframe in which the rows represent patients and columns are SNPs. Entries
     are 0, 1, or 2 that count how many times that patient had the SNP.
    """
    universe = set(df.columns[1:])
    if not all(DBSNP_RE.match(dbsnp_id) for dbsnp_id in universe):
        raise ValueError('Some invalid SNP labels used')

    # TODO needs caching or ability to load from file
    dbsnp_impacted = {
        result['dbsnp']['rsid']: check_impacted(result)
        for result in batched_query(universe)
    }

    pathway_to_snps = {
        (db, db_id, db_label): snps
        for db in DATABASES
        for (db_id, db_label), snps in load_pathway_to_snps(db).items()
    }

    pathway_to_patient_to_score = {}

    for pathway, pathway_snps in pathway_to_snps.items():
        # 1. Get all SNPs in the pathway that are also in the universe
        relevant_pathway_snps = universe.intersection(pathway_snps)
        n_relevant_pathway_snps = len(relevant_pathway_snps)

        # 2. Get all damaged SNPs from this list
        impactful_relevant_pathway_snps = {
            relevant_snp
            for relevant_snp in relevant_pathway_snps
            if relevant_snp in dbsnp_impacted
        }

        pathway_to_patient_to_score[pathway] = {
            sum(
                count  # this is either 0, 1, or 2s
                for patient_dbsnp_id, count in zip(df.columns[1:], patient_snps)
                if patient_dbsnp_id in impactful_relevant_pathway_snps
            ) / n_relevant_pathway_snps
            for patient, *patient_snps in df.values
        }

    return pathway_to_patient_to_score
