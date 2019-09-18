# -*- coding: utf-8 -*-

"""Full pipeline for getting patient-pathway score dataframe."""

import logging
import re

import pandas as pd
from tqdm import tqdm

from .api import batched_query
from .constants import DATABASES
from .lookup_predictions import load_pathway_to_snps
from .report_impactful import check_impacted

__all__ = [
    'get_pathway_to_patient_to_score',
]

DBSNP_RE = re.compile(r'^rs\d+$')

logger = logging.getLogger(__name__)


def get_pathway_to_patient_to_score(df: pd.DataFrame):
    """Get the patient to pathway score dataframe.

    :param df: A dataframe in which the rows represent patients and columns are SNPs. Entries
     are 0, 1, or 2 that count how many times that patient had the SNP.
    """
    universe = set(df.columns[1:])
    if not all(DBSNP_RE.match(dbsnp_id) for dbsnp_id in universe):
        raise ValueError('Some invalid SNP labels used')
    logging.info(f'Data has {len(universe)} SNPs as columns')

    logger.info('Looking up SNPs in MyVariant')
    dbsnp_impacted = {
        result['dbsnp']['rsid']: check_impacted(result)
        for result in batched_query(universe)
    }

    logger.info(f'Building mappings for {", ".join(DATABASES)}')
    pathway_to_snps = {
        (db, db_id, db_label): snps
        for db in DATABASES
        for (db_id, db_label), snps in load_pathway_to_snps(db).items()
    }

    logger.info('Calculating pathway scores')
    pathway_to_patient_to_score = {}
    for pathway, pathway_snps in tqdm(pathway_to_snps.items(), desc='Scoring pathway'):
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
