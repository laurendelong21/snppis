# -*- coding: utf-8 -*-

"""Full pipeline for getting patient-pathway score dataframe."""

import json
import logging
import os
import re
from typing import Dict, List, Mapping, Set

import pandas as pd
from tqdm import tqdm

from .api import batched_query
from .constants import DATABASES, PathwayTuple
from .lookup_predictions import load_pathway_to_snps
from .report_impactful import check_impacted

__all__ = [
    'get_pathway_to_patient_to_score_df',
]

DBSNP_RE = re.compile(r'^rs\d+$')

logger = logging.getLogger(__name__)


def get_pathway_to_patient_to_score_df(df: pd.DataFrame, output_directory: str) -> pd.DataFrame:
    """Get the patient to pathway score dataframe.

    :param df: A dataframe in which the rows represent patients and columns are SNPs. Entries
     are 0, 1, or 2 that count how many times that patient had the SNP.
    :param output_directory: The place where results are stored
    """
    universe = set(df.columns[1:])
    if not all(DBSNP_RE.match(dbsnp_id) for dbsnp_id in universe):
        raise ValueError('Some invalid SNP labels used')
    logging.info(f'Data has {len(universe)} SNPs as columns')

    logger.info('Getting SNPs info')
    dbsnp_impacted: Mapping[str, bool] = {
        result['dbsnp']['rsid']: check_impacted(result)
        for result in batched_query(universe)
    }

    logger.info(f'Building mappings for {", ".join(DATABASES)}')
    pathway_to_snps: Mapping[PathwayTuple, List[str]] = {
        pathway: snps
        for db in DATABASES
        for pathway, snps in load_pathway_to_snps(db).items()
    }

    logger.info('Calculating pathway scores')

    logger.info('Filtering SNPs in each pathway to those in the universe')
    pathway_to_relevant: Dict[PathwayTuple, Set[str]] = {}
    for pathway, snps in pathway_to_snps.items():
        intersection = universe.intersection(snps)
        if not intersection:
            logger.debug('No SNPs in %s:%s ! %s', *pathway)
            continue
        pathway_to_relevant[pathway] = intersection

    logger.info('Get all damaged SNPs from this list')
    pathway_to_impactful_relevant: Mapping[PathwayTuple, Set[str]] = {
        pathway: {
            snp
            for snp in snps
            if snp in dbsnp_impacted
        }
        for pathway, snps in pathway_to_relevant.items()
    }

    with open(os.path.join(output_directory, 'impactful_relevant.json'), 'w') as file:
        json.dump(
            [
                dict(db=db, id=db_id, label=db_label, snps=sorted(snps))
                for (db, db_id, db_label), snps in pathway_to_impactful_relevant.items()
                if snps
            ],
            file,
            indent=2,
        )

    pathway_to_snps_it = tqdm(pathway_to_snps.items(), desc='Scoring pathways')

    r = []
    for (db, db_id, db_label), pathway_snps in pathway_to_snps_it:
        pathway_curie = f'{db}:{db_id}'

        try:
            relevant_pathway_snps: Set[str] = pathway_to_relevant[db, db_id, db_label]
            pathway_to_snps_it.write(f'Found: {(db, db_id, db_label)}')
        except KeyError:
            pathway_to_snps_it.write(f'Could not find: {(db, db_id, db_label)}')
            continue

        impactful_relevant_pathway_snps: Set[str] = pathway_to_impactful_relevant[db, db_id, db_label]
        n_relevant_pathway_snps = len(relevant_pathway_snps)
        for patient, *patient_snps in tqdm(df.values, leave=False, desc=f'{db}:{db_id}'):
            score = sum(
                count  # this is either 0, 1, or 2s
                for patient_dbsnp_id, count in zip(df.columns[1:], patient_snps)
                if patient_dbsnp_id in impactful_relevant_pathway_snps
            ) / n_relevant_pathway_snps
            r.append((patient, pathway_curie, score))

    pathway_df = pd.DataFrame(r, columns=['patient', 'pathway', 'score'])
    pathway_df_path = os.path.join(output_directory, 'scores.tsv')
    logger.info(f"Outputting scores to {pathway_df_path}")
    pathway_df.to_csv(pathway_df_path, sep='\t', index=False)

    return pathway_df
