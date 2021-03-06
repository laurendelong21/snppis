# -*- coding: utf-8 -*-

"""Look up the SNP from each pathway in MyVariant."""

import json
import logging
import os
from typing import Mapping, Set

import click

from snppis.api import batched_query
from snppis.constants import DATABASES, MAPPINGS, PathwayTuple

logger = logging.getLogger(__name__)


@click.command()
def main():
    """Generate scores for each pathway."""
    for db in DATABASES:
        lookup_predictions(db)


def load_pathway_to_snps(db: str) -> Mapping[PathwayTuple, Set[str]]:
    """Load mapping from pathway to SNPs for the given database."""
    path = os.path.join(MAPPINGS, f'{db}.json')
    logger.info(f'Loading {db} mappings from {path}')
    with open(path) as file:
        pathway_to_snp = json.load(file)

    logger.info(f'Reorganizing {db} mappings')
    return {
        (
            db,
            entry['pathway']['identifier'],
            entry['pathway']['name'],
        ): {
            snp['identifier']
            for snp in entry['snps']
        }
        for entry in pathway_to_snp
    }


def lookup_predictions(db: str):
    """Look up the SNP from each pathway in MyVariant."""
    pathway_to_snps = load_pathway_to_snps(db)

    dbsnp_ids: Set[str] = {
        snp
        for pathway, snps in pathway_to_snps.items()
        for snp in snps
    }

    click.echo(f'Got {len(dbsnp_ids)} SNPs from {db}')

    results = list(batched_query(dbsnp_ids))
    path = os.path.join(MAPPINGS, f'{db}_snp_scores.json')
    with open(path, 'w') as file:
        json.dump(results, file, indent=2)


if __name__ == '__main__':
    main()
