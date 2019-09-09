# -*- coding: utf-8 -*-

"""Report the impactful SNPs in each pathway."""

import json
import os

import click

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES = os.path.join(HERE, 'mappings')

POLYPHEN2_DAMAGING = {"damaging", "possibly damaging", "probably damaging"}
SIFT_DAMAGING = {"deleterious", "deleterious - low confidence"}


@click.command()
def main():
    """Generate scores for each pathway."""
    for db in ('kegg', 'wikipathways', 'reactome'):
        report_impactful(db)


def report_impactful(db: str):
    """Report the impactful SNPs in the given database."""
    with open(os.path.join(RESOURCES, f'{db}_snp_scores.json')) as file:
        scores = json.load(file)

    for entry in scores:
        dbsnp_id = entry['dbsnp']['rsid']
        cadd = entry.get('cadd')
        if cadd is None:
            continue

        polyphen = cadd.get('polyphen', {})
        if isinstance(polyphen, list):
            polyphen = polyphen[0]

        sift = cadd.get('sift', {})
        if isinstance(sift, list):
            sift = sift[0]

        impact = polyphen.get('cat') in POLYPHEN2_DAMAGING or sift.get('cat') in SIFT_DAMAGING

        click.echo(f"{db}\t{dbsnp_id}\t{'⚠️' if impact else ''}")


if __name__ == '__main__':
    main()
