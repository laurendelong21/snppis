# -*- coding: utf-8 -*-

"""Report the impactful SNPs in each pathway."""

import json
import os

import click

from snppis.constants import DATABASES, MAPPINGS

POLYPHEN2_DAMAGING = {"damaging", "possibly damaging", "probably damaging"}
SIFT_DAMAGING = {"deleterious", "deleterious - low confidence"}


@click.command()
def main():
    """Generate scores for each pathway."""
    for db in DATABASES:
        report_impactful(db)


def report_impactful(db: str):
    """Report the impactful SNPs in the given database."""
    with open(os.path.join(MAPPINGS, f'{db}_snp_scores.json')) as file:
        scores = json.load(file)

    for entry in scores:
        dbsnp_id = entry['dbsnp']['rsid']
        impact = check_impacted(entry)
        click.echo(f"{db}\t{dbsnp_id}\t{'⚠️' if impact else ''}")


def check_impacted(entry) -> bool:
    """Check if an entry from MyVariant has been reported as damaged in either PolyPhen2 or SIFT."""
    cadd = entry.get('cadd')
    if cadd is None:
        return False

    polyphen = cadd.get('polyphen', {})
    if isinstance(polyphen, list):
        polyphen = polyphen[0]

    sift = cadd.get('sift', {})
    if isinstance(sift, list):
        sift = sift[0]

    return (
        polyphen.get('cat') in POLYPHEN2_DAMAGING
        or sift.get('cat') in SIFT_DAMAGING
    )


if __name__ == '__main__':
    main()
