# -*- coding: utf-8 -*-

"""Command line interface for ``snppis``."""

import click

from snppis.api import get_polyphen_scores


@click.command()
@click.argument('dbsnp_ids')
def main(dbsnp_ids: str):
    """Get polyphen2 scores."""
    dbsnp_ids = dbsnp_ids.split(',')
    scores = get_polyphen_scores(dbsnp_ids)
    click.echo(scores)


if __name__ == '__main__':
    main()
