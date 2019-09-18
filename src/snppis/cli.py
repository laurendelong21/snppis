# -*- coding: utf-8 -*-

"""Command line interface for :mod:`snppis`."""

import logging
import sys
from typing import TextIO

import click
import pandas as pd

from .pipeline import get_pathway_to_patient_to_score

logger = logging.getLogger(__name__)


@click.command()
@click.option('-f', '--file', type=click.File(), required=True)
@click.option('-o', '--output', type=click.File('w'), required=True)
@click.option('-s', '--sep', default=',')
@click.option('-v', '--debug', is_flag=True)
def main(file: TextIO, output: TextIO, sep: str, debug: bool):
    """Run the pipeline on a given file."""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=level)
    logging.getLogger('snppis').setLevel(level)

    click.secho(f'Loading data from {file.name}', fg='blue', err=True)
    df = pd.read_csv(file, sep=sep)

    click.secho('Getting patient/pathway scores', fg='blue', err=True)
    pathway_df = get_pathway_to_patient_to_score(df)

    pathway_df.to_csv(output, sep=sep, index=False)


if __name__ == '__main__':
    main()
