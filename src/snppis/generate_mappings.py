# -*- coding: utf-8 -*-

"""Map pathways to SNPs using PheWAS Catalog."""

import json
import os
from collections import defaultdict
from typing import Mapping, Set

import click
import pandas as pd
from compath_utils import CompathManager
from tqdm import tqdm

import bio2bel_kegg
import bio2bel_phewascatalog
import bio2bel_reactome
import bio2bel_wikipathways
from snppis.constants import RESOURCES


@click.command()
def main():
    """Generate the mappings file in the resources directory."""
    click.echo('Getting PheWAS Catalog data')
    phewascatalog_df = get_phewascatalog_df()

    click.echo('Getting gene to SNP mapping from PheWAS Catalog')
    gene_to_snps = get_gene_to_snps(phewascatalog_df)

    managers = [
        bio2bel_kegg.Manager(),
        bio2bel_wikipathways.Manager(),
        bio2bel_reactome.Manager(),
    ]

    for manager in managers:
        click.echo(f'Generating pathway/gene mapping for {manager.module_name}')
        pathway_to_gene = get_pathway_to_gene(manager)

        click.echo(f'Generating pathway/snp mapping for {manager.module_name}')
        pathway_to_snp = get_pathway_to_snp(pathway_to_gene, gene_to_snps)

        click.echo(f'Writing pathway/snp mapping TSV for {manager.module_name}')
        pathway_to_snp_df = get_pathway_to_snp_df(pathway_to_snp, manager.module_name)
        tsv_path = os.path.join(RESOURCES, f'{manager.module_name}.tsv')
        pathway_to_snp_df.to_csv(tsv_path, sep='\t', index=False)

        click.echo(f'Writing pathway/snp mapping JSON for {manager.module_name}')
        pathway_to_snp_json = get_pathway_to_snp_json(pathway_to_snp, manager.module_name)
        json_path = os.path.join(RESOURCES, f'{manager.module_name}.json')
        with open(json_path, 'w') as file:
            json.dump(pathway_to_snp_json, file, indent=2)


def get_phewascatalog_df() -> pd.DataFrame:
    """Get the PheWAS Catalog as a pandas DataFrame."""
    phewascatalog_df = bio2bel_phewascatalog.parser.get_df()
    phewascatalog_df = phewascatalog_df[phewascatalog_df.gene_name.notna()]
    return phewascatalog_df


def get_gene_to_snps(phewascatalog_df: pd.DataFrame) -> Mapping[str, Set[str]]:
    """Get a mapping from gene to set of dbSNP identifiers."""
    gene_to_snps = defaultdict(set)
    for snp, gene_symbol in phewascatalog_df[['snp', 'gene_name']].values:
        gene_to_snps[gene_symbol].add(snp)
    return dict(gene_to_snps)


def get_pathway_to_gene(manager: CompathManager):
    """Get a pathway to set of HGNC gene symbol mapping."""
    pathway_to_gene = defaultdict(set)

    for pathway in tqdm(manager.get_all_pathways(), desc='Getting pathways/genes'):
        for protein in pathway.proteins:
            pathway_to_gene[pathway].add(protein.hgnc_symbol)

    return dict(pathway_to_gene)


def get_pathway_to_snp(pathway_to_gene, gene_to_snps):
    """Combine the mappings to relate pathways to SNPs.

    This could be further extended to weight pathway-SNP associations
    by count, or to normalize by the frequency of each SNP being
    mapped to multiple genes.
    """
    pathway_to_snp = defaultdict(set)

    for pathway, gene_symbols in pathway_to_gene.items():
        for gene_symbol in gene_symbols:
            pathway_to_snp[pathway].update(gene_to_snps.get(gene_symbol, set()))

    return pathway_to_snp


def get_pathway_to_snp_df(pathway_to_snp, db) -> pd.DataFrame:
    """Get a data frame with pathway/SNP info."""
    return pd.DataFrame(
        [
            (
                db,
                getattr(pathway, f'{db}_id'),
                pathway.name,
                ','.join(snps),
            )
            for pathway, snps in pathway_to_snp.items()
        ],
        columns=['pathway_namespace', 'pathway_identifier', 'pathway_name', 'snps'],
    )


def get_pathway_to_snp_json(pathway_to_snp, db):
    """Get a JSON object with pathway/SNP info."""
    return [
        {
            'pathway': {
                'namespace': db,
                'identifier': getattr(pathway, f'{db}_id'),
                'name': pathway.name,
            },
            'snps': [
                {
                    'namespace': 'dbsnp',
                    'identifier': snp,
                }
                for snp in snps
            ],
        }
        for pathway, snps in pathway_to_snp.items()
    ]


if __name__ == '__main__':
    main()
