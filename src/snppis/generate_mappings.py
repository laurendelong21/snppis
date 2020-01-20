# -*- coding: utf-8 -*-

"""Map pathways to SNPs using dbSNP."""

import json
import os
from collections import defaultdict
from typing import Mapping, Set

import click
import pandas as pd
from compath_utils import CompathManager
from tqdm import tqdm

import bio2bel_kegg
import bio2bel_reactome
import bio2bel_wikipathways
from snppis.constants import MAPPINGS

#import urllib.request
import json
import bz2
import itertools
import time


@click.command()
def main():
    """Generate the mappings file in the resources directory."""

    click.echo('Getting gene to SNP mapping from dbSNP')
    gene_to_snps = make_snp_dict(gene_to_snps=defaultdict(set))

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
        tsv_path = os.path.join(MAPPINGS, f'{manager.module_name}.tsv')
        pathway_to_snp_df.to_csv(tsv_path, sep='\t', index=False)

        click.echo(f'Writing pathway/snp mapping JSON for {manager.module_name}')
        pathway_to_snp_json = get_pathway_to_snp_json(pathway_to_snp, manager.module_name)
        json_path = os.path.join(MAPPINGS, f'{manager.module_name}.json')
        with open(json_path, 'w') as file:
            json.dump(pathway_to_snp_json, file, indent=2)


def make_snp_dict(gene_to_snps):
    """Get a mapping from gene to set of dbSNP identifiers."""
    chroms = list(range(1, 23)) + ['X', 'Y', 'MT']
    for chrom in chroms:  
        t0 = time.time()
        # Here we begin the downloading of JSON files from the dbSNP database:
        # this is not allowed on the Fraunhofer cluster
    
        #url = 'ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chr{}.json.bz2'.format(chrom)
        path = '/home/llong/Downloads/refsnp/refsnp-chr{}.json.bz2'.format(chrom)
        #if not os.path.exists(path):
            #print('Beginning file download of chromosome {} with urllib2...'.format(chrom))
            #urllib.request.urlretrieve(url, path)
            #print('...Finished file download of chromosome {} with urllib2.'.format(chrom))
    
    
        # Here we parse through the files:
        print('Now decompressing and reading JSON.bz2 files from chromosome {} with *bz2* and *json* ...'.format(chrom))
        with bz2.BZ2File(path, 'rb') as f_in:
            for line in f_in:
                rs_obj = json.loads(line.decode('utf-8'))
                dbsnp_id = 'rs' + rs_obj['refsnp_id']  # the dbsnp id
    
                all_ann_list_raw = rs_obj['primary_snapshot_data'][
                    'allele_annotations']  # these are the assembly annotations
    
                if len(all_ann_list_raw) >= 2:  # if it has sufficient info
                    assembl_ann_list_raw = all_ann_list_raw[1]['assembly_annotation']  # against each assembly
                    if len(assembl_ann_list_raw) != 0:  # if it contains gene info
                        gene_list_raw = assembl_ann_list_raw[0][
                            'genes']  # and each of the genes affected within each assembly
                        if len(gene_list_raw) > 0:
                            # Here I start extracting gene info:
                            for x, y, z in itertools.product(range(len(all_ann_list_raw)),
                                                             range(len(assembl_ann_list_raw)),
                                                             range(len(gene_list_raw))):
                                symbol = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['locus']
    
                                gene_to_snps[symbol].add(dbsnp_id)
        t1 = time.time()
        totaltime = t1 - t0
        print("Finished writing gene:snp pairs from chromosome {} to the dictionary in {} seconds.".format(chrom, totaltime))
   
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
