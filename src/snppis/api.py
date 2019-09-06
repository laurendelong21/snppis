# -*- coding: utf-8 -*-

"""Acquisition and exploration of SNP-related Bio2BEL packages."""

from typing import List

import myvariant
import pandas as pd

mv = myvariant.MyVariantInfo()


def get_polyphen_scores(dbsnp_ids: List[str]) -> pd.DataFrame:
    """Get polyphen2 scores for a list of SNPs."""
    query = ' or '.join(
        f'dbsnp.rsid:{dbsnp_id}'
        for dbsnp_id in dbsnp_ids
    )
    df = mv.query(
        query,
        fields='dbsnp.rsid, dbnsfp.polyphen2.hdiv',
        as_dataframe=True,
    )
    df = df[['dbsnp.rsid', 'dbnsfp.polyphen2.hdiv.pred']].drop_duplicates()
    df.rename(columns={
        'dbsnp.rsid': 'dbsnp_id',
        'dbnsfp.polyphen2.hdiv.pred': 'polyphen_prediction',
    }, inplace=True)
    return df
