# -*- coding: utf-8 -*-

"""Acquisition and exploration of SNP-related Bio2BEL packages."""

import itertools as itt
from typing import Any, Iterable, Mapping

from myvariant import MyVariantInfo

__all__ = [
    'batched_query',
    'query',
]

mv = MyVariantInfo()


def batched_query(dbsnp_ids: Iterable[str], n: int = 10) -> Iterable[Mapping[str, Any]]:
    """Submit a batch query and iterate over the "hits" entries."""
    assert n <= 10
    for dbsnp_ids_chunk in _grouper(n, dbsnp_ids):
        yield from query(dbsnp_ids_chunk)


def _grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itt.islice(it, n))
        if not chunk:
            return
        yield chunk


def query(dbsnp_ids: Iterable[str], as_dataframe: bool = False) -> Iterable[Mapping[str, Any]]:
    """Query a set of SNPs by their dbSNP identifiers."""
    q = ' or '.join(
        f'dbsnp.rsid:{dbsnp_id}'
        for dbsnp_id in dbsnp_ids
    )
    result = mv.query(
        q,
        fields='dbsnp.rsid, cadd.polyphen, cadd.sift, wellderly.polyphen',
        as_dataframe=as_dataframe,
    )
    return result['hits']
