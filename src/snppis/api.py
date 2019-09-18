# -*- coding: utf-8 -*-

"""Acquisition and exploration of SNP-related Bio2BEL packages."""

import itertools as itt
import json
import os
from typing import Any, Iterable, Mapping

from myvariant import MyVariantInfo

from .constants import MYVARIANT_CACHE

__all__ = [
    'batched_query',
    'query',
]

mv = MyVariantInfo()


def batched_query(dbsnp_ids: Iterable[str], n: int = 10) -> Iterable[Mapping[str, Any]]:
    """Submit a batch query and iterate over the "hits" entries."""
    dbsnp_ids = list(dbsnp_ids)

    uncached_dbsnp_ids = set()
    cache = {}
    for dbsnp_id in dbsnp_ids:
        cache_path = os.path.join(MYVARIANT_CACHE, f'{dbsnp_id}.json')
        if not os.path.exists(cache_path):
            uncached_dbsnp_ids.add(dbsnp_id)
            continue
        with open(cache_path) as file:
            cache[dbsnp_id] = json.load(file)

    assert n <= 10
    for dbsnp_ids_chunk in _grouper(n, uncached_dbsnp_ids):
        for result in query(dbsnp_ids_chunk):
            dbsnp_id = result['dbsnp']['rsid']
            with open(os.path.join(MYVARIANT_CACHE, f'{dbsnp_id}.json'), 'w') as file:
                json.dump(result, file, indent=2)
            yield result


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
