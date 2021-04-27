import json
import sys
from typing import List

import click as click

from marker_search import blast_utils, match_utils

porA_length = 1188


# 'Present' if at least one non-disrupted copy is present.
def determine_status(matches: List):
    if len(matches) == 0:
        return 'Not found'
    status = 'Present'
    for match in matches:
        if match.isDisrupted:
            status = 'Pseudogene'
        elif match.coverage > 0.8:
            return 'Present'
    return status


@click.command()
@click.option('--fasta', type=click.Path(exists=True), help='Input FASTA file path')
@click.option('--blast_db', type=click.Path(exists=True), help='BLAST DB directory')
def run_search(fasta: str, blast_db):
    evalue = 1e-20
    matches = blast_utils.run_blast(query=fasta, blast_db=blast_db + "/markers", evalue=evalue)
    # All porA at the moment.
    classified_hits = match_utils.classify_matches(matches['porA'], porA_length)

    result = dict()
    result['status'] = determine_status(classified_hits)
    result['matches'] = classified_hits

    print(json.dumps(result, default=lambda x: x.__dict__), file=sys.stdout)


if __name__ == '__main__':
    import logging

    logging.basicConfig()
    run_search()
