import uuid
from collections import defaultdict
from typing import List

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

# Considered present if any match is complete and not disrupted.
from marker_search.match_utils import Match


def determine_status(matches: List[Match]) -> str:
    if len(matches) == 0:
        return "Not found"
    status = "Present"
    for match in matches:
        if match.isDisrupted or match.coverage != 1.0:
            status = "Incomplete"
        else:
            return "Present"
    return status


def overlaps(coords1: tuple, coords2: tuple, threshold=60):
    return min(coords1[1], coords2[1]) - max(coords1[0], coords2[0]) >= threshold


def process_contig(contig_id: str, alignments: list) -> dict:
    excluded = set()
    contig_keep = defaultdict(dict)

    for query in range(0, len(alignments)):
        selected = list()
        query_alignment = alignments[query]

        title = query_alignment.title.split(" ")[0]

        for hsp in query_alignment.hsps:
            name = title + "_" + str(hsp.query_start)
            if name in excluded:
                continue
            for test in range(query, len(alignments)):
                test_ali = alignments[test]
                test_title = test_ali.title.split(" ")[0]
                for test_hsp in test_ali.hsps:
                    test_name = test_title + "_" + str(test_hsp.query_start)
                    if test_name in excluded:
                        continue
                    if name == test_name:
                        continue
                    if overlaps(
                        (hsp.query_start, hsp.query_end),
                        (test_hsp.query_start, test_hsp.query_end),
                    ):
                        if (
                            hsp.align_length - hsp.gaps
                            < test_hsp.align_length - test_hsp.gaps
                        ):
                            excluded.add(test_name)
                            continue
                        elif hsp.identities >= test_hsp.identities:
                            excluded.add(test_name)
                        else:
                            excluded.add(name)
                            break
                    else:
                        continue
            if name not in excluded:
                selected.append(hsp)
        if len(selected) != 0:
            contig_keep[title][contig_id] = selected
    return contig_keep


def remove_overlaps(blast_records) -> dict:
    record_list = list(blast_records)
    kept = dict()

    for contig_search in record_list:
        if len(contig_search.alignments) == 0:
            continue
        kept.update(process_contig(contig_search.query, contig_search.alignments))
    return kept


def run_blast(query, blast_db, evalue):
    out = "/tmp/" + uuid.uuid4().hex
    blastn_cline = NcbiblastnCommandline(
        query=str(query), db=blast_db, evalue=evalue, outfmt=5, out=out
    )
    blastn_cline()
    with open(out, "r") as r_fh:
        selected_records = remove_overlaps(NCBIXML.parse(r_fh))
    # os.remove(out)
    return selected_records
