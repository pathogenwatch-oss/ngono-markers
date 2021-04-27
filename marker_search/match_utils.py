import re
from typing import List
from Bio.Seq import Seq

indel_finder = re.compile(r'[-]+')


class Match:
    def __init__(self, contig_id: str, query_start: int, query_end: int, ref_start: int, ref_end: int, frame: int,
                 is_forward: bool, coverage: float, is_disrupted: bool, is_exact: bool, identity: float):
        self.contigId = contig_id
        self.queryStart = query_start
        self.queryEnd = query_end
        self.refStart = ref_start
        self.refEnd = ref_end
        self.frame = frame
        self.isForward = is_forward
        self.coverage = coverage
        self.isDisrupted = is_disrupted
        self.isExact = is_exact
        self.identity = identity


def at_least_one_exact_match(hits):
    return 0 < len(only_exact_matches(hits))


def classify_matches(family_matches: dict, ref_length: int) -> List[Match]:
    matches = list()
    for contig_id, hsps in family_matches.items():
        for hsp in hsps:
            ref_start, ref_end = sorted([hsp.sbjct_start, hsp.sbjct_end])
            is_forward = 0 < hsp.frame[1]
            coverage = float(ref_end - ref_start + 1) / float(ref_length)
            is_exact = coverage == 1.0 and hsp.query == hsp.sbjct
            query_seq = hsp.query if is_forward else str(Seq(hsp.query).reverse_complement())
            is_disrupted = find_frameshift(query_seq, hsp.sbjct) or find_premature_stop(query_seq, hsp.frame[0],
                                                                                        ref_end < ref_length)
            percent_identity = round((float(hsp.identities) / float(ref_length)) * 100, 2)
            matches.append(
                Match(contig_id, hsp.query_start, hsp.query_end, ref_start, ref_end, hsp.frame[0],
                      is_forward, coverage, is_disrupted, is_exact, percent_identity))
    return matches


def find_frameshift(query: str, sbjct: str) -> bool:
    insert_matches = [insert for insert in indel_finder.findall(query) if len(insert) % 3 != 0]
    if len(insert_matches) != 0:
        return True
    deletion_matches = [deletion for deletion in indel_finder.findall(sbjct) if len(deletion) % 3 != 0]
    if len(deletion_matches) != 0:
        return True


def find_premature_stop(dna: str, frame: int, includes_end: bool) -> bool:
    dna = dna.replace('-', '')[frame - 1:]
    end_offset = (len(dna) % 3) * -1
    if end_offset != 0:
        includes_end = False
    dna = dna[0:end_offset]
    coding_seq = Seq(dna)
    translation = coding_seq.translate()
    if translation.count('*') < 1:
        return False
    else:
        end_check = 1 if includes_end else 0
        return str(translation).index('*') < len(translation) - end_check


def only_exact_matches(hits):
    return [hit for hit in hits if hit.isExact]
