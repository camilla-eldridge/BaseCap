#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: camilla eldridge
"""

''' Trims low qc bases, primers(if found) and output trimmed sequence and its '.qual file '''

import regex
import os
import sys
from typing import List, Tuple, Optional


def reverse_comp(seq: str) -> str:
    """ Return reverse complement """
    result = ""
    seq = seq.lower()
    comp = seq.replace("a", "T").replace("t", "A").replace("c", "G").replace("g", "C").lower()
    result = "".join(comp[::-1])
    return result


def find_variable_seq(sequence: str, primer: str, error_n: int) -> Optional[regex.Match]:
    """ Return variable string match allowing n errors (del &/or subst) """
    p = str("(" + primer.lower() + ")" + "{e<=" + str(error_n) + "}")
    q = regex.search(p, sequence.lower())
    return q


def index_phd(phd: str) -> str:
    """ Index phd file - adds column of base pos """
    indexed = ""
    phd = filter(None, phd.split("\n"))

    for x, value in enumerate(phd, 1):
        indexed = indexed + str(x) + " " + str(value) + "\n"

    return indexed


class Main(object):
    """   trim low quality bases and primers from .Phd file, output qual and trimmed sequence """

    def __init__(self, phd_file: str, ID: str, primers: str, threshold: int) -> None:
        self.phd_file = phd_file
        self.ID = ID
        self.primers = primers
        self.threshold = threshold
        self.sequence: str = ""
        self.trim_pos: str = ""
        self.check: int = 1
        self.trim_start: int = 0
        self.trim_stop: int = 0
        self.start: int = 0
        self.stop: int = 0
        self.primer_pos: Optional[regex.Match] = None
        self.search: Optional[regex.Match] = None
        self.primz: List[str] = []
        self.trimmed_seq: str = ""
        self.trimmed_phd: List[str] = []
        self.seq_cap: str = ""
        self.scores: str = ""
        self.phd2: str = ""

    def phd_to_sequence(self) -> str:
        """ Get sequence from phd_file """
        with open(self.phd_file) as self.phd1:
            self.phd2 = self.phd1.read()
            self.phd = index_phd(self.phd2.split("BEGIN_DNA")[1].split("END_DNA")[0].rstrip())
            self.sequence = self.sequence + "".join(self.phd.split()[1::4])
        return self.sequence

    def trim_pos_ttuner(self) -> Tuple[int, int, int]:
        """ Find trim positions from base caller ttuner """
        self.trim_pos = ""

        for line in self.phd2.split("\n"):
            if "TRIM:" in line:
                self.trim_pos = self.trim_pos + line

        self.trim_pos = self.trim_pos.split()[1:3]

        self.trim_start = int(self.trim_pos[0])
        self.trim_stop = int(self.trim_pos[1])

        if self.trim_stop < 0 and self.trim_start < 0:
            self.check = 0
        else:
            self.check = 1

        return self.trim_start, self.trim_stop, self.check

    def trim_check(self) -> Tuple[int, int]:
        """ See if trimming is recommended by ttuner """
        if self.trim_start < 0:
            self.trim_start = 0
        if self.trim_stop < 0:
            self.trim_stop = int(len(self.sequence))
        return self.trim_start, self.trim_stop

    def reverse_primers(self) -> List[str]:
        """ Reverse complement both primers """
        with open(self.primers) as primerz:
            primer_list = primerz.read().split("\n")[1::2]
            self.primz = primer_list + [reverse_comp(prim_seq) for prim_seq in primer_list]
        return self.primz

    def primer_search(self) -> Optional[regex.Match]:
        """ Search for each primer (forward and reverse) - expecting one hit """
        for p in self.primz:
            self.search = find_variable_seq(self.sequence, p, 3)
            if self.search is not None:
                self.primer_pos = self.search
                return self.primer_pos
        return None

    def locate_primer(self) -> Tuple[int, int]:
        """ Locate which primer hit (assuming orient not known beforehand) """
        midpoint = int(len(self.sequence)) // 2
        if self.primer_pos.span()[0] > midpoint:
            self.primer_start = 0
            self.primer_stop = int(self.primer_pos.span()[0])
        else:
            self.primer_start = int(self.primer_pos.span()[1])
            self.primer_stop = int(len(self.sequence))
        return self.primer_start, self.primer_stop

    def compare_trim(self) -> Tuple[int, int]:
        """ Compare primer and trim positions """
        self.start = max(self.trim_start, self.primer_start)
        self.stop = min(self.trim_stop, self.primer_stop)
        return self.start, self.stop

    def trim_phd_and_seq(self) -> Tuple[str, List[str]]:
        """ Trim phd and sequence """
        self.trimmed_seq = self.sequence[int(self.start):int(self.stop)]
        self.trimmed_phd = self.phd.split("\n")[int(self.start):int(self.stop)]
        return self.trimmed_seq, self.trimmed_phd

    def cap_the_bases(self) -> Tuple[str, str]:
        """ Capitalise the low scoring bases and get scores for qual """
        self.trimmed_phd2 = [ele for ele in self.trimmed_phd if ele.strip()]
        self.seq_cap = ""
        self.scores = ""

        for j in self.trimmed_phd2:
            parts = j.split()
            if len(parts) >= 3:
                score = int(parts[2])
                base = parts[1]
                self.scores += f"{score} "
                self.seq_cap += base.upper() if score < self.threshold else base.lower()
        return self.seq_cap, self.scores

    def write_out(self) -> None:
        """ Write out trimmed and capped sequence """
        with open(f"{self.ID}.qual", "w") as qual_file:
            qual_file.write(f">{self.ID}\n{self.scores}")
        with open(f"{self.ID}.fasta", "w") as seq_file:
            seq_file.write(f">{self.ID}\n{self.seq_cap}")

    def call_methods(self) -> str:
        """ Call initial methods """
        methods1 = ["phd_to_sequence", "trim_pos_ttuner", "reverse_primers", "primer_search", "trim_check"]
        for method in methods1:
            getattr(self, method)()

        if self.primer_pos is None:
            self.start, self.stop = self.trim_start, self.trim_stop
        else:
            self.locate_primer()
            self.compare_trim()

        if self.check == 0:
            self.trimmed_seq = self.sequence
            self.trimmed_phd = self.phd.split("\n")
        else:
            self.trim_phd_and_seq()

        self.cap_the_bases()
        self.write_out()

        search_output = self.search.group() if self.search else "None"
        trimmed_length = len(self.trimmed_seq) if self.check else len(self.sequence)

        summary_vars = [
            self.ID, self.check, self.start, self.stop,
            len(self.sequence), trimmed_length,
            self.trim_start, self.trim_stop, search_output
        ]
        summary = ",".join(map(str, summary_vars))
        return summary


if __name__ == "__main__":
    Directory: str = sys.argv[1]
    Primers: str = sys.argv[2]
    Threshold: int = int(sys.argv[3])

    header_list: List[str] = ["ID", "trimmed(0/1)", "start", "stop", "read_len", "trimmed_read_len", "trim_start", "trim_stop", "primer_hit"]

    print(",".join(header_list))

    for filename in os.listdir(Directory):
        f = os.path.join(Directory, filename)
        ID = str(filename.split(".")[0])

        if filename.endswith(".phd.1"):
            Z = Main(f, ID, Primers, Threshold)
            print(Z.call_methods())
