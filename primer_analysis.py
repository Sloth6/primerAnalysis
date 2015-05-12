#!/usr/bin/python
"""
post_merged takes a fastq file and set of forward/reverse primers and
calculates the distributions of each primer. It can also produce a
.has_primers.fastq file of all the sequence with the forward and reverse primer.
"""
import re, json, sys
# from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Alphabet import generic_dna#, generic_protein
from Bio.Seq import Seq

class Primer(object):
    """A single primer sequence"""
    def __init__(self, seq, index, reverse=False):
        if reverse:
            self.index = 1
            self.seq = str(Seq(seq, generic_dna).reverse_complement())
        else:
            self.seq = seq
        self.frame = (index-1) % 3
        self.stripped = seq[3 - self.frame:]

    def __repr__(self):
        return self.seq

class PrimerGroup(object):
    """
    A collection of forward and reverse primers. A sequence 'matches' a
    primer group if it contains any one of the forward primers and any one
    of the reverse primers.
    """
    def __init__(self, primers):
        self.forwards = self.make_all_primers(primers["forwards"])
        self.reverses = [Primer(primers["reverses"][0],1, True)]
        #self.make_all_primers(primers["reverses"], True)
        self.regex = self.make_regex()

    def make_all_primers(self, raw_primers, reverse=False):
        """Expands all degenerate primers"""
        nt_map = {
            'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'],
            'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T']
        }
        primers = []
        # primer_map = dict()
        while len(raw_primers) > 0:
            had_degen = False
            raw_primer = raw_primers.pop(0)
            for nt in raw_primer[0]:
                if nt in nt_map:
                    had_degen = True
                    for map_to in nt_map[nt]:
                        raw_primer[0] = raw_primer[0].replace(nt, map_to, 1)
                        raw_primers.append(raw_primer)
                    break
            if not had_degen:
                primers.append(Primer(raw_primer[0], raw_primer[1], reverse))
        return primers

    def make_regex(self):
        """
        Create a regular expression that will match a sequence if it has one of
        the forward primers before one of the reverse primerse somewhere in the
        sequence.
        """
        forwards = self.forwards
        re_string = '^.*((' + (')|('.join([f.seq for f in forwards])) + \
            ')).*'+ str(self.reverses[0])+'.*$'
        return re.compile(re_string)

    def forward_primer(self, seq):
        """
        Confirm the sequence matched againsed our regular expression. If it
        does then return the forward primer that matched to it.
        """
        match = self.regex.match(seq)
        if match:
            return match.group(1)
        else:
            return None

# def encode_protein(seq, forward_primer, reverse_primer):
#     """
#     Translate the sequence to protein. Remove segments before and after
#     the primers. 
#     """
#     if forward_primer is None:
#         return None

#     start = seq.index(forward_primer)
#     end = seq.index(reverse_primer)+len(reverse_primer)
#     subseq = seq[start: end]
#     protein = Seq(subseq, generic_dna).translate(to_stop=True)
#     return protein
    # if len(protein) >= 100:
    #     return protein
    # else:
    #     return None
    # print protein

def process_seq(seq, primerGroup, primerCount):
    forward_primer = primerGroup.forward_primer(seq)
    # protein = encode_protein(seq, forward_primer, primerGroup.reverses[0])
    if forward_primer:# and protein:
        if forward_primer in primerCount:
            primerCount[forward_primer] += 1
        else:
            primerCount[forward_primer] = 1

def calculate_primer_distribution(path, primerGroup):
    """
    Iterate all sequences and see if the match to the primer group.
    If they do, record which forward primer.
    """
    primerCount = dict(zip(primerGroup.forwards, [0]*len(primerGroup.forwards)))
    for _, seq, _ in FastqGeneralIterator(open(path)):
        one_seq(seq, primerGroup, primerCount)
    return primerCount

def main(argv):
    """
    Load the primer group and ensure the passed argv param matches any one 
    from file.
    """
    [path, primerGroup] = argv
    with open('primers.json') as all_primers_file:
        all_primers = json.load(all_primers_file)
        for group in all_primers['primer_groups']:
            if group["name"] == primerGroup:
                return calculate_primer_distribution(path, PrimerGroup(group))
        return "Please select one of the primer groups defined in primers.json"

if __name__ == "__main__":
    print main(sys.argv[1:])
