#!/usr/bin/python
"""docstring for post_merged"""
import re, json, sys
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Alphabet import generic_dna#, generic_protein
from Bio.Seq import Seq

class Primer(object):
    """docstring for Primer"""
    def __init__(self, seq, index, reverse=False):
        if reverse:
            self.seq = str(Seq(seq, generic_dna).reverse_complement())
        else:
            self.seq = seq
        self.frame = (index-1) % 3
        self.stripped = seq[3 - self.frame:]
    
    def __repr__(self):
        return self.seq

class PrimerGroup(object):
    """docstring for PrimerGroup"""
    def __init__(self, primers):
        self.forwards = self.make_all_primers(primers["forwards"])
        self.reverses = map(lambda s:
            str(Seq(s, generic_dna).reverse_complement()),
            primers["reverses"]
        )
        self.regex = self.makeRegex()

    def make_all_primers(self, raw_primers, reverse=False):
        """Expands all degenerate primers"""
        nt_map = {
            'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'],
            'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T']
        }
        # open_seqs = map(lambda f: f[0][3 - (f[1]-1)%3:], raw_primers)
        primers = []
        # primer_map = dict()
        while len(raw_primers) > 0:
            had_degen = False
            raw_primer = raw_primers.pop(0)
            for c in raw_primer[0]:
                if c in nt_map:
                    had_degen = True
                    for map_to in nt_map[c]:
                        raw_primer[0] = raw_primer[0].replace(c, map_to, 1)
                        raw_primers.append(raw_primer)
                    break
            if not had_degen:
                primers.append(Primer(raw_primer[0], raw_primer[1], reverse))
        return primers

    def makeRegex(self):
        """
        Create a regular expression that will match a sequence if it has one of
        the forward primers before one of the reverse primerse somewhere in the
        sequence.
        """
        forwards = self.forwards
        re_string = '^.*(('+(')|('.join(map(lambda f: f.seq, forwards)))+')).*'+str(self.reverses[0])+'.*$'
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

def encode_protein(seq, forward_primer, reverse_primer):
    if forward_primer is None:
        return None

    start = seq.index(forward_primer)
    end = seq.index(reverse_primer)+len(reverse_primer)
    subseq = seq[start: end]
    protein = Seq(subseq, generic_dna).translate(to_stop=True)
    return protein
    # if len(protein) >= 100:
    #     return protein
    # else:
    #     return None
    # print protein


def calculate_primer_distribution(path, primerGroup):
    primer_count = dict(zip(primerGroup.forwards, [0]*len(primerGroup.forwards)))

    for title, seq, qual in FastqGeneralIterator(open(path)):
        forward_primer = primerGroup.forward_primer(seq)
        protein = encode_protein(seq, forward_primer, primerGroup.reverses[0])

        if forward_primer and protein:
            if forward_primer in primer_count:
                primer_count[forward_primer] += 1
            else:
                primer_count[forward_primer] = 1

    print primer_count

def main(argv):
    [ path, primerGroup ] = argv
    with open('primers.json') as all_primers_file:
        all_primers = json.load(all_primers_file)
        for group in all_primers['primer_groups']:
            if group["name"] == primerGroup:
                return calculate_primer_distribution(path, PrimerGroup(group))
        print "Please select one of the primer groups defined in primers.json"

if __name__ == "__main__":
   main(sys.argv[1:])

# path_prefix = '/Users/joelsimon/Documents/Immufind/SequenceAnalysis/data/read_1/FlashMerge/'
# main(path_prefix+'2.extendedFrags.fastq', 'igl')
# main(path_prefix+'2.extendedFrags.fastq', igl_primers)