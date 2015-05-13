#!/usr/bin/python
"""
post_merged takes a fastq file and set of forward/reverse primers and
calculates the distributions of each primer. It can also produce a
.has_primers.fastq file of all the sequence with the forward and reverse primer.
"""
import re, json, sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

class PrimerGroup(object):
    """
    A collection of forward and reverse primers. A sequence 'matches' a
    primer group if it contains any one of the forward primers and any one
    of the reverse primers.
    """
    def __init__(self, json_data):
        self.forwards = self.make_all_primers(json_data["forwards"])
        self.reverses = self.make_all_primers(json_data["reverses"], True)
        self.name = json_data["name"]
        self.regex = self.make_regex()

    def make_all_primers(self, raw_primers, reverse=False):
        """Expands all degenerate primers"""
        nt_map = {
            'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'],
            'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T']
        }
        primers = []
        primer_map = dict()
        while len(raw_primers) > 0:
            had_degen = False
            raw_primer = raw_primers.pop(0)
            for nt in raw_primer:
                if nt in nt_map:
                    had_degen = True
                    for map_to in nt_map[nt]:
                        raw_primers.append(raw_primer.replace(nt, map_to, 1))
                    break
            if not had_degen:
                if reverse:
                    raw_primer = str(Seq(raw_primer, generic_dna).reverse_complement())
                primers.append(raw_primer)
        return primers

    def make_regex(self):
        """
        Create a regular expression that will match a sequence if it has one of
        the forward primers before one of the reverse primerse somewhere in the
        sequence.
        """
        forwards_str = ')|('.join(self.forwards)
        reverses_str = ')|('.join(self.reverses)
        re_str = '^.*((' + forwards_str +')).*((' + reverses_str + ')).*$'
        return re.compile(re_str)

    def match(self, seq):
        """
        Confirm the sequence matched againsed our regular expression. If it
        does then return the forward primer that matched to it.
        """
        match = self.regex.match(seq)
        if match:
            return match.group(1)
        else:
            return False

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

def process_seq(seq, primer_groups, primer_count):
    for p_group in primer_groups:
        match = p_group.match(seq)
        if p_group.name not in primer_count:
            primer_count[p_group.name] = dict()
        if match:
            if match in primer_count[p_group.name]:
                primer_count[p_group.name][match] += 1
            else:
                primer_count[p_group.name][match] = 1
            return True
    return False


def output(primer_count, total, matches):
    print 'Total reads:', total
    print "Reads with primers:", matches, "("+str(matches/float(total))+")"
    print ""
    for name, values in primer_count.iteritems():
        print name
        for forward, count in values.iteritems():
            print '\t', forward, ':', count
        print ''
    print "#"*80

def main(argv):
    """
    Load the primer group and ensure that every passed primer groups match
    one from file.
    """
    if len(argv) != 2:
        print "Usage: primer_analysis.py <primerPath.json> <myRead.fastq>"
        return

    [path, primers_path] = argv
    print "#"*80
    print "Primer analysis version 0.0.1"
    print "\t", 'File:', path
    print "\t", 'Primers:', primers_path
    with open(primers_path) as all_primers_file:
        try:
            primers_json = json.load(all_primers_file)['primer_groups']
        except Exception, e:
            print "Failed to parse JSON file. Ensure it is formatted correctly."
            return
        
        primer_groups = [PrimerGroup(group) for group in primers_json]
        primer_count = dict()
        # Iterate all sequences and see if the match to the primer group.
        # If they do, record which forward primer.
        total = 0
        matches = 0
        for _, seq, _ in FastqGeneralIterator(open(path)):
            total += 1
            if process_seq(seq, primer_groups, primer_count):
                matches += 1
            
        return output(primer_count, total, matches)

if __name__ == "__main__":
    main(sys.argv[1:])
