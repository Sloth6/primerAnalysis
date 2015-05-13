import json
import primer_analysis as pa

forwards = ["AAAAAA", "TTTTTT", "GGGGGG", "AYAAYA"]
reverse = "CCGGGAA"
noise = "ATCGATCG"

sequences = [
	noise+noise, #no match
	noise+reverse+noise, #reverse, no forward
	noise+forwards[0]+noise, #forward, no reverse
	noise+forwards[0]+noise+reverse+noise, #should match AAAAAA
	noise+forwards[1]+noise+reverse+noise, #should match TTTTTT
	noise+forwards[2]+noise+reverse+noise, #should match GGGGGG
	
	# Test degenerate primer for all combinations. Should match each one. 
	noise+"ACAACA"+noise+reverse+noise,
	noise+"ACAATA"+noise+reverse+noise,
	noise+"ATAACA"+noise+reverse+noise,
	noise+"ATAATA"+noise+reverse+noise,
]

def test_process_seq():
	with open('test_primers.json') as primers_file:
		primers_json = json.load(primers_file)["primer_groups"]
		primer_groups = [pa.PrimerGroup(group) for group in primers_json]
		primer_count = dict()
		total = 0
		for seq in sequences:
			total += 1
			pa.process_seq(seq, primer_groups, primer_count)
		print pa.output(primer_count, total, total-1)
		# assert(primer_count[forwards[0]] == 1)
		# assert(primer_count[forwards[1]] == 1)
		# assert(primer_count[forwards[2]] == 1)
		# assert(primer_count[forwards[3]] == 4)

if __name__ == "__main__":
    test_process_seq()