import json
import primer_analysis

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
		primerJSON = json.load(primers_file)["primer_groups"][0]
		primerGroup = primer_analysis.PrimerGroup(primerJSON)
		primerCount = dict()
		for seq in sequences:
			primer_analysis.process_seq(seq, primerGroup, primerCount)
		print primerCount
		assert(primerCount[forwards[0]] == 1)
		assert(primerCount[forwards[1]] == 1)
		assert(primerCount[forwards[2]] == 1)
		# assert(primerCount[forwards[3]] == 4)

if __name__ == "__main__":
    test_process_seq()