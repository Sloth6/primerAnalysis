import json
import primer_analysis

forwards = ["AAAAAA", "TTTTTT", "GGGGGG"]
reverse = "CCGGGAA"
noise = "ATCGATCG"

sequences = [
	"ATCGATCGATCG", #no match
	noise+reverse+noise, #reverse, no forward
	noise+forwards[0]+noise, #forward, no reverse
	noise+forwards[0]+noise+reverse+noise, #should match AAAAAA
	noise+forwards[1]+noise+reverse+noise, #should match TTTTTT
	noise+forwards[2]+noise+reverse+noise, #should match GGGGGG
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