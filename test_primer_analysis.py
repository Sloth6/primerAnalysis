import json
import primer_analysis as pa

def test_process_seq():
	with open('./test_primers.json') as primers_file:
		primers_json = json.load(primers_file)["primer_groups"]
		primer_groups = [pa.PrimerGroup(group) for group in primers_json]

		noise = ''#"ATCGATCGTTGAC"
		test_sequences = [
			noise+noise, #no match
			noise+primers_json[0]["reverses"][0]+noise, #reverse, no forward
			noise+primers_json[0]["forwards"][0]+noise, #forward, no reverse

			noise+primers_json[0]["forwards"][0]+noise+"CCGGGAA"+noise,
			noise+primers_json[0]["forwards"][1]+noise+"CCGGGAA"+noise,
			noise+primers_json[0]["forwards"][2]+noise+"CCGGGAA"+noise,

			# Test degenerate primer for all combinations. Should match each one. 
			noise+"ACAACA"+noise+"CCGGGAA"+noise,
			noise+"ACAATA"+noise+"CCGGGAA"+noise,
			noise+"ATAACA"+noise+"CCGGGAA"+noise,
			noise+"ATAATA"+noise+"CCGGGAA"+noise,
		]
		
		total = 0
		matches = 0
		for seq in test_sequences:
			total += 1
			if pa.process_seq(seq, primer_groups):
				matches += 1
		pa.output(primer_groups, total, matches)
		# assert(primer_groups[0].counts[forwards[0]] == 1)
		# assert(primer_groups[0].counts[forwards[1]] == 1)
		# assert(primer_groups[0].counts[forwards[2]] == 1)
		# assert(primer_groups[0].counts[forwards[3]] == 0)

if __name__ == "__main__":
    test_process_seq()