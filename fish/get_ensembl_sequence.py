import requests, sys, re, difflib
from biomart import BiomartServer

# can update these if not working (biomart mirror may be down or change)
server = 'https://rest.ensembl.org'
biomart_server = BiomartServer( "http://useast.ensembl.org/biomart" )
biomart_genes = biomart_server.datasets['hsapiens_gene_ensembl']

#python /home/data/Shared/shared_datasets/sm_fish/code/get_ensembl_sequence.py

def get_json_response(server, ext):
	headers = { "Content-Type" : "application/json"}
	r = requests.get(server+ext, headers=headers)
	if not r.ok:
		r.raise_for_status()
		sys.exit()
	return r.json()

def get_all_common_substrings(strings):
	if not strings:
		return set()
	
	# Sort by length to reduce search space
	strings = sorted(strings, key=len)
	shortest = strings[0]
	n = len(shortest)
	common = set()
	
	# Try all substrings of the shortest string
	for i in range(n):
		for j in range(i + 1, n + 1):
			substr = shortest[i:j]
			if all(substr in s for s in strings[1:]):
				common.add(substr)
	
	# Filter to only maximal substrings (not fully contained in any longer common substrings)
	maximal_common = set()
	for s in common:
		if not any((s != other and s in other) for other in common):
			maximal_common.add(s)
	
	return maximal_common


def main():
	
	gene = input("Please enter a gene name (example NR2F2): ")
	gene_info = get_json_response(server, f'/lookup/symbol/homo_sapiens/{gene}?expand=1;mane=1')
	
	gene_id = gene_info['id']
	canonical_transcript_id = gene_info['canonical_transcript'].split('.')[0]
	print('Successfully found gene ID:', gene_id)
	
	print()
	
	print("Please select which type of sequence to query:")
	print("1. Genomic (Includes entire sequence including introns)")
	print("2. cDNA (Includes exons and UTRs)")
	print("3. CDS (Includes only exons)")
	
	while True:
		choice = input('')
		if choice == '1':
			seq_type = 'genomic'
			break
		elif choice == '2':
			seq_type = 'cdna'
			break
		elif choice == '3':
			seq_type = 'cds'
			break
		else:
			print("Invalid choice. Please try again.")
	
	transcripts = get_json_response(server, f'/sequence/id/{gene_id}?type={seq_type};multiple_sequences=1')
	
	# Get additional data on transcripts
	response = biomart_genes.search({
			'attributes': ['ensembl_transcript_id', 'transcript_tsl', 'transcript_gencode_basic'],
			'filters': {'ensembl_transcript_id': [transcript['id'] for transcript in transcripts]}
		}, 
		header = 1 
	)
	
	biomart_list = [item.split('\t') for item in response.text.split('\n')]
	biomart_dict = {}
	for item in biomart_list:
		if len(item) < 2:
			continue
		match = re.search('(?<=tsl)\d+', item[1])
		if match:
			biomart_dict[item[0]] = {'tsl': int(match.group(0)), 'gencode_basic': item[2]}
			
	print("Please select which transcripts to include:")
	all_sequences = []
	for transcript in transcripts:
		
		id = transcript['id']
		seq_length = len(transcript['seq'])
		tsl = biomart_dict[id]['tsl']
		
		if biomart_dict[id]['gencode_basic'] == '':
			protein_coding = False
		else:
			protein_coding = True
		
		if id == canonical_transcript_id:
			is_canonical = True
		else:
			is_canonical = False
		
		print(f'Transcript ID: {id}, Canonical Transcript: {is_canonical}, Length: {seq_length}, TSL: {tsl}, Protein Coding: {protein_coding}')
		print('Include? (y/n)')
		
		while True:
			choice = input('')
			if choice == 'y':
				print(f'Selected {id}')
				all_sequences.append(transcript['seq'])
				break
			elif choice == 'n':
				break
			else:
				print("Invalid choice. Please try again.")
			
	print()
	
	if len(all_sequences) > 1:
		
		print('Finding common contiguous sequence(s)...')
		
		all_common_substrings = get_all_common_substrings(all_sequences)
		
		for substring in all_common_substrings:
			# contiguous sequences that are shorter than a probe's length are not relevant
			if len(substring) > 10:
				print('Length:', len(substring))
				print(substring)
	else:
		print('Length:', len(all_sequences[0]))
		print(all_sequences[0])
	
if __name__=="__main__":
	main()
