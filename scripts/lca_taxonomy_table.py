#!/usr/bin/env python3

from argparse import ArgumentParser
from collections import defaultdict

def parse_options():
	parser = ArgumentParser(description='Resolve taxonomy with LCA for multiple hits')
	parser.add_argument('-i', '--input', required=True, help='path to vsearch input file')
	parser.add_argument('-o', '--output', required=True, help='name of output file')
	parser.add_argument('--lca-threshold', type=float, default=0.0,
	                    help='Only apply LCA if top hits are within this percent identity difference (default: 0.0 = exact ties only)')
	return parser.parse_args()


def importing(file):
	"""Import VSEARCH blast6out format file"""
	input_list = []
	with open(file) as vsearch:
		for line in vsearch.readlines():
			if line.strip():  # Skip empty lines
				input_list.append(line.strip().split('\t'))
	return input_list


def parse_taxonomy(tax_string):
	"""Parse PR2 taxonomy format: tax=d:Domain,k:Kingdom,p:Phylum,..."""
	tax_string = tax_string.replace('tax=', '')
	ranks = ['Domain', 'Supergroup', 'Division', 'Subdivision', 'Class', 'Order', 'Family', 'Genus', 'Species']

	# Parse the taxonomy parts
	parts = tax_string.split(',')
	tax_dict = {}

	for part in parts:
		if ':' in part:
			rank_code, value = part.split(':', 1)
			tax_dict[rank_code] = value if value and value != '-' else None

	# Map to standard rank names (PR2 format uses d, k, p, c, o, f, g, s or similar)
	# Adjust based on actual number of taxonomic levels
	taxonomy = []
	for part in parts:
		if ':' in part:
			value = part.split(':', 1)[1]
			taxonomy.append(value if value and value != '-' else None)

	# Pad or trim to match expected number of ranks
	while len(taxonomy) < len(ranks):
		taxonomy.append(None)

	return taxonomy[:len(ranks)]


def group_hits_by_query(imported_list):
	"""Group all hits by query sequence ID"""
	hits_by_query = defaultdict(list)

	for rec in imported_list:
		if len(rec) >= 4:
			query_id = rec[0]
			subject_id = rec[1]
			pident = float(rec[2])
			length = rec[3]

			# Parse subject ID and taxonomy (format: seqid;tax=...)
			if ';' in subject_id:
				seqacc, seqtax = subject_id.split(';', 1)
			else:
				seqacc = subject_id
				seqtax = 'tax=' + ','.join(['-'] * 9)  # Unknown taxonomy

			taxonomy = parse_taxonomy(seqtax)

			hits_by_query[query_id].append({
				'seqacc': seqacc,
				'taxonomy': taxonomy,
				'pident': pident,
				'length': length
			})

	return hits_by_query


def find_lca(taxonomies):
	"""Find the lowest common ancestor from a list of taxonomies"""
	if not taxonomies:
		return [None] * 9

	if len(taxonomies) == 1:
		return taxonomies[0]

	# Find the deepest level where all taxonomies agree
	lca = []
	num_ranks = len(taxonomies[0])

	for rank_idx in range(num_ranks):
		# Get all values at this rank
		values_at_rank = set()
		for tax in taxonomies:
			if tax[rank_idx] is not None:
				values_at_rank.add(tax[rank_idx])

		# If all agree (only one unique value), keep it
		if len(values_at_rank) == 1:
			lca.append(list(values_at_rank)[0])
		# If there's disagreement, use the parent level and stop
		else:
			# Mark this and all deeper levels as unresolved
			lca.append(None)
			for remaining in range(rank_idx + 1, num_ranks):
				lca.append(None)
			break

	return lca


def resolve_taxonomy_lca(hits_by_query, lca_threshold):
	"""Resolve taxonomy using LCA for ties"""
	resolved = []

	for query_id, hits in hits_by_query.items():
		if not hits:
			continue

		# Sort by percent identity (descending)
		hits.sort(key=lambda x: x['pident'], reverse=True)

		# Get the best percent identity
		best_pident = hits[0]['pident']

		# Find all hits within the LCA threshold
		top_hits = [h for h in hits if h['pident'] >= best_pident - lca_threshold]

		# If only one top hit, use it directly
		if len(top_hits) == 1:
			hit = top_hits[0]
			resolved.append({
				'qseqid': query_id,
				'seqacc': hit['seqacc'],
				'taxonomy': hit['taxonomy'],
				'pident': hit['pident'],
				'length': hit['length'],
				'lca_applied': False,
				'num_hits': 1
			})
		else:
			# Multiple top hits - apply LCA
			taxonomies = [h['taxonomy'] for h in top_hits]
			lca_taxonomy = find_lca(taxonomies)

			# Use the first hit's metadata but with LCA taxonomy
			resolved.append({
				'qseqid': query_id,
				'seqacc': f"LCA_{len(top_hits)}_hits",
				'taxonomy': lca_taxonomy,
				'pident': best_pident,
				'length': top_hits[0]['length'],
				'lca_applied': True,
				'num_hits': len(top_hits)
			})

	return resolved


def format_taxonomy(taxonomy):
	"""Convert taxonomy list to tab-separated string, replacing None with '-'"""
	return '\t'.join([t if t is not None else '-' for t in taxonomy])


def saving_output(resolved_data, output_file):
	"""Save the resolved taxonomy table"""
	with open(output_file, 'w') as output:
		# Write header
		header = 'OTU\tPident\tAccession\tDomain\tSupergroup\tDivision\tSubdivision\tClass\tOrder\tFamily\tGenus\tSpecies\tLength\tLCA_Applied\tNum_Hits\n'
		output.write(header)

		# Write data
		for entry in resolved_data:
			line = (
				f"{entry['qseqid']}\t"
				f"{entry['pident']:.1f}\t"
				f"{entry['seqacc']}\t"
				f"{format_taxonomy(entry['taxonomy'])}\t"
				f"{entry['length']}\t"
				f"{'Yes' if entry['lca_applied'] else 'No'}\t"
				f"{entry['num_hits']}\n"
			)
			output.write(line)


def main():
	options = parse_options()

	# Import VSEARCH output
	imported_file = importing(options.input)

	# Group hits by query
	hits_by_query = group_hits_by_query(imported_file)

	# Resolve taxonomy using LCA
	resolved_data = resolve_taxonomy_lca(hits_by_query, options.lca_threshold)

	# Save output
	saving_output(resolved_data, options.output)

	# Print summary statistics
	total_queries = len(resolved_data)
	lca_applied = sum(1 for entry in resolved_data if entry['lca_applied'])

	print(f"Taxonomy resolution complete:")
	print(f"  Total OTUs: {total_queries}")
	print(f"  OTUs with unique best hit: {total_queries - lca_applied}")
	print(f"  OTUs resolved with LCA: {lca_applied} ({100*lca_applied/total_queries:.1f}%)")


if __name__ == '__main__':
	main()
