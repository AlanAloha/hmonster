#!/usr/bin/env python3
import argparse
from readfasta import read_record

parser = argparse.ArgumentParser(description='Simulate splice-site-generating point mutaitons')
parser.add_argument('fasta', type=str, metavar='<fasta>', help='input path to fasta file')
arg = parser.parse_args()

DON = 'GT'
ACC = 'AG'
for idn, seq in read_record(arg.fasta):
	for i in range(0, len(seq)-1, 2):
		# donor site
			cur_id  = None
			cur_seq = None
			if seq[i] == 'G' and seq[i+1] != 'T':
				cur_id  = f'{idn} {i+2} {seq[i+1]}~T'
				cur_seq = f'{seq[:i+1]}T{seq[i+2:]}'
			if seq[i] != 'G' and seq[i+1] == 'T':
				cur_id  = f'{idn} {i+1} {seq[i]}~G'
				cur_seq = f'{seq[:i]}G{seq[i+1:]}'
		# acceptor site
			if seq[i] == 'A' and seq[i+1] != 'G':
				cur_id  = f'{idn} {i+2} {seq[i+1]}~G'
				cur_seq = f'{seq[:i+1]}G{seq[i+2:]}'
			if seq[i] != 'A' and seq[i+1] == 'G':
				cur_id  = f'{idn} {i+1} {seq[i]}~A'
				cur_seq = f'{seq[:i]}A{seq[i+1:]}'
		# out
			if cur_id is not None:
				print(f'>{cur_id}')
				for i in range(0, len(cur_seq), 80):
					print(cur_seq[i:i+80])
				
				
				
