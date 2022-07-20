import argparse
from readfasta import read_record

parser = argparse.ArgumentParser(description='select mutant isoforms in frame of annotated original isoform')
parser.add_argument('mutfa', type=str, metavar='<fasta>', help='input path to mutant fasta file')
parser.add_argument('mutiso', type=str, metavar='<isoform>', help='input path to mutant.isoform')
parser.add_argument('fa', type=str, metavar='<fasta>', help='input path to og fasta file')
parser.add_argument('gff', type=str, metavar='<gff>', help='input path to gff')
arg = parser.parse_args()

def get_cdss(exons):
	rnaseq = ''
	orf    = ''
	cdss   = []
	for exon in exons:
		rnaseq += mutseq[exon[0]:exon[1]+1]
	print(f'>rnaseq\n{rnaseq}')
	
	for i in range(0, len(rnaseq)-2):
		start = (rnaseq[i:i+3] == 'ATG')
		if start:
			cur_orf = ''
			for j in range(i, len(rnaseq)-2, 3):
				cur_orf += rnaseq[j:j+3]
				if rnaseq[j:j+3] in stop_codons:
					if len(cur_orf) > len(orf):
						orf = cur_orf
						orf_pos = [i,j+3]
					break
	print(f'>orf\n{orf}')
	for exon in exons: print(exon, end = ' ')
	print(f'\n{orf_pos}')
	print(f'orf length: {len(orf)}')
	
	r = len(orf)
	for idx, exon in enumerate(exons):
		if idx == 0:
			if exon[0] + orf_pos[0] + r > exon[1]:
				cdss.append([exon[0]+orf_pos[0],exon[1]])
				r = r - (exon[1] - (exon[0] + orf_pos[0])) -1
			else:
				cdss.append([exon[0]+orf_pos[0],exon[0]+orf_pos[0]+r]-1)
				break
		else:
			if exon[0] + r > exon[1]:
				cdss.append(exon)
				r = r - (exon[1] - exon[0]) -1
			else:
				cdss.append([exon[0],exon[0]+r-1])
				break
	
	for cds in cdss: print(cds, end = ' ')
	print('')
	
	test = ''
	for cds in cdss: test += mutseq[cds[0]:cds[1]+1]
	print(test)

stop_codons = ['TAA', 'TAG', 'TGA']
mutseq = ''
orgseq = ''
for idn, seq in read_record(arg.mutfa): mutseq += seq
for idn, seq in read_record(arg.fa):    orgseq += seq

gff_cdss = []
with open(arg.gff) as fh:
	while True:
		line = fh.readline()
		if line == '': break
		line = line.rstrip()
		if 'WormBase' not in line: continue
		fields = line.split()
		if fields[2] == 'CDS':
			gff_cdss.append([int(fields[3])-1,int(fields[4])-1])

print(f'>mutseq\n{mutseq}')

with open(arg.mutiso) as fh:
	exons = []
	while True:
		line = fh.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('#'): continue
		fields = line.split()
		if len(fields) < 9: continue
		if fields[2] == 'mRNA':
			if len(exons) > 0:
				orf = get_cdss(exons)
				break ########################
			exons = []
		if fields[2] == 'exon':
			exons.append([int(fields[3])-1, int(fields[4])-1])
	#orf = get_cdss(exons)

