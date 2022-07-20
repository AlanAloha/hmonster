import argparse
from readfasta import read_record

parser = argparse.ArgumentParser(description='select mutant isoforms in frame of annotated original isoform')
parser.add_argument('mutfa', type=str, metavar='<fasta>', help='input path to mutant fasta file')
parser.add_argument('mutiso', type=str, metavar='<isoform>', help='input path to mutant.isoform')
parser.add_argument('fa', type=str, metavar='<fasta>', help='input path to og fasta file')
parser.add_argument('gff', type=str, metavar='<gff>', help='input path to gff')
parser.add_argument('-t', '--threshold', type = int, metavar='<int>', required = False,
	default=0.01, help='threshold for selecting isoform candidates (default: 0.01)')
arg = parser.parse_args()

def get_cdss(exons):
	rnaseq = ''
	orf    = ''
	cdss   = []
	for exon in exons:
		rnaseq += mutseq[exon[0]:exon[1]+1]
	#print(f'>rnaseq\n{rnaseq}')
	
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
	'''
	print(f'>orf\n{orf}')
	for exon in exons: print(exon, end = ' ')
	print(f'\n{orf_pos}')
	print(f'orf length: {len(orf)}')
	'''
	tot = 0
	exon_start = None
	orf_start = orf_pos[0]
	for idx, exon in enumerate(exons):
		exlen = exon[1]-exon[0]
		tot += exlen
		if tot > orf_pos[0]:
			exon_start = idx
			break
		orf_start -= exlen + 1
	
	r = len(orf)
	for idx, exon in enumerate(exons[exon_start:]):
		if idx == 0:
			if exon[0] + orf_start + r > exon[1]:
				cdss.append([exon[0]+orf_start,exon[1]])
				r = r - (exon[1] - (exon[0] + orf_start)) -1
			else:
				cdss.append([exon[0]+orf_start,exon[0]+orf_start+r-1])
				break
		else:
			if exon[0] + r > exon[1]:
				cdss.append(exon)
				r = r - (exon[1] - exon[0]) -1
			else:
				cdss.append([exon[0],exon[0]+r-1])
				break
	
	'''
	for cds in cdss: print(cds, end = ' ')
	print('')
	
	test = ''
	for cds in cdss: test += mutseq[cds[0]:cds[1]+1]
	print(test)
	'''
	return cdss

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

#print(f'>mutseq\n{mutseq}')

mut_pos = None
mut_typ = None
isoforms = []
with open(arg.mutiso) as fh:
	isoform = []
	while True:
		line = fh.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('#'):
			if line.startswith('# name'):
				mut_pos = int(line.split()[5])
				mut_typ = line.split()[6]
			continue
		fields = line.split()
		if len(fields) > 8 and fields[2] == 'mRNA' and float(fields[5]) > arg.threshold:
			isoform.append(line)
		if len(fields) < 1 and len(isoform) > 1:
			isoforms.append(isoform)
			isoform = []
		if len(isoform) > 0 and fields[2] != 'mRNA': isoform.append(line)


for cds in gff_cdss: print(cds, end=' ')
print('')
for cds in gff_cdss: print(f'{orgseq[cds[0]:cds[1]+1]}', end = '')
print('')	


for isoform in isoforms:
	exons = []
	introns = []
	for line in isoform:
		fields = line.split()
		if fields[2] == 'exon':
			exons.append([int(fields[3])-1, int(fields[4])-1])
		if fields[2] == 'intron':
			introns.append(int(fields[3])-1)
			introns.append(int(fields[4])-1)
	if mut_pos in introns:
		mut_cdss = get_cdss(exons)
		for cds in mut_cdss: print(cds, end = ' ')
		print('')
		for cds in mut_cdss: print(f'{mutseq[cds[0]:cds[1]+1]}', end = '')
		print('')
	
