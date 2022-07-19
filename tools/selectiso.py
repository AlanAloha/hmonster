import json
import collections

mutants = {}
true_mutants = {}
with open('cmpiso_all') as fh:
	for line in fh.readlines():
		dist = float(line.split()[0])
		mut  = line.split()[1]
		gene = mut.split('_')[2]
		if dist > 0.5:
			if gene not in mutants: mutants[gene] = {}
			mutants[gene][mut] = str(dist)
			true_mut = False
			with open(f'mutant_isoforms/{gene}/{mut}.isoform') as mf:
				idn = mf.readline().split()
				pos = int(idn[5])
				typ = idn[6]
				for line in mf.readlines():
					if line.startswith('#'): continue
					if len(line.split()) < 8: continue
					if line.split()[2] != 'intron': continue
					if typ == 'd' and int(line.split()[3]) == pos+1: true_mut = True
					if typ == 'a' and int(line.split()[4]) == pos+1: true_mut = True
			if true_mut:
				if gene not in true_mutants: true_mutants[gene] = {}
				true_mutants[gene][mut] = str(dist)
						
#mutants = dict(sorted(mutants.items(), key=lambda item: item[1], reverse=True))
#true_mutants = dict(sorted(true_mutants.items(), key=lambda item: item[1], reverse=True))
#print(f'mutants: {len(mutants)} genes, {sum(mutants.values())} mutants')
#print(f'true mutants: {len(true_mutants)} genes, {sum(true_mutants.values())} mutants')
print(json.dumps(true_mutants, indent=4))
