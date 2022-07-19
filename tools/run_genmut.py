import os

for f in os.listdir('data/mini_genes'):
	if f.endswith('.fa'):
		os.system(f'./genmut data/mini_genes/{f} -o mutants')
