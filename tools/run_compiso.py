import os

for og in os.listdir('og_isoforms'):
	og_base = os.path.basename(og)
	og_name = os.path.splitext(og_base)[0]
	for mut in os.listdir(f'mutant_isoforms/{og_name}'):
		os.system(f'cmpiso_hm og_isoforms/{og} mutant_isoforms/{og_name}/{mut} >> cmpiso_out/{og_name}.cmp')
