import os

models = '--dpwm /home/azhang/genomikon/isoformer/data/donor.pwm --apwm /home/azhang/genomikon/isoformer/data/acceptor.pwm --emm /home/azhang/genomikon/isoformer/data/exon.mm --imm /home/azhang/genomikon/isoformer/data/intron.mm --elen /home/azhang/genomikon/isoformer/data/exon.len --ilen /home/azhang/genomikon/isoformer/data/intron.len'


for f in os.listdir('data/mini_genes'):
	if f.endswith('.gff3'): continue
	base = os.path.splitext(f)[0]
	os.system(f'isoformer {models} data/mini_genes/{f} --introns data/mini_genes/{base}.gff3 > og_isoforms/{base}.isoform')

'''
for dirc in os.listdir('mutants'):
	os.makedirs(f'mutant_isoforms/{dirc}', exist_ok=True)
	for f in os.listdir(f'mutants/{dirc}'):
		#continue
		base = os.path.splitext(f)[0]
		os.system(f'isoformer_hm {models} mutants/{dirc}/{f} --introns data/mini_genes/{dirc}.gff3 > mutant_isoforms/{dirc}/{base}.isoform')
'''	
