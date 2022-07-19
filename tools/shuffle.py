import random
import sys
from readfasta import read_record

pep    = ''
pep_id = ''
for name, seq in read_record(sys.argv[1]):
	pep += seq
	pep_id = name

pep = list(pep)
random.shuffle(pep)
pep = ''.join(pep)

print(f'>{pep_id}')
for i in range (0, len(pep), 60):
	print(pep[i:i+60])
