import argparse

parser = argparse.ArgumentParser(description='Make Snakefiles')
parser.add_argument('pdb', type=str, metavar='<path>', help='input pdb file')
parser.add_argument('-w', type=int, metavar='<int>', required=False,
	help='window size for sliding window average (default: 5)', default=5)
arg = parser.parse_args()

def get_avg(scores):
	return sum(scores)/len(scores)
	
def get_sw_avg(scores):
	for i in range(0, len(scores)-arg.w):
		
	
scores = []
with open(arg.pdb) as fh:
	for line in fh.readlines():
		line = line.split()
		if len(line) < 11: continue
		score = float(line[10])
		scores.append(score)

print(get_avg(scores))
