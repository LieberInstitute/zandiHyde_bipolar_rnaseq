
import sys

infile = sys.argv[1]
outfile = "./"+infile.split('/')[-1].replace(".count",".junc")

with open(infile, 'r') as inf, open(outfile, 'w') as out:
	for line in inf.readlines():
		s = line.rstrip().split('\t')
		newline = s[0]+'\t'+s[1]+'\t'+s[2]+'\t'+"."+'\t'+s[4]+'\t'+s[3]+"\n"
		out.write(newline)

