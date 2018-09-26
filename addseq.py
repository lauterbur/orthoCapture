#!/bin/python

import sys

with open(sys.argv[1], "r") as f:
	for i in range(244): # number of lines in file divided by two
		header = f.readline()
		sequence = f.readline()
		print(header[1:-4])
		if header[0] is not ">" or sequence[0] is ">":
			print("Lines out of order")
			sys.exit()
		else:
#			if header[1:-4] in "CMC2_3_":
			new = header[1:-4]+"Mustelid.fa"
			print(new)
			with open(new, "a") as n:
				n.write(header[:-4]+"Ailurus\n")
				n.write(sequence+"\n")
#			else:
#				print(header+" not in file")
