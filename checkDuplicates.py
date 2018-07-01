#about: goes through file "filename" and checks if any two lines are same, if yes, prints them
#usage: python checkDuplicates.py filename

import sys
args = sys.argv[1:]
infilename = args[0]
print("Running check Duplicates on file: "+infilename)
f = open(infilename,"r")
l = f.read().split('\n')
tmp = {}
for i in range(0,len(l)):
	if l[i] in tmp:
		print(l[i])
	else:
		tmp[l[i]] = 0
f.close()
