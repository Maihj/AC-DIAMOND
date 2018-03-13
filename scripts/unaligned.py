#!/usr/bin/python

'''
Usage: python unaligned.py ALIGNMENT_RESULT QUERY_FILE
'''

import sys

if len(sys.argv) != 3:
    print "Usage: python unaligned.py ALIGNMENT_RESULT QUERY_FILE"
    sys.exit()

overlap = 0
dic = {}
aln_result = sys.argv[1]
all_queries = sys.argv[2]

file1 = file(aln_result, 'r')
file2 = file(all_queries, 'r')

for line in file1.readlines():
    aln = line.split()
    key = str(aln[0])
    if key not in dic:
        dic[key] = 1

id = 0
for line in file2.readlines():
    aln1 = line.split()
    key = str(aln1[0])
    if key[0] == ">":
        if key[1:] not in dic:
            print line[:-1]
            id = 1
        else:
            id = 0
    else:
        if id == 1:
            print line[:-1]
