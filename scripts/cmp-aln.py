#!/usr/bin/python

'''
Usage: python cmp-aln.py BLASTX_ALIGNMENT_RESULT OTHER_TOOL_ALIGNMENT_RESULT
'''

import sys

if len(sys.argv) != 3:
    print "Usage: python cmp-aln.py BLASTX_ALIGNMENT_RESULT OTHER_TOOL_ALIGNMENT_RESUL"
    sys.exit()

overlap = 0
dic = {}
blastx_file_name = sys.argv[1]
otherTool_file_name = sys.argv[2]

file1 = file(blastx_file_name, 'r')
file2 = file(otherTool_file_name, 'r')

for line in file2.readlines():
    aln = line.split()
    # filter evalue > 0.001
    if float(line.split()[10]) > float(0.001): continue

    key = str(aln[0])+"&"+str(aln[1])
    if key in dic:
        if int(aln[6]) > int(aln[7]):
            dic[key].append([int(aln[7]), int(aln[6]), int(aln[8]), int(aln[9])])
        else:
            dic[key].append([int(aln[6]), int(aln[7]), int(aln[8]), int(aln[9])])
    else:
        if int(aln[6]) > int(aln[7]):
            dic[key] = [[int(aln[7]), int(aln[6]), int(aln[8]), int(aln[9])]]
        else:
            dic[key] = [[int(aln[6]), int(aln[7]), int(aln[8]), int(aln[9])]]
            
count = 0
aln1 = []
query = {}
for line in file1.readlines():
    if float(line.split()[10]) > float(0.001): continue
    
    count += 1
    if int(line.split()[6]) < int(line.split()[7]):
        aln1.append([line.split()[0], line.split()[6], line.split()[7], line.split()[1], line.split()[8], line.split()[9], 0])
    else:
        aln1.append([line.split()[0], line.split()[7], line.split()[6], line.split()[1], line.split()[8], line.split()[9], 0])
    
# compared with BLASTX, the number of overlapped alignments in other tool
for i in range(0, len(aln1)):
    id = 0
    
    key = str(aln1[i][0])+"&"+str(aln1[i][3])
    key2 = str(aln1[i][0])
    if key in dic:
        for aln2 in dic.get(key):
            if int(aln1[i][1]) > aln2[1]: continue
            if int(aln1[i][2]) < aln2[0]: continue
            if int(aln1[i][4]) > aln2[3]: continue
            if int(aln1[i][5]) < aln2[2]: continue
            overlap += 1
            aln1[i][6] = 1
            id = 1
            break
        
        if id == 1:
            if key2 not in query:
                query[key2] = 1
                
print "Sensitivity: " + str(float(overlap) / float(count) * 100) + "%"

file1.close()
file2.close()
