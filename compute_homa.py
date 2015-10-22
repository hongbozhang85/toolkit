#!/usr/bin/python

import sys

ALPHA = 257.7
ROPT = 1.388

def distance(atom_1, atom_2, coords):
    atom_1 = int(atom_1)
    atom_2 = int(atom_2)
    atom_1_coords = [float(coords[atom_1 - 1].split()[1]), \
                     float(coords[atom_1 - 1].split()[2]), \
                     float(coords[atom_1 - 1].split()[3])]
    atom_2_coords = [float(coords[atom_2 - 1].split()[1]), \
                     float(coords[atom_2 - 1].split()[2]), \
                     float(coords[atom_2 - 1].split()[3])]
    return ((atom_1_coords[0] - atom_2_coords[0])**2 + \
            (atom_1_coords[1] - atom_2_coords[1])**2 + \
            (atom_1_coords[2] - atom_2_coords[2])**2)**(.5)

if len(sys.argv) < 3:
    print "Correct Input Format: coords.xyz .ring"
    exit()

coord_filename = sys.argv[1]
ring_filename = sys.argv[2]

atom_coords = []
ring_index = []
homa_values = []

coord_file = open(coord_filename)
for line in coord_file:
    if len(line.split()) < 4: continue
    else:
        atom_coords += [line]
coord_file.close()

ring_file = open(ring_filename)
for line in ring_file:
    ring_index += [line]
ring_file.close()

for ring in ring_index:
    ring_distances = []
    n = len(ring.split())
    sum = 0
    i = 0
    while i < n:
        if i == (n - 1):
            ring_distances += [distance(ring.split()[i],ring.split()[0], atom_coords)]
        else:
            ring_distances += [distance(ring.split()[i],ring.split()[i+1], atom_coords)]
        i = i + 1

    for dist in ring_distances:
        sum += (ROPT - dist)**2
    homa_values += [1.0 - (ALPHA/n)*sum]

homa_aver=0

for element in homa_values:
    homa_aver += element 

Nring=len(ring_index)
homa_aver_ring = homa_aver/Nring


output_file = open(ring_filename.rstrip('.ring') + ".out", "w")
for element in homa_values:
    output_file.write(str(element) + '\n')

output_file.write('\n')
output_file.write('HOMA cumulative    ' + str(homa_aver) + '\n')
output_file.write('HOMA average    ' + str(homa_aver_ring) + '\n')

output_file.close()
