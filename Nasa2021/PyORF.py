import itertools
import os
import sys
#import re

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#pattern = re.compile(r'(ATG)((?:[ATGC]{3})*?)(TAG|TGA|TAA)')

#test_sequence = "ACATAGCTGGAATATCGTAAGAGATTACACGGTAATGACGCACCGGCCATCTTTCGGTCCATTGAAACCATGTGCTGTGTTGAGCATTCGACCGTTTAAATAAGGTACACGGTACTTCACTTGAATTCTCTCTCAAGCTTGGATCGGCCTTGCGCGTACATTGGTTATAAGAATCTGCGCTCGTTGTTTCGGTCAGATTTAACGGGCTGACGGGACCGACCTGAAAACCTCCCTTCAACATTTCGTAAAGGCCTAATGGGAGCGTTGCCTTCTCCTTGTAGAGCTTTACGTCTGTGTGTACCTGCCCCCAAGGTTACCAAAGACAGCCACATCCAGGGTATCAGTTCTCTGCTCAAAGGTAAAGGTGGTGATTTACGCGTGAACGAGGTGTTCGATCGCCCTGCTCATCGGTATCGAGTTCGAAGGAATGATCCAGTGACTGTCAAAATTGCTTCTTGTTCGGCGCCCATTAGCTACCTGTTCAGATCGATCCTCACCTGTTTCCGCATGCTACGGAGCGGCTGGGCGGGTATCAATCTCCGTGTCACTAGCGAATCTCGGAGAACCGAATTCCGATCAACCAAACTGAGTTGGTTTCAACTCGCAGGCACTGTGGTACAGATTATTGCCCATGATTTCCAACCTCGATTTGAACCAGAATCCGGAGGGATGGGGAGTGAAGTCCTGAATACTTTGATTATCCGCTTACCTCGTAGTGGGACTTCCCACCTCGGCCGAGTGACAGCGTACGTATGGATAATGTACCTTACGATGAATCAAGGTTAGGACAAGCCGAAGGGATATACGGCGTTCATAGTCGCGCGATAACTCCGGCTGTAGTACTAGTTTTTAAGGCCATACGCACACAGTGGCGCTCGTCGGTCTGGCGTCTGCAGTTTACAGTTCGGATAGTGAATCAAGTTGTCCTGCCTCTGCAAGGTATGGAGCGTTCGCGATTCCACCGTGTGCCCCACCTTGGCCACAAACGGAAACGGGGTAG"

def chunkstring(string, length, frame=1):
    frame = frame-1
    string = string[frame::]
    return list((string[0+i:length+i] for i in range(0, len(string), length)))

def find_codons(sequence, start_codon = "ATG", stop_codons = ("TAG", "TAA", "TGA")):
  found_start = []
  found_stop = []
  in_sequence = False

  for i in range(len(sequence)):
    if sequence[i] == start_codon and not in_sequence:
      found_start.append(i)
      in_sequence = True
    if sequence[i] in stop_codons and in_sequence:
      found_stop.append(i)
      in_sequence = False
  found_codons = list(itertools.zip_longest(found_start, found_stop))
  return found_codons


def orf_length(pair):
  if pair[1] is None:
    return 0
  else:
    return pair[1] - pair[0]

def biggest_orf(orf_list, frame):
  pair = max(orf_list, key=orf_length)
  correct_pair = calculate_indexes(pair, frame)
  pair = correct_pair
  return (pair[1] - pair[0]+1), pair, frame

def calculate_indexes(pair, frame):
  CODON = 3
  EMPTY = (0,1)
  if None in pair:
    return EMPTY
  if frame==1:
    return (pair[0]*3+frame, pair[1]*3+CODON)
  if frame==2:
    return (pair[0]*3+frame, pair[1]*3+1+CODON)
  if frame==3:
    return (pair[0]*3+frame, pair[1]*3+2+CODON)

def max_orf_by_frame(sequence, frame):
  codon_list = chunkstring(sequence,3,frame)
  orfs = find_codons(codon_list)
  max_orf = biggest_orf(orfs, frame)
  return max_orf

def max_orfs_by_frame(sequence):
  frames_with_orfs = []
  all_maxes = []
  for i in range(1,4):
    found_codons = find_codons(chunkstring(sequence, 3, i))
    if found_codons:
      frames_with_orfs.append(i)
  for frame in frames_with_orfs:
    all_maxes.append(max_orf_by_frame(sequence, frame))
  return all_maxes

def get_biggest_orf_all_frames(sequence):
  all_maxes = max_orfs_by_frame(sequence)
  biggest_max = max(all_maxes, key= lambda x: x[0])
  return biggest_max


def validate_orfs_in_frame(filename, frame):
  with open(filename, 'r') as file:
    all_lines = file.read().splitlines()
    sequences = all_lines[1::2]
    valid = True
    for i in range(len(sequences)):
      #line = (i+1) *2
      #print(f"Last line read: {line} in {filename}")
      orf_attrib = get_biggest_orf_all_frames(sequences[i])
      if orf_attrib[2] != frame:
        valid = False
    return valid

def print_all_orfs(sequence):
  frame_1 = chunkstring(sequence, 3, 1)
  frame_2 = chunkstring(sequence, 3, 2)
  frame_3 = chunkstring(sequence, 3, 3)

  print("All ORFs found in frame 1: ", find_codons(frame_1))
  print("All ORFs found in frame 2: ", find_codons(frame_2))
  print("All ORFs found in frame 3: ", find_codons(frame_3))

print("ORFs generated in Frame 1 are all Valid: ",validate_orfs_in_frame("frame1.fasta", 1) )
print("ORFs generated in Frame 2 are all Valid: ",validate_orfs_in_frame("frame2.fasta", 2) )
print("ORFs generated in Frame 3 are all Valid: ",validate_orfs_in_frame("frame3.fasta", 3) )
