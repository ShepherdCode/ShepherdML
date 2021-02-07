import random
import itertools
import os
import sys
import argparse

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
Pairs = ["A","C", "G", "T"]

Codons = itertools.product(Pairs, repeat=3)

all_codons = ["".join(list(codon)) for codon in Codons]

start = "ATG"
stop = ["TAA", "TAG", "TGA"]
'''
 Appends a random A, G, or C for the input frame
 e.g. appends 1 character if in frame 2, and 2 if in frame 3
'''
def shift_frame(input_seq,frame = 2):
  output = input_seq
  if frame in (1,2,3):
    for i in range(1%frame):
      output.insert(0, random.choice(("A","G","C")))
    return output
  else:
    raise ValueError("Frame Must Be 1, 2 or 3. Frame Entered: " +frame)
      

'''
  Places random codons length times, and uses choices from no_stops if coding, 
  and all_codons if not coding
'''
def codon_placer(length, coding = True):
  output = []
  no_stops = [codon for codon in all_codons if codon not in stop]

  if coding == True:
    for i in range(length):
      output.append(random.choice(no_stops))
    return output
  else:
    for i in range(length):
      output.append(random.choice(all_codons))
    return output



"""
  returns random indexes a start and stop codon should be placed
  arbitrarily chooses a variation within 1/3rd the 2nd and 4th Quintile (Is that even a word? Am I making this up?)
"""
def get_index_placement(total_codons):
  quintile = total_codons // 5
  variation = quintile // 3
  start_placement = quintile + random.randint(-variation, variation)
  stop_placement = (quintile*4) + random.randint(-variation, variation)
  return start_placement, stop_placement


"""
  Generates a random (hypothesized) coding or non-coding sequence of length characters in frame
  length: Number of characters the sequence should be
  coding: Whether or not the sequence should have stop codons placed, or randomly generated (True, False respectively)
  frame: The frame the sequence should be in which determines how many bases are appended to the sequences start. Must be 1, 2 or 3
"""  
def generate_seq(length, coding = False, frame = 1):
  codons_to_place = (length//3)

  if coding and frame in (1, 2 ,3):
    codons_to_place = codons_to_place - 2
    pre_stop_placed = codon_placer(codons_to_place)
    start_index,stop_index = get_index_placement(codons_to_place)
    pre_stop_placed.insert(start_index, start)
    pre_stop_placed.insert(stop_index, random.choice(stop))
    output = shift_frame(pre_stop_placed, frame)
    output = ''.join(output)
    return output[0:length] , coding, frame, (start_index, stop_index)

  elif not coding and frame in (1,2,3):
    placed_codons = codon_placer(codons_to_place, coding)
    output = shift_frame(placed_codons, frame)
    output = ''.join(output)
    return (output[0:length] , coding, frame)
  else:
    raise ValueError("Frame must be 1, 2 or 3")




"""
  #Todo Implement argparse

  Lines (int) is the amount of sequences a file should contain
  Coding (boolean) is the whether the sequence generated is mRNA(True) or lncRNA(False)
  Frame(int) Can be 1 2 or 3, and is the frame the sequence should be generated/ read in
  outputfile_name(string) The name of the outputfile. must end in .fasta to be compatible.
"""
def CLI_GEN(lines, coding, frame, outputfile_name):

  lines = int(lines)
  coding = bool(coding)
  frame = int(frame)

  # The minimum number of bases in a sequence
  MIN = 800

  # The Maximum number of bases in a sequence
  MAX = 1000

  if coding:
    with open(outputfile_name, 'w') as file:
      headerframe = f">GENF{frame}"
      fastaheader = [headerframe, ".1"]
      for i in range(1, lines+1):
        padded = f'{i:010}'
        sequence, target_coding, target_frame, placed_indices  = generate_seq(random.randint(MIN, MAX+1), coding, frame)
        file.write(f"{(padded).join(fastaheader)} Coding {target_coding} Frame {target_frame} Start_Index {placed_indices[0]} Stop_Index {placed_indices[1]}" +"\n")
        file.write(sequence+"\n")  
    return True
  else:
    with open(outputfile_name, 'w') as file:
      headerframe = f">GENF{frame}"
      fastaheader = [headerframe, ".1"]
      for i in range(1, lines+1):
        padded = f'{i:010}'
        sequence, target_coding, target_frame, placed_indices  = generate_seq(random.randint(MIN, MAX+1), coding, frame)
        file.write(f"{(padded).join(fastaheader)} Coding {target_coding} Frame {target_frame}" +"\n")
        file.write(sequence+"\n")  
    return False

CLI_GEN(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])