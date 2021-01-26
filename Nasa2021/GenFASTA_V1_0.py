import random
import itertools
import os
import sys

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
Pairs = ["A","C", "G", "T"]

Codons = itertools.product(Pairs, repeat=3)

all_codons = ["".join(list(codon)) for codon in Codons]

start = "ATG"
stop = ["TAA", "TAG", "TGA"]

# Code Smell Alert!
"""
  a function to place codons in each frame in such a way that avoids stop codons
  by detecting if the previous codon could possibly create a stop codon in a given frame,
  and choosing an option that will not lead to a stop codon. For Frame 1 this is straight-forward,
  as you can use all codons but stop codons. In frame 2, codons ending in TG or TA can not precede A or G.
  In frame 3, codons ending with T can not precede AA, AG or GA.

  returns an array of codons in the specified frame without a stop codon
"""
def codon_placer(length, frame=1):
  output = []
  stop_ending_frame_2 = ["TG", "TA"]
  stop_ending_frame_3 = "T"

  safe_choice_frame_1 = [codon for codon in all_codons if codon not in stop]
  safe_choice_frame_2 = [codon for codon in all_codons if codon[0] not in ("A", "G")]
  safe_choice_frame_3 = [codon for codon in all_codons if codon[0:3] not in ("AA", "AG", "GA")]

  if frame==1:
    for i in range(length):
      output.append(random.choice(safe_choice_frame_1))
    return output

  if frame==2:
    for i in range(length):
      if i == 0:
        output.append(random.choice(all_codons))
      else:
        prev_codon = output[i-1]
        possible_stop = prev_codon in stop_ending_frame_2
        if possible_stop:
          output.append(random.choice(safe_choice_frame_2))
        else:
          output.append(random.choice(all_codons))
    return output

  if frame==3:
    for i in range(length):
      if i == 0:
        output.append(random.choice(all_codons))
      else:
        prev_codon = output[i-1]
        possible_stop = prev_codon[2] == stop_ending_frame_3
        if possible_stop:
          output.append(random.choice(safe_choice_frame_3))
        else:
          output.append(random.choice(all_codons))
    return output






"""
  returns random indexes a start and stop codon should be placed
  arbitrarily chooses a variation within 1/3rd the 2nd and 4th Quintile (Is that even a word? Am I making this up?)
"""
def get_index_placement(total_codons):
  quintile = total_codons//5
  variation = quintile // 3
  start_placement = quintile + random.randint(-variation, variation)
  stop_placement = (quintile*4) + random.randint(-variation, variation)
  return start_placement, stop_placement

def get_pair_start(frame = 2):
  first_choice_frame_2 = [codon for codon in all_codons if codon[1:3] == "AT"]
  second_choice_frame_2 = [codon for codon in all_codons if codon[0] == "G"]

  first_choice_frame_3 = [codon for codon in all_codons if codon[2] == "A"]
  second_choice_frame_3 = [codon for codon in all_codons if codon[0:2] =="TG"]

  if frame == 2:
    return random.choice(first_choice_frame_2), random.choice(second_choice_frame_2)
  elif frame == 3:
    return random.choice(first_choice_frame_3), random.choice(second_choice_frame_3)
  else:
    return "ATG"
  
# These should be constants defined outside of these functions and used as "lookup tables" so these don't get repeatedly constructed
# Every call. That's the idea behind using list builders of choices so much anyway.
# But for now, while I get it working, I appreciate them being right inside the function for the sake of cohesiveness
def get_pair_stop(frame = 2):
  first_choice_frame_2 = [codon for codon in all_codons if codon[1:3] in ("TG, TA")]
  second_choice_frame_2_TG = [codon for codon in all_codons if codon[0] == "A"]
  second_choice_frame_2_TA = [codon for codon in all_codons if codon[0] in ("A", "G")]

  first_choice_frame_3 = [codon for codon in all_codons if codon[2] == "T"]
  second_choice_frame_3 = [codon for codon in all_codons if codon[0:2] in ("AA", "AG", "GA")]

  if frame == 2:
    first = random.choice(first_choice_frame_2)
    if first[1:3] == "TG":
      second = random.choice(second_choice_frame_2_TG)
      return first, second
    else:
      second = random.choice(second_choice_frame_2_TA)
      return first, second

  if frame == 3:
    return random.choice(first_choice_frame_3), random.choice(second_choice_frame_3)
  
def generate_seq(length, coding = False, frame = 1):
  codons_to_place = (length//3)
  if coding and frame == 1:
    codons_to_place = codons_to_place - 2
    pre_stop_placed = codon_placer(codons_to_place, frame)
    start_index,stop_index = get_index_placement(codons_to_place)
    pre_stop_placed.insert(start_index, start)
    pre_stop_placed.insert(stop_index, random.choice(stop))
    output = ''.join(pre_stop_placed)
    return output[0:length]
  elif coding and frame == 2:
    codons_to_place = codons_to_place - 4
    pre_stop_placed = codon_placer(codons_to_place, frame)
    start_index,stop_index = get_index_placement(codons_to_place)
    _AT, G__ = get_pair_start(2)
    pre_stop_placed.insert(start_index, _AT)
    pre_stop_placed.insert(start_index+1, G__)
    _TA_TG, A__G__ = get_pair_stop(2)
    pre_stop_placed.insert(stop_index, _TA_TG)
    pre_stop_placed.insert(stop_index+1, A__G__)
    output = ''.join(pre_stop_placed)
    return output[0:length] 
  elif coding and frame == 3:
    codons_to_place = codons_to_place - 4
    pre_stop_placed = codon_placer(codons_to_place, frame)
    start_index,stop_index = get_index_placement(codons_to_place)
    __A, TG_ = get_pair_start(3)
    pre_stop_placed.insert(start_index, __A)
    pre_stop_placed.insert(start_index+1, TG_)
    __T, AA_AG_GA_ = get_pair_stop(3)
    pre_stop_placed.insert(stop_index, __T)
    pre_stop_placed.insert(stop_index+1, AA_AG_GA_)
    output = ''.join(pre_stop_placed)
    return output[0:length]
  else:
    output = []
    for i in range(codons_to_place):
      output.append(random.choice(all_codons))
    return "".join(output)[0:length]





def CLI_GEN(lines, coding, frame, outputfile_name):

  lines = int(lines)
  coding = bool(coding)
  frame = int(frame)

  with open(outputfile_name, 'w') as file:
    headerframe = f"GENF{frame}"
    fastaheader = [headerframe, ".1"]
    for i in range(1, lines+1):
      padded = f'{i:010}'
      file.write((padded).join(fastaheader)+"\n")
      file.write(generate_seq(1000, coding, frame)+"\n")  
  return True

CLI_GEN(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])