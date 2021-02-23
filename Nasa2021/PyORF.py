import re

pattern = re.compile(r'(ATG)((?:[ATGC]{3})*?)(TAG|TGA|TAA)')
test_sequence = ""


def orfs(dna):
    return list(set(pattern.findall(dna)))


def orf_indexes(sequence):
  index_matches = []
  for match in pattern.finditer(sequence):
    start = match.start()
    end = match.end()
    length = match.end() - match.start()
    index_matches.append((start, end, length))
  return index_matches

def get_biggest(indexes):
    Biggest_ORF = sorted(indexes, key = lambda orf: orf[2], reverse = True)[0]
    return Biggest_ORF

def contains_orf(sequence, min, max):
  indexes = orf_indexes(sequence)
  start, stop, length = get_biggest(indexes)
  if length >= min and length <= max:
    return True
  else:
    return False
    
def get_sequence_attributes(sequence):
  seq_len = len(sequence)
  ORF_START, ORF_STOP, ORF_LEN = get_biggest(orf_indexes(sequence))
  UTR_5 = 0
  UTR_3 = 0
  if ORF_START != 0:
      UTR_5 = ORF_START
  if ORF_STOP + 3 != seq_len:
      UTR_3 = seq_len - (ORF_STOP + 3)
  return seq_len, UTR_5, ORF_LEN, UTR_3

is_valid = contains_orf(test_sequence, 200, 2000)
print(str(is_valid))
seq_len, UTR_5, ORF_LEN, UTR_3 = get_sequence_attributes(test_sequence)

print(f"SeqLen {seq_len} 5UTR {UTR_5} ORF {ORF_LEN} 3UTR {UTR_3}")

