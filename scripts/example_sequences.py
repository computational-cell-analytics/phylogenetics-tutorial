from pathlib import Path
from Bio import SeqIO

path = Path("../data/primate_mt_cytb.fasta")
records = list(SeqIO.parse(path, "fasta"))

start = 10
length = 10

human = str(records[0].seq)
chimp = str(records[1].seq)

print("Human:")
print(human[start:(start + length)])

print("Chimp:")
print(chimp[start:(start + length)])
