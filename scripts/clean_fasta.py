import os, sys

import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

iupac = {'R' : ['A', 'G'],
         'Y' : ['C', 'T'],
         'S' : ['G', 'C'],
         'W' : ['A', 'T'],
         'K' : ['G', 'T'],
         'M' : ['A', 'C'],
         'B' : ['C', 'G', 'T'],
         'D' : ['A', 'G', 'T'],
         'H' : ['A', 'C', 'T'],
         'V' : ['A', 'C', 'G']}

def main():
    fa_path = sys.argv[1]

    for record in SeqIO.parse(fa_path, "fasta"):
        if record.id not in [str(i) for i in range(1,23)] + ['X', 'Y']:
            continue
        n = len(record)
        changed = 0
        new_seq = ""
        for i in range(0, n):
            if record[i] in ['A', 'C', 'G', 'T', 'N']:
                new_seq += record[i]
            else:
                new_seq += iupac[record[i]][0]
                changed +=1
        new_record = SeqRecord(Seq(new_seq),
                               id=record.id, name=record.name,
                               description=record.description)
        print("{}: {}/{} chars changed".format(record.id, changed, n), file=sys.stderr)
        SeqIO.write(new_record, sys.stdout, "fasta")

if __name__ == '__main__':
    main()
