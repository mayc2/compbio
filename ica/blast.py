#-------------------------------------------------------------------------------
# Name:        blast.py
# Purpose:     simple BLAST demo
#
# Author:      Bill Thompson
#
# Created:     11/01/2015
# Copyright:   (c) Bill Thompson 2015
# Licence:     GPL
#-------------------------------------------------------------------------------

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def ReadFasta(filename):
    """
    ReadFasta - read a single fasta sequence from filename
    returns
    a BioPython SeqRecord
    """
    
    rec = list(SeqIO.parse(filename, 'fasta'))
    return rec[0]


def Main():
    e_value_cutoff = 1.0e-06
    fasta_file = 'dino1.fasta'

    seq_rec = ReadFasta(fasta_file)

    # go to NCBI and BLAST
    result_handle = NCBIWWW.qblast('blastn', 'nt', seq_rec.seq)

    # loop through alignmnents and HSPs
    blast_records = NCBIXML.parse(result_handle)
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:                
                if hsp.expect < e_value_cutoff:
                    print('****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print('e value:', hsp.expect)
                    print(hsp.query[0:75] + '...')
                    print(hsp.match[0:75] + '...')
                    print(hsp.sbjct[0:75] + '...')
    

if __name__ == '__main__':
    Main()
