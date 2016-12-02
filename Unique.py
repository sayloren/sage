"""
Script for checking if the new hg19 UCEs are present in the hg18 lifted over to hg19 UCEs

Input: Two fasta files, one for the New UCEs and one for the Old UCEs

Output: Prints out the Old and New UCEs being checked against, and whether or not they are the same

To Do:
-Read in as bedfile, use pybedtools to get the fasta sequence
-Instead of just checking that the two strings are not the same, check that one dosen't overlap at a different start point

Wren Saylor
November 28 2016

"""

import argparse
from Bio import SeqIO

def get_args():
	# Catch command-line arguments
	parser = argparse.ArgumentParser(description="Script for checking if the new hg19 UCEs are present in the hg18 lifted over to hg19 UCEs")
	parser.add_argument('-n', '--new', type=str, required=True,
                        help="An input UCE file with the new UCEs to be checked")
	parser.add_argument('-o', '--old', type=str, required=True,
                        help="An input UCE file with the hg18 UCEs to be checked against")    
	return parser.parse_args()
    
def getNewFasta(New,Old):
	NewSequences = SeqIO.parse(open(New),'fasta')
	OldSequences = SeqIO.parse(open(Old),'fasta')
	for fasta in NewSequences:
		nameNew, sequenceNew = fasta.id, fasta.seq
		strSeqNew = str(sequenceNew)
		for fasta in OldSequences:
			nameOld, sequnceOld = fasta.id, fasta.seq
			strSeqOld = str(sequnceOld)
			if strSeqNew == strSeqOld:
				print(nameNew, nameOld, "not unique")
			else:
				print (nameNew, nameOld, "unique")

	
def main():
    args = get_args()

    getNewFasta(args.new,args.old)    

if __name__ == "__main__":
     main()
