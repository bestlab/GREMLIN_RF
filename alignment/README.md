Uses hhblit to compute the alignments. The files will be read in from [fasta](fasta).

A typical use of hhblits would look like:

    hhblits -i input.seq -d hhdir -oa3m alignment.a3m.aln -n 4 -diff inf -cov 75 -e 0.0000000001 

Fix the symbolic links so hhblits points to a valid executable and hhdir is a valid hhblits directory of alignments.

./hhblits -i fasta/1a3a.fasta -d hhdir/uniprot20_2013_03 -oa3m alignment.a3m.aln -n 4 -diff inf -cov 75 -e 0.0000000001
Search results will be written to fasta/1a3a.hhr
