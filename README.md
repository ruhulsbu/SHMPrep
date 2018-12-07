# SHMPrep

Use the following options to run SHMPrep:

- f: To specify the name of META file (Ex: -f meta.txt)
- n: To specify the number of Threads (Ex: -n2)
- r: To enable Referenceless Alignment (Ex: -r)
- b: To process BARCODED sequences for Referenceless Alignment (Ex: -b)
- l: To apply Local Alignment for overlapping reads (Ex: -l)
- e: To exclude the Alignment with INDELs (Ex: -e)
- q: To print the final output in FASTQ format (Ex: -q)
- p: To print BARCODES found in the output read name (Ex: -p)
- s: If enabled then Alignment will include BARCODE/PRIMERs (Ex: -s)
- d: If enabled then both CONSCOUNT and DUPCOUNT will be shown (Ex: -d)
- Q: Optional filter for BASIC Correction by per base quality (Ex: -Q20)
- M: Optional filter for Mean Quality of the sequence (Ex: -M25)
- X: To specify the INDEL Penalty for alignment (Ex: -X9)
- C: to filter out Alignments with CONSCOUNT less than given number (Ex: -C5)
- D: To filter out Alignments with DUPCOUNT less than given number (Ex: -D10)
- F: To set the length of Forward Adapter (Ex: -F24) as Prefix in Read1
- O: To set the length of Reverse Adapter (Ex: -O24) as Prefix in Read2
- K: To set the size of KBAND in Global Alignment (Ex: -K10)
- P: To specify the minimum similarity score for PRIMERs (Ex: -P70)
- L: To specify Alignment Identity to apply Smith Waterman to the Overlapped region (Ex: -L50)
- G: To specify Alignment Identity to apply Global Alignment to a Consensus (Ex: -G90)
- h: To print this help menu


