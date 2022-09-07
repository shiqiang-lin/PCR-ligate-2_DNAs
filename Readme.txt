Notes to the script

1. Requirements
1) macOS Monterey 12.3.1 or higher (Apple Inc., CA, USA);
2) Python 3.10.4;
3) Biopython 1.79;
4) Overlap_mutagenic_primer_v1.py;
5) ECdnaQ.fasta and TBdnaN.fasta;
6) IDLE (Integrated Development and Learning Environment) (www.python.org) for modifying the script when necessary.

2. Running method
1) Make a new directory, copy the Overlap_mutagenic_primer_v1.py, ECdnaQ.fasta and TBdnaN.fasta to the directory;
2) open a terminal and cd to the new directory in 1);
3) input the following command and press enter, wait for a few seconds, a new file storing the primers will appear in the directory.

python3.10 Overlap_ligation_primer_v1.py ECdnaQ.fasta TBdnaN.fasta

3. Adjusting the parameters of the script
Here are the parameters of the script. The number in the left indicates the position of the parameter in the script.

line 90: Tm_favorite = 63
line 119: cnd1 = abs(Tm_primer5 - Tm_overlap3)
line 120: cnd2 = abs(Tm_overlap5 - Tm_primer3)
line 121: cnd3 = abs(Tm_primer5 - Tm_primer3)
line 122: cnd4 = abs((Tm_primer5 + Tm_primer3)/2 - Tm_overlap - 3)
line 123: cnd5 = abs(Tm_overlap - Tm_favorite)
line 124: cnd6 = abs((Tm_primer5 + Tm_primer3)/2 - Tm_favorite)


The script relies on the sum of cnd1 to cnd6 (line 119-124) to sort the usefulness of the primers. The ideal primers should have a low value for the sum of cnd1, cnd2, cnd3, cnd4, cnd5 and cnd6.

Commonly, the annealing temperature for PCR is 3 degrees less than the theoretical annealing temperature. In our script, the Tm_favorite is set 63 degrees Celsius (line 90 of the script), therefore, the ideal annealing temperature for PCR running is 60 degrees Celsius. This temperature is not too high or too low, which is good for most PCR running. Therefore, the idea of cnd5 that the Tm difference between Tm_overlap and Tm_favorite is as small as possible. The idea of cnd6 is that the Tm difference between the (Tm_primer5 and Tm_primer3)/2 and Tm_favorite is as small as possible. However, the value of Tm_favorite can be adjusted if users feel it necessary.

The idea of cnd4 (line 122) is that the difference between (Tm_primer5 + Tm_primer3)/2 and Tm_overlap should be close to 3 degrees Celsius. For a successful overlap PCR, the bridging temperature (Tm_overlap) should be a little lower than the Tm_primer5 and Tm_primer3, which ensures a smooth bridging. Here we set the difference 3 degrees Celsius and it can be set 2 or 4 or other number degrees Celsius if users feel it necessary.

4. Trouble shooting
To run the script successfully, users need to ensure the following:
1) The fasta files contain DNA sequences only in ‘A’, ‘T’, ‘C’ and ‘G’ and not any other letters;
2) In the command, input the FIRST fasta file, input a blank, and then input the SECOND fasta file. The order of the input fasta files is important.





