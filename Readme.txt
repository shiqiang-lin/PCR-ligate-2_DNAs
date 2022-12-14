Notes to the script

1. Requirements
1) macOS Monterey 12.3.1 or higher (Apple Inc., CA, USA);
2) Python 3.10.4;
3) Biopython 1.79;
4) Overlap_mutagenic_primer_v1.py;
5) ECdnaQ.fasta and TBdnaN.fasta;
6) IDLE (Integrated Development and Learning Environment) (www.python.org) for modifying the script when necessary.

2. Running method
1) Make sure that the Biopython 1.79 has already been installed to the python 3.10.4. If not, 
then open a terminal and input the following command and press enter.

pip3.10 install biopython

For more details about the Biopython 1.79, the readers can visit the link 
https://biopython.org/wiki/Download.

2) Make a new directory, copy the Overlap_mutagenic_primer_v1.py, ECdnaQ.fasta and TBdnaN.fasta to the directory;
3) open a terminal and cd to the new directory in 1);
4) input the following command and press enter, wait for a few seconds, a new file storing the primers will appear in the directory.

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


The Tm value of the overlapping region between the two primers should be moderate. We set the favorite Tm value 63°C (line 90 of the script). Commonly, the annealing temperature for PCR is 3 degrees less than the theoretical annealing temperature. Therefore, the ideal annealing temperature for PCR running is 60°C. This temperature is not too high nor too low, which is good for most PCR running. The idea of cnd5 is that the Tm difference between Tm_overlap and Tm_favorite is as small as possible. The idea of cnd6 is that the Tm difference between the (Tm_primer5 and Tm_primer3)/2 and Tm_favorite is as small as possible. However, the value of Tm_favorite can be adjusted if users feel necessary.

The idea of cnd4 (line 122) is that the difference between (Tm_primer5 + Tm_primer3)/2 and Tm_overlap should be close to 3°C. For a successful overlap PCR, the bridging temperature (Tm_overlap) should be a little lower than the Tm_primer5 and Tm_primer3, which ensures a smooth bridging. Here we set the difference 3°C and it can be set 2 or 4 or other number °C if users feel necessary.

4. Trouble shooting
To run the script successfully, users need to ensure the following:
1) The fasta files contain DNA sequences only in ‘A’, ‘T’, ‘C’ and ‘G’ and not any other letters;
2) In the command, input the FIRST fasta file, input a blank, and then input the SECOND fasta file. The order of the input fasta files is important.


5. Additional note
We also have another script that is useful for gene cloning in the lab. The GitHub Repo link is 
https://github.com/shiqiang-lin/sdm-mutagenesis and the reference is 

Xiaofang Huang, Liangting Xu, Chuyun Bi, Lili Zhao, Limei Zhang, Xuanyang Chen, Shiqian Qi, Shiqiang Lin. Designing overlap extension PCR primers for protein mutagenesis: a programmatic approach. In: Currin, A., Swainston, N. (eds) Directed Evolution. Methods in Molecular Biology, vol 2461. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-2152-3_1.

Thank you!




 
