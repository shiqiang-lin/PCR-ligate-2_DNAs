#
#usage:
#python3.10 Overlap_ligation_primer_v1.py first_DNA.fasta second_DNA.fasta
#


import sys
import os,shutil
import string
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC,MeltingTemp

#get current path
current_path = os.getcwd()
print("Current path is %s." % current_path)
print("\n")


#get the fasta files of the two fragments
first_DNA_file = sys.argv[1]
second_DNA_file = sys.argv[2]



#get 1st DNA sequence and print
first_DNA_sequence_str = str(SeqIO.read(first_DNA_file, "fasta").seq)

if len(first_DNA_sequence_str) > 0:
    print("The first DNA sequence is")   
    for j in range(0,len(first_DNA_sequence_str),100):
        DNA_string_100_per_line_str = first_DNA_sequence_str[j:j+100]
        print(DNA_string_100_per_line_str)
    print("The length of the first sequence is %s bps." % len(first_DNA_sequence_str))
    print("\n")
    
else:
    print("The fisrt DNA sequence not found.")
    print("Please reinput the first DNA sequence in the txt file. Thank you!")
    sys.exit()
    

#get 2nd DNA sequence and print
second_DNA_sequence_str = str(SeqIO.read(second_DNA_file, "fasta").seq)

if len(second_DNA_sequence_str) > 0:
    print("The second DNA sequence is")   
    for j in range(0,len(second_DNA_sequence_str),100):
        DNA_string_100_per_line_str = second_DNA_sequence_str[j:j+100]
        print(DNA_string_100_per_line_str)
    print("The length of the second sequence is %s bps." % len(second_DNA_sequence_str))
    print("\n")
    
else:
    print("The second DNA sequence not found.")
    print("Please reinput the second DNA sequence in the txt file. Thank you!")
    sys.exit()


#get forward gene primers, with length 16-30bps, stored in a list of strs
gene_primers_5_list = []
for i in range(14,31):                                                      #search area can be adjusted if necessary.
    gene_5_str = first_DNA_sequence_str[0:i]      
    gene_primers_5_list.append(gene_5_str)


#get reverse gene primers, with lenght 16-30bps, stored in a list of strs
gene_primers_3_list = []

for i in range(14,31):                                                      #search area can be adjusted if necessary.
    gene_3_str = str(Seq(second_DNA_sequence_str).reverse_complement())[0:i]      
    gene_primers_3_list.append(gene_3_str)

"""
The forward and reverse gene primers have now been stored in lists. Then, we need
to caculate the position of conjunction between the orginial two DNA fragments.
"""

#get ligated DNA sequence = 
ligated_DNA_sequence_str = first_DNA_sequence_str + second_DNA_sequence_str


#nested loops              
primer_sixer_list = []                     #element is primer Seq,Seq,Seq,Seq,Tm_overlap,K
primer_sixers_list = []                    #element is primer_quintet_Seq

p = len(first_DNA_sequence_str)   # int type    
        
Tm_favorite = 63        
for m in range(len(gene_primers_5_list)):
    for n in range(len(gene_primers_3_list)):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        gene5_Seq = Seq(gene_primers_5_list[m])
                        gene3_Seq = Seq(gene_primers_3_list[n])


                        overlap_for_left_str = ligated_DNA_sequence_str[p-9-i:p]
                        overlap_for_right_str = ligated_DNA_sequence_str[p:p+15+j]
                        overlap5_Seq = Seq(overlap_for_left_str + overlap_for_right_str)

                        overlap_rev_left_str = ligated_DNA_sequence_str[p-15-k:p]
                        overlap_rev_right_str = ligated_DNA_sequence_str[p:p+9+l]
                        overlap3_Seq = Seq(overlap_rev_left_str + overlap_rev_right_str).reverse_complement()

                        overlap_Seq = Seq(ligated_DNA_sequence_str[p-9-i:p+9+l])
                        
                        Tm_primer5 = round(MeltingTemp.Tm_NN(gene5_Seq),2)
                        Tm_primer3 = round(MeltingTemp.Tm_NN(gene3_Seq),2)
                        Tm_overlap5 = round(MeltingTemp.Tm_NN(overlap5_Seq),2)
                        Tm_overlap3 = round(MeltingTemp.Tm_NN(overlap3_Seq),2)
                        Tm_overlap = round(MeltingTemp.Tm_NN(overlap_Seq),2)
        

                
                        cnd1 = abs(Tm_primer5 - Tm_overlap3) 
                        cnd2 = abs(Tm_overlap5 - Tm_primer3)
                        cnd3 = abs(Tm_primer5 - Tm_primer3)
                        cnd4 = abs((Tm_primer5 + Tm_primer3)/2 - Tm_overlap - 3)
                        cnd5 = abs(Tm_overlap - Tm_favorite)
                        cnd6 = abs((Tm_primer5 + Tm_primer3)/2 - Tm_favorite)
                        #here set difference = 3, but can be changed if needed

                         
                        primer_sixer_list.append(gene5_Seq)
                        primer_sixer_list.append(gene3_Seq)
                        primer_sixer_list.append(overlap5_Seq)
                        primer_sixer_list.append(overlap3_Seq)
                        primer_sixer_list.append(Tm_overlap)
                        
 
                        K = 100*(cnd1 + cnd2 + cnd3 + cnd4 + cnd5 + cnd6)
                        K = int(K)
                        primer_sixer_list.append(K)
                        
                        primer_sixers_list.append(primer_sixer_list)
                        primer_sixer_list = []


#sort primer_sixers_list by K value first then by the overlap Tm value.
primer_sixers_list_sorted = sorted(primer_sixers_list, key = itemgetter(5,4))

#print the original primer list
print("The number of original primers is %s." %len(primer_sixers_list))
print("The original primers are shown in the following")
for i in range(len(primer_sixers_list)):
    for j in range(6):
        print(primer_sixers_list[i][j])
        
#print the sorted primer list
print("The number of sorted primers is %s." %len(primer_sixers_list_sorted))
print("The sorted primers are shown in the following")
for i in range(len(primer_sixers_list_sorted)):
    for j in range(6):
        print(primer_sixers_list_sorted[i][j])
    
           
"""
This section print each combination of tetrad primers to txt file, plus necessary parameters of each primer,
including name, sequence, length, GC, Tm , overlap Tm and K value.
"""
    
     
file_name_str = sys.argv[1].split('.')[0] + '_' + sys.argv[2].split('.')[0] + "_ligation_primers" + ".txt"
primer_file = open(file_name_str,'w')
print("primer".ljust(20,' '),\
      "sequence".ljust(40,' '),\
      "length".ljust(10,' '),\
      "GC".ljust(10,' '),\
      "Tm".ljust(10,' '),\
      "Overlap_Tm".ljust(15,' '),\
      "K Value".ljust(10,' '),\
      sep="",\
      file=primer_file
      )

for i in range(len(primer_sixers_list_sorted)):
    for j in range(4):
        print_primer_sequence_str = primer_sixers_list_sorted[i][j]
            
        if j == 0:
            print("\n",file=primer_file)
            S0_print_primer_line_str = "gene_5"
            S5_print_primer_line_str = ""
        elif j == 1:
            S0_print_primer_line_str = "gene_3"
            S5_print_primer_line_str = ""
        elif j == 2:
            S0_print_primer_line_str = "overlap_5"
            S5_print_primer_line_str = str(primer_sixers_list_sorted[i][4])
        else:
            S0_print_primer_line_str = "overlap_3"
            S5_print_primer_line_str = str(primer_sixers_list_sorted[i][4])

        S1_print_primer_line_str = str(print_primer_sequence_str)
        S2_print_primer_line_str = str(len(S1_print_primer_line_str))
        S3_print_primer_line_str = str(round(GC(S1_print_primer_line_str),2))
        S4_print_primer_line_str = str(round(MeltingTemp.Tm_NN(S1_print_primer_line_str),2))
        
        S6_print_primer_line_str = str(primer_sixers_list_sorted[i][5])

        print(S0_print_primer_line_str.ljust(20,' '),\
              S1_print_primer_line_str.ljust(40,' '),\
              S2_print_primer_line_str.ljust(10,' '),\
              S3_print_primer_line_str.ljust(10,' '),\
              S4_print_primer_line_str.ljust(10,' '),\
              S5_print_primer_line_str.ljust(15,' '),\
              S6_print_primer_line_str.ljust(10,' '),\
              sep="",\
              file=primer_file
              )
            
primer_file.close()
    




