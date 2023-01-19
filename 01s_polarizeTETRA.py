#  
# Luis Leal (2021)
#
# Written in Python 3



### Script used to polarize allotetraploid sequence in a MSA.
### It is assumed that sequences in input MSA are aligned.
###

### Input arguments
### 1. MSA (fasta format)
### 2. Number of species (including outgroup)
### 3. Reference sequence ID
### 4. Allopolyploid sequence ID



#print('\n Parsing started ...')




######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time
from collections import Counter
from pdb import set_trace as bp                           # when using break points during debugging >> usage: bp()






######################################################## READ INPUT ARGUMENTS


try:
    fhand = open(sys.argv[1], 'r')                  			# input fasta file
except:
    print('\n Error: input file missing.')
    exit() 


try:
    Nspecies = sys.argv[2]                  				# Number of species
except:
    print('\n Error: provide the following arguments:')
    print('01s_polarizeTETRA.py [input.fasta] [No. species] [refSequence ID] [polyploid ID] \n')
    exit() 


try:
    RefSeqID = sys.argv[3]                  				# Reference sequence ID (polarizing sequence)
except:
    print('\n Error: provide the following arguments:')
    print('01s_polarizeTETRA.py [input.fasta] [No. species] [refSequence ID] [polyploid ID] \n')
    exit() 


try:
    POLYPID = sys.argv[4]                  				# Allopolyploid sequence ID  
except:
    print('\n Error: provide the following arguments:')
    print('01s_polarizeTETRA.py [input.fasta] [No. species] [refSequence ID] [polyploid ID] \n')
    exit() 















### AUXILIARY FUNCTION: code heterozygotic sites according to IUPAC nomenclature
### See also: http://www.iqtree.org/doc/Frequently-Asked-Questions#how-does-iq-tree-treat-gapmissingambiguous-characters

def extendedCODE(input_list):
   species_nuc_uniq = list(set(input_list))  # get unique nucleotides
   species_nuc_uniq.sort()		     # sort list alphabetically

   outvalue = ''
   if species_nuc_uniq == ['A'] : outvalue = 'A'
   if species_nuc_uniq == ['C'] : outvalue = 'C'
   if species_nuc_uniq == ['G'] : outvalue = 'G'
   if species_nuc_uniq == ['T'] : outvalue = 'T'
   if species_nuc_uniq == ['A','C'] : outvalue = 'M'
   if species_nuc_uniq == ['A','G'] : outvalue = 'R'
   if species_nuc_uniq == ['A','T'] : outvalue = 'W'
   if species_nuc_uniq == ['C','G'] : outvalue = 'S'
   if species_nuc_uniq == ['C','T'] : outvalue = 'Y'
   if species_nuc_uniq == ['G','T'] : outvalue = 'K'
   if species_nuc_uniq == ['A','C','G'] : outvalue = 'V'
   if species_nuc_uniq == ['A','C','T'] : outvalue = 'H'
   if species_nuc_uniq == ['A','G','T'] : outvalue = 'D'
   if species_nuc_uniq == ['C','G','T'] : outvalue = 'B'
   if species_nuc_uniq == ['A','C','G','T'] : outvalue = 'N'

   return outvalue


### end function




### AUXILIARY FUNCTION: expand IUPAC codes

def unfoldCODE(species_nuc):

   outvalue = ''
   if species_nuc == 'A' : outvalue = ['A']
   if species_nuc == 'C' : outvalue = ['C']
   if species_nuc == 'G' : outvalue = ['G']
   if species_nuc == 'T' : outvalue = ['T']
   if species_nuc == 'M' : outvalue = ['A','C']
   if species_nuc == 'R' : outvalue = ['A','G']
   if species_nuc == 'W' : outvalue = ['A','T']
   if species_nuc == 'S' : outvalue = ['C','G']
   if species_nuc == 'Y' : outvalue = ['C','T']
   if species_nuc == 'K' : outvalue = ['G','T']
   if species_nuc == 'V' : outvalue = ['A','C','G']
   if species_nuc == 'H' : outvalue = ['A','C','T']
   if species_nuc == 'D' : outvalue = ['A','G','T']
   if species_nuc == 'B' : outvalue = ['C','G','T']
   if species_nuc == 'N' : outvalue = ['N']
   if outvalue == '' :
      print('\n Error: sequences contain non-IUPAC nucleotide character\n')
      exit()

   return outvalue

### end function










######################################################## Read MSA

INfasta_species = list()
INfasta_seq = list()
FLAG_aux = 0
for line in fhand:
    aux_NC = line.replace("\n", "")
    if FLAG_aux == 0 :
      INfasta_species.append(aux_NC)
      FLAG_aux = 1
    else :
      INfasta_seq.append(aux_NC)
      FLAG_aux = 0


INfasta_species_CLEAN = INfasta_species
INfasta_seq_CLEAN = INfasta_seq


# save species name and associated sequence to dictionary
INfasta_dict = dict()
for i in range(len(INfasta_species_CLEAN)) :
   INfasta_dict[INfasta_species_CLEAN[i]] = INfasta_seq_CLEAN[i]

#for i in range(len(INfasta_species_CLEAN)) :
#   print(INfasta_species_CLEAN[i], INfasta_dict[INfasta_species_CLEAN[i]])

#print(len(INfasta_dict))



# species list
species_list = list()
for tag in INfasta_species_CLEAN :
   tag = tag.replace(">", "")                   #remove '>' from sequence name
   species_list.append(tag)	

species_list = list(set(species_list))		# unique elements only


# Check number of sequences
DIC_LEN=len(INfasta_dict)
if DIC_LEN != int(Nspecies) :
   print('\n Error: check number of sequences in fasta file\n')
   exit()

if len(species_list) != int(Nspecies) :
   print('\n Error: check number of sequences in fasta file\n')
   exit()



# Check all sequences have the same length
FLAG_aux = 0
for i in range(len(INfasta_species_CLEAN)) :
   len_seq = len(INfasta_dict[INfasta_species_CLEAN[i]])
   #print(len_seq)
   if FLAG_aux == 0 :
      len_ref = len_seq
      FLAG_aux = 1
   else :
      if len_seq != len_ref :
         print('\n Error: sequences have different lenghts\n')
         exit()

          






###################################################### Unfold allotetraploid and ref sequences

## Reference sequence
REFG_ID = '>' + RefSeqID
REFseq = list()
for i in INfasta_dict[REFG_ID] :
   nuc_ste = unfoldCODE(i)
   REFseq.append(nuc_ste)
   #print(i,nuc_ste)



## Allotetraploid sequence
POLY_ID = '>' + POLYPID
TETRAseq = list()

for i in INfasta_dict[POLY_ID] :
   nuc_ste = unfoldCODE(i)
   TETRAseq.append(nuc_ste)
   #print(i,nuc_ste)


## Add allotetraploid and ref sequence to new dictionary
DICT_seq = dict()	
DICT_seq[POLYPID] = TETRAseq
DICT_seq[RefSeqID] = REFseq


#for key, value in DICT_seq.items():
#   print(key,len(DICT_seq[key]))











###################################################### Create polarized sequence for ALLOTETRAPLOID


DICT_POLAR_ALT = dict()	        # dictionary containing polarized allopolyploid sequence (also all other sequences, but which stay as in the input MSA)


### Add non-allopolyploid sequences to dictionary
species_list_aux = species_list[:]	# copy species list
species_list_aux.remove(POLYPID)	# remove allopolyploid

for item in species_list_aux :
   aux_ID = '>' + item
   DICT_POLAR_ALT[item] = INfasta_dict[aux_ID]



### Allopolyploid

ALT_list = list()
species_seq = DICT_seq[POLYPID]
for i in range(len(REFseq)) :		# for each nucleotide in reference sequence
   ref_nuc_ls = REFseq[i]
   species_nuc = list(species_seq[i])
   #print(ref_nuc_ls, species_nuc)
   #
   if species_nuc[0] == 'N' :		# if polyploid is masked
      aux_nuc = 'N' 
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc) 
   elif ref_nuc_ls[0] == 'N' :		# if reference is masked
      aux_nuc = 'N' 
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif (species_nuc == ['A'] or species_nuc == ['C'] or species_nuc == ['G'] or species_nuc == ['T']) :		# polyploid if fixed
      ALT_list.append(species_nuc[0])
      #print(ref_nuc_ls, species_nuc, species_nuc[0])
   elif species_nuc == ref_nuc_ls :						# if reference alleles were found in the polyploid and there are no other alleles present
      aux_nuc = extendedCODE(species_nuc)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif all(x0 in species_nuc for x0 in ref_nuc_ls) :				# if all reference alleles were found also in polyploid, but there are other alleles present as well
      aux_nuc0 = [pp for pp in species_nuc if pp not in ref_nuc_ls]		# remove reference alleles from polyploid allele set
      if aux_nuc0 == [] : aux_nuc0 = species_nuc
      aux_nuc = extendedCODE(aux_nuc0)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif any(x0 in ref_nuc_ls for x0 in species_nuc) :				# if at least one of the reference alleles is found also in the polyploid
      aux_nuc0 = [pp for pp in species_nuc if pp not in ref_nuc_ls]		# remove reference alleles from polyploid allele set, if present
      if aux_nuc0 == [] : aux_nuc0 = species_nuc
      aux_nuc = extendedCODE(aux_nuc0)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   else :									# reference allele not found in polyploid
      aux_nuc = extendedCODE(species_nuc)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)



# update dictionary
DICT_POLAR_ALT[POLYPID]=''.join(ALT_list)

#for key, value in DICT_POLAR_ALT.items():
#   print(key,len(DICT_POLAR_ALT[key]))





## Count differences between polarized and non-polarized allotetraploid sequence

tetra_RAW = INfasta_dict[POLY_ID]
tetra_ALT = DICT_POLAR_ALT[POLYPID]

counter_ALT_RAW = 0

for i in range(len(tetra_RAW)) :		
   hap_RAW = tetra_RAW[i]
   hap_ALT = tetra_ALT[i]
   if hap_RAW != hap_ALT :
      counter_ALT_RAW += 1
      #print(hap_RAW,hap_ALT)


#print(counter_ALT_RAW)







##### Save polarized MSA to file

straux1=str(sys.argv[1])
straux1 = straux1[:(len(straux1)-3)]
outfile_ALT = straux1 + '-ALT.fasta'
outfile_stats1 = straux1 + '.stats'

outfile2 = open(outfile_ALT, 'w')   
for key, value in DICT_POLAR_ALT.items() :
    outfile2.write('>' + key)
    outfile2.write('\n')
    outfile2.write(value)
    outfile2.write('\n')

outfile3 = open(outfile_stats1, 'w')   
outfile3.write('Number of nucleotide differences in polyploid sequence (before vs after polarization):')
outfile3.write('\n')
outfile3.write(str(counter_ALT_RAW))
outfile3.write('\n')

outfile2.close()
outfile3.close()



 
      
