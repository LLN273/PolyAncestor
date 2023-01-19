# GenePolar
Pipeline used to identify the phylogenetic position of the parental species of a tetraploid species.

Citation:
Leal J.L., Milesi P., Salojärvi J., Lascoux M. 2022. Phylogenetic Analysis of Allotetraploid Species Using Polarized Genomic Sequences. *Syst. Biol.* (forthcoming)

____

# TUTORIAL

In this tutorial, we provide a step-by-step guide on how to use genomic polarization [1] to identify the phylogenetic position of the two parental species of an allotetraploid in an iterative process. We provide an example for a phylogeny containing fifteen diploid species and one allotetraploid species.

It is assumed you are working in a Unix or macOS environment.

## 1. Software required 

Genomic polarization script, which you can download from our GitHub page [here](https://github.com/LLN273/GenePolar/blob/main/01s_polarizeTETRA.py).

IQ-TREE2 [2] for performing phylogenetic inference for each locus: [http://www.iqtree.org/](http://www.iqtree.org/)

ASTRAL [3] for estimating the species tree: [https://github.com/smirarab/ASTRAL](https://github.com/smirarab/ASTRAL)

Newick Utilities [4] for contracting branches with low bootstrap support: [https://github.com/tjunier/newick_utils](https://github.com/tjunier/newick_utils)

You also need to have [Python3](https://www.python.org/) and [R](https://www.r-project.org/) installed.

## 2. Input data

We will be using a small dataset that contains multiple sequence alignments (MSAs) for 50 loci that have evolved under moderate levels of incomplete lineage sorting (ILS). Each MSA contains nucleotide sequences for fifteen diploid species and one allotetraploid species ("T-9-14") whose ancient parental species, now extinct, are phylogenetically located near species "10" and "16". Each species included in the MSA, including the allotetraploid, is represented by a single nucleic acid sequence. Heterozygotic sites are expected to be coded according to the IUPAC nomenclature *but the input data does not need to be pre-phased.* It is assumed sequences are aligned.

You can download the test data [here](https://github.com/LLN273/GenePolar/tree/main/aux/MSA_test_data.zip). 

Note: When processing your own data, it is expected that you will have MSAs for hundreds or even thousands of loci. There are many ways to generate MSAs based on a VCF file. We provide some scripts that can be of help [here](https://github.com/LLN273/Genomic_polarization_allotetraploids/tree/main/02_Generate_MSAs).
If you use our pipeline, notice that we use separate VCF files for each sample, and that we also have separate VCF files for SNPs and indels. Also, it will help if you have an annotated reference genome but this is not strictly necessary. See our paper [1] for a more detailed discussion and suggestions on what to do if you don't have a suitable reference genome.

## 3. Polarize allotetraploid sequence (first iteration)

Copy the genomic polarization script (01s_polarizeTETRA.py) to the folder where you saved the MSAs. 

This is the general command used to polarize the allotetraploid sequence: 

```bash
python3 01s_polarizeTETRA.py [input.fasta] [No. species] [refSequence ID] [polyploid ID]
```

Input arguments:
- [input.fasta]: MSA in fasta format. Sequences are assumed to be aligned and coded using IUPAC nomenclature. 
  If your MSA includes gaps (-), recode them as masked sites (N) before running this script.                   
- [No. species]: Number of species in MSA.
- [refSequence ID]: Reference sequence ID (polarizing sequence), as shown in MSA.
- [polyploid ID]: Allopolyploid sequence ID, as shown in MSA.

NOTE: When running this script over a very large number of loci, you may occasionally get the following error message "Error: check number of sequences in fasta file". This happens if some of your MSAs are missing some sequences. This could have happened for a variety of reasons. For instance, trimAl, a program used to clean MSAs, will remove any sequence that contains only masked sites (NNNNN). If only a few of your MSAs are missing sequences, you can ignore the error messages you get when running the polarization script.

During the first iteration, the reference sequence used to polarize the allotetraploid ("T-9-14") is selected randomly among the diploid species included in the MSA. For this example, we will use species "11" as the reference sequence. There are sixteen sequences in the MSA. For locus 0001, the command used to polarize the allopolyploid is therefore:

```bash
python3 01s_polarizeTETRA.py locus_0001.fa 16 "11" "T-9-14"
```

To polarize all 50 loci using one single command, one can use the following code:

```bash
for k in `seq -w 01 1 50`; do
   python3 01s_polarizeTETRA.py locus_00${k}.fa 16 "11" "T-9-14"
done
```

This should take just a few seconds to finish.

The command above produces 50 new MSAs (files ending in "-ALT.fasta") where the allopolyploid sequence has been polarized using species "11" as the reference sequence.

Create a new folder named "Iteration_1" and move the polarized MSAs into it;

```bash
mkdir Iteration_1
cd Iteration_1
mv ../*-ALT.fasta .
```

## 4. Phylogenetic inference for each individual locus

During the next step, we use IQ-TREE2 [2], a widely used tree inference software package based on the maximum-likelihood criterion, to infer phylogenetic trees for each locus. Running IQ-TREE2 properly and efficiently requires some familiarity with its many options. Here, and in order to speed up things, we will run it with (mostly) default options. For more advanced usage, check out the tutorials available on the [IQ-TREE2 website](http://www.iqtree.org/doc/) or our [scripts](https://github.com/LLN273/Genomic_polarization_allotetraploids/tree/main/04_Phylogenetic_Analysis) on GitHub.

The following command will run IQ-TREE2 for all 50 loci, using a "HKY+F+G4" model, 1000 ultrafast bootstrap replicates, and setting species "0" as the outgroup. You must provide the path to the folder where iqtree has been installed in your computer. This should take 3-4 minutes to finish.

```bash
for k in `seq -w 01 1 50`; do
      /path/to/iqtree/iqtree \
      -bb 1000 \
      -s locus_00${k}-ALT.fasta \
      -m "HKY+F+G4" \
      -o "0"
	  
done
```

IQ-TREE2 produces several output files. The ones of interest to us have the ".iqtree" extension and contain the inferred consensus trees. You should have 50 such files.

## 5. Generate species tree

We will use ASTRAL [3], a super-tree inference method statistically consistent under the multi-species coalescent model, to produce a species tree based on the single-loci phylogenies produced above.

Before running ASTRAL, we must prepare the input file, which is a file containing all the consensus trees inferred using IQ-TREE2. For the current example, this file should contain 50 gene trees.

You can generate ASTRAL's input file using the following command:

```bash
rm -f astral_infile.txt
for k in `seq -w 01 1 50`; do
   cat locus_00${k}-ALT.fasta.iqtree | grep -A 2 'Consensus tree in newick format:' | tail -1 >> astral_infile.txt
done
```

This produces a file called "astral_infile.txt ", which in our case contains 50 lines, one for each gene tree.

It is usually recommended for branches with very low bootstrap support to be contracted before running ASTRAL. This can be done using Newick Utilities. In the example below, we will contract branches with less than 30% bootstrap support (you must update the path to Newick Utilities).

```bash
/path/to/newick_utils/src/nw_ed astral_infile.txt 'i & b<=30' o > astral_infile_BS30.newick
```

Finally, we run ASTRAL in exact mode (-x; meaning that the entire tree space is tested. This option should only be used if the number of species included in the phylogeny is relatively small, e.g. <20). We will also request an output tree with branch quartet support values (-t 8). This should take a couple of minutes.

```bash
java -jar /path/to/astral/astral.5.7.8.jar \
              -i astral_infile_BS30.newick \
              -o astral_outfile_BS30.quartet.t8.newick \
              -t 8 \
              -x
```

You can visualize the output tree using [TreeGraph](http://treegraph.bioinfweb.info) or [FigTree](http://tree.bio.ed.ac.uk/software/Figtree/). The inferred species tree after the first iteration is shown in Fig. 1. The polarized allotetraploid pairs with the clade containing species "3" and "10".

<img src="https://github.com/LLN273/GenePolar/tree/main/aux/Fig_1_tutorial.png" width="500" />

Figure 1 | Species tree inferred after the first iteration.

We can also check each individual gene tree (IQ-TREE2 output) and compute the frequency with which the polarized allotetraploid pairs with any of the diploid species included in the phylogeny. We use the R script "05_Polyploid_pairing_analysis.R" to perform this analysis (available [here](https://github.com/LLN273/Genomic_polarization_allotetraploids/tree/main/04_Phylogenetic_Analysis)). The general command to run this script is:

```bash
Rscript --no-save 05_Polyploid_pairing_analysis.R [input.file] [polyploid ID] [outgroup ID] [output folder] [output file]
```

Copy this script to your working folder and then run:

```bash
Rscript --no-save 05_Polyploid_pairing_analysis.R \
                            "astral_infile.txt" \
                            "T-9-14" \
                            "0" \
                             ./ \
                            "sister_ID_analysis_GENE_TREES.txt"
```

The output file reveals that for 56% of all gene trees the polarized allotetraploid pairs either with species "10" (32%) or species "3" (24%):

```bash
cat sister_ID_analysis_GENE_TREES.txt
```

| Species   |	Frequency	| Normalized frequency |
| --------- | --------- | -------------------- |
| 10		    | 16			  | 0.32                 |
| 3		      | 12			  | 0.24                 |
| 10_3		  | 7			    | 0.14                 |
| 16		    | 5			    | 0.1                  |
| 4		      | 4			    | 0.08                 |
| 10_3_4	  | 1			    | 0.02                 |
| 10_3_4_6	|1			    | 0.02                 |
| 10_4	    | 1				  | 0.02                 |
| 13	      | 1				  | 0.02                 |
| 13_15	    | 1				  | 0.02                 |
| 13_15_16	| 1			    | 0.02                 |


## 6. Second iteration

During the second iteration, one of {"3", "10"} is used as the new reference sequence. In the present example species "10" gets the most support based on the analysis of individual gene trees (as seen in the table above). However, for didactic purposes we will select species "3" instead. Go back to the folder where you saved your raw (unpolarized) MSAs and run the following command (notice the reference sequence ID has been updated):

```bash
for k in `seq -w 01 1 50`; do
   python3 01s_polarizeTETRA.py locus_00${k}.fa 16 "3" "T-9-14"
done
```

Move the polarized MSAs to a new folder ("Iteration_2") and then repeat steps 4 and 5. The resulting species tree is shown in Figure 2. The polarized allotetraploid now pairs with species "16", sister species to one of the polyploid's parental species.

<img src="https://github.com/LLN273/GenePolar/tree/main/aux/Fig_2_tutorial.png" width="500" />

Figure 2 | Species trees inferred after the first two iterations.

The polyploid pairing frequencies after the second iteration are:

| Species	 | Frequency | Normalized frequency  |
| -------- | --------- | --------------------- |
| 16		   | 30			   | 0.6                   |
| 13_15		 | 9			   | 0.18                  |
| 13_15_16 | 5			   | 0.1                   |
| 15		   | 3			   | 0.06                  |
| 15_16	   | 2	       | 0.04                  |
| 13	     | 1	       | 0.02                  |

## 7. Third iteration

During the third iteration, species "16" is used as the reference sequence. In the resulting phylogeny, the allotetraploid pairs with species "10" (Fig. 3), sister species to the second ancient parental species.

<img src="https://github.com/LLN273/GenePolar/tree/main/aux/Fig_3_tutorial.png" width="800" />

Figure 3 | Species trees inferred after the first three iterations

The polyploid pairing frequencies after the third iteration are:

| Species	 | Frequency | Normalized frequency  |
| -------  | --------- | --------------------- |
| 10		   | 21			   | 0.42                  |
| 3		     | 12			   | 0.24                  |
| 10_3		 | 7			   | 0.14                  |
| 4		     | 5			   | 0.1                   |
| 10_3_4	 | 2			   | 0.04                  |
| 10_3_4_6 | 2			   | 0.04                  |
| 10_4		 | 1			   | 0.02                  |

## 8. Fourth iteration

During the fourth iteration, species "10" was used as the reference sequence. The polarized allotetraploid pairs again with species "16" (Fig. 4), signaling that the iterative procedure has converged. The analysis indicates that the two species closest to the allopolyploid's ancient parental species are species "10" and "16" (from among those included in the analysis).

<img src="https://github.com/LLN273/GenePolar/tree/main/aux/Fig_4_tutorial.png" width="1000" />

Figure 4 | Species trees inferred after the first four iterations. The location of the first parental species was found after two iterations (near species '16') while convergence is attained by the fourth iteration (iteration 5 would have produced a result similar to that shown for iteration 3).


The polyploid pairing frequencies after the fourth iteration are:

| Species	 | Frequency | Normalized frequency |
| -------- | --------- | -------------------- |
| 16		   | 30			   | 0.6                  |
| 13_15		 | 9			   | 0.18                 |
| 13_15_16 | 5			   | 0.1                  |
| 15		   | 3         | 0.06                 |
| 15_16		 | 2			   | 0.04                 |
| 13		   | 1			   | 0.02                 |



### References

[1] Leal J.L., Milesi P., Salojärvi J., Lascoux M. 2022. Phylogenetic Analysis of Allotetraploid Species Using Polarized Genomic Sequences. *Syst. Biol.* (forthcoming)

[2] Minh B.Q., Schmidt H.A., Chernomor O., Schrempf D., Woodhams M.D., von Haeseler A., Lanfear R. 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol. Biol. Evol.* 37:1530–1534.

[3] Zhang C., Rabiee M., Sayyari E., Mirarab S. 2018. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. *BMC Bioinf.* 19:153.

[4] Junier, T. and Zdobnov, E. M. 2010. The Newick Utilities: High-throughput phylogenetic tree processing in the UNIX Shell. *Bioinformatics*. 26: 1669–1670.
