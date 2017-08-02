This set of MATLAB code creates a pipeline for performing a genome-wide association scan (GWAS) on phenotype data while controlling for genetic relatedness using a random effect of genetic background.  Here the implementation is geared towards GWAS performed on common and recombinant inbred mouse strain phenotypes, including behavioral and mRNA expression (i.e., eQTL) phenotypes, such that strain genetic background is included the model as a random effect.  Models are fitted using FaST-LMM (Factored Spectrally Transformed Linear Mixed Model; see refs). The leave-one-chromosome-out approach is taken to mitigate linkage disequilibrium contamination of results while computing the relatedness matrix.  This pipeline is made for analysis on Windows but could easily be adapted for Unix platforms or to use the (more frequently updated) Python version of FaST-LMM.  The results can be visualized using Haploview, and  summary output files providing the best results (i.e, lowest P-values, and computed false-discovery rate Q-values for each phenotype) are produced.

These scripts require installation of Windows compiled FaST-LMM to C:\Fast-LMM and PLINK to C:\Fast-LMM\Plink.  
An example phenotype file format is given in Files\Example_data.txt, and an auxiliary excel macro written in VBA to output this file from a master excel data file is given in Files\Example_data.xlsm (note regarding macro: read code carefully, this macro is customized to a dataset I used and specifically excludes certain subjects from the output).  See the FaST-LMM reference for detail regarding the appropriate format of input phenotype files, genotype files, covariate files, etc. Note that as per FaST-LMM guidelines, a value of -9 in a phenotype file indicates missing data.

An example analysis pipeline is given in run_gwas.m; the general processes involves the following functions in the following order:

BC_no_log() - (optional) box-cox transforms phenotype data
extractStrainsGenotypes() - extracts individual strains from your master PLINK genomic data
mergeStrainFiles() - creates/merges a genome file representative of your sample
filterSNP_RS_Vals() - creates 2 files for each chromosome, one with only the SNPs in that chromosome and one containing all other SNPs
fastLMM() - perform GWAS by calling FaST-LMM, using leave-one-chromosome-out linkage disequilibrium contamination control
mergeResults()  - merges the results from each chromosome into a master file, calculating Q-values / false-discovery rate statistics, and generates summary files across phenotypes for lowest Ps/Qs

Functions are described in greater detail in the header of the function, including explanations of all arguments and outputs; additional auxiliary functions are present which are also described headers.


References: 

FaST linear mixed models for genome-wide association studies
Christoph Lippert, Jennifer Listgarten,	Ying Liu, Carl M Kadie,	Robert I Davidson & David Heckerman
Nature Methods 8, 833–835 (2011) doi:10.1038/nmeth.1681

Plink, Whole genome association analysis toolset, http://zzz.bwh.harvard.edu/plink/faq.shtml