BC_no_log('Files\Phenotype_data') % box-cox transformation of phenotype data

extractStrainsGenotypes('your_genome_file'); % extracts individual strains from your master genomic data
mergeStrainFiles('Phenotype_data','output_genome',0,1,0,0); % merges a genome file representative of your sample
filterSNP_RS_Vals('output_genome'); % creates 2 files for each chromosome, one with only the SNPs in that chromosome and one containing all other SNPs

phenotypes = fastLMM('output_genome', 'untransformed', 'Phenotype_data', [1:5], 0); % perform GWAS using leave-one-chromosome-out linkage disequilibrium contamination control, on phenotypes [1:5]; can take coviarate as optional last argument
mergeResults('untransformed',phenotypes,1,0,0); % merges the results from each chromosome into a master file, calculating Q-values / false-discovery rate statistics, and generates summary files across phenotypes for lowest Ps/Qs

phenotypes = fastLMM('output_genome', 'bc_transformed', 'Phenotype_data_BCtrans', [1:5], 0);  % same opertions, now on our box-cox transformed data
mergeResults('bc_transformed',phenotypes,1,0,0); % same opertions, now on our box-cox transformed data
