BC_no_log('Files\Phenotype_data')
extractStrainsGenotypes('your_genome_file');
mergeStrainFiles('AllPhenotypes_Genotyped','HMDPinst',0,0,0);
filterSNP_RS_Vals('HMDPinst');
phenotypes = fastLMM('HMDPinst', 'untransAll', 'AllPhenotypes_Genotyped', [1:43]); %can take cov as last argument
mergeResults('untransAll',phenotypes,1);
%crossPhenotypeQ(phenotypes,'All');
clear phenotypes;

phenotypes = fastLMM('HMDPinst', 'bctransAll', 'AllPhenotypes_Genotyped_BCtrans', [1:22]); %can take cov as last argument
mergeResults('bctransAll',phenotypes,1);
%crossPhenotypeQ(phenotypes,'All');
clear phenotypes;

mergeStrainFiles('AllPhenotypes_Genotyped_Reached50','HMDPinst_Reached50',0,1,0);
filterSNP_RS_Vals('HMDPinst_Reached50');
phenotypes = fastLMM('HMDPinst_Reached50', 'untrans_Reached50', 'AllPhenotypes_Genotyped_Reached50', [1:55]); %can take cov as last argument
mergeResults('untrans_Reached50',phenotypes,1);
%crossPhenotypeQ(phenotypes,'All');
clear phenotypes;

phenotypes = fastLMM('HMDPinst_Reached50', 'bctrans_Reached50', 'AllPhenotypes_Genotyped_Reached50_BCtrans', [1:33]); %can take cov as last argument
mergeResults('bctrans_Reached50',phenotypes,1);
%crossPhenotypeQ(phenotypes,'All');
clear phenotypes;
