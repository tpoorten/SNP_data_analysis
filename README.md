SNP_data_analysis
=================
Convert vcf file to table format with script: 'vcf_to_genotypeTab_meta.Rscript'
 - usage example: ./vcf_to_genotypeTab_meta.Rscript -i=test.vcf 
 - input: vcf file (e.g. GATK vcf file)
		- optional arguments: 
				'-acgt' : include this flag to output nucleotide genotypes  (default: on)	
				'-meta' : include this flag to output metadata tables (default: on)	
 - output: tab-delimited table with genotypes
			AND tab-delimited tables with genotype metadata:
				AD: allele depth		DP: read depth
				GQ: genotype quality	PL: phred likelihood				
 - the genotypes dataframe can be used in downstream analyses
 - this script is meant to be run as an executable Rscript
 - genotype format:
	 0 represents homozygous reference allele
	 1 represents heterozygous
	 2 represents homozygous alternative allele
	  NA represents missing data


PCA plot with script: 'snps_PCA.Rscript'
 - purpose - Make Principal Components Analysis (PCA) plot from SNP genotypes
 - usage example: ./snps_PCA.Rscript -i=test.GT.tab  OR  ./snps_PCA.Rscript test.GT.tab
 - input: GT table where genotypes are coded 0,1,2, ...
 - output: PCA plot
