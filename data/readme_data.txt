data.mat
	Genotype and flowering time phenotypes obtained from Atwell et. al. (2010).
        For more info: https://www.ncbi.nlm.nih.gov/pubmed/20336072
	
        -Rare variants (MAF < 0.1) are filtered out.
        -Number of SNPs are reduced to 173219 from 214051.

        X             - binary genotype matrix 199 x 173219 (samples vs SNPs)
        Y             - phenotype matrix 199 x 17 (samples vs phenotypes) - contains NaN values
        genotypes     - genotype matrix 214051 x 2 - 1st column = "0" in X, 2nd column = "1" in X
        MAF           - minor allele frequency vector 173219 x 1 = (# of ones in X) / (# of zeros in X)
        phenotypes    - flowering time phenotypes - 1st column : phenotype_name, 2nd column : index in Atwell (2010) phenotypes
        snp           - SNP information matrix 173219 x 2 - 1st column : chromosome, 2nd column : position
        samples       - samples information matrix 199 x 2 - Contains 2 different unique index for each sample

network.mat
    Contains Genomic Membership (GM) SNP-SNP network
        where the SNPs within the same gene are connected as a clique.
                        
        network       - a 173219 x 173219 binary sparse matrix.