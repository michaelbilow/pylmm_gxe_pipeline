Instructions for running this:

The ./raw and ./raw_pheno folders contain the original data from JaeHoon

git clone https://github.com/michaelbilow/pylmm_zarlab (commit 7bf525c)

pip install numpy scipy pandas matplotlib --user

1) Preprocess the data into the following form:
  - Plink study.bed, stud.bim, study.fam
  - study.pheno
    - Phenotypes as columns
    - Also, IID, FID, and sex to match the fam file
    - Covariates starting with "Cov"
    - Envirnoment as column named Env
    - Everything should be coded numerically
  - put the .pheno, .bim, .bed, .fam in the raw_data folder


2) Run plink (run_plink.py) on that output, send to the clean_plink 

3) Quantile normalize the phenotypes within groups and between groups (quantile_normalize.py)

5a) Clean and separate the phenotypes, each into its own file (gcta_prepare.py)

5b) Create .gxe files for GCTA to compute variance components (gcta_prepare.py)

5c) Move the files from 5a and 5b to the appropriate folders (gcta_prepare.py)

6) run gcta to make the .grm matrix and .hsq variance components (run_gcta.py)

7) Generate the genetic kinship matrix and the GxE kinship from the grm files (grm_to_kin.py)

7a) Also generate a kinship matrix with pylmm?

8) Generate a blank kinship matrix (blank_kinship.py)

9) Generate a mixed Genetic-GxE kinship matrix using the .hsq files from GCTA (make_2RE_kinship.py)

9) Run pylmm (Pylmm version from pylmm_zarlab, commit 7bf525c)) for GxE using the kinship matrix (run_pylmm.py)

10) run analyze_pylmm_results on the output to get the lambdas

11) run read_hsq to collect the variance components

