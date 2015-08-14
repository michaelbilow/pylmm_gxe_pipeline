import subprocess

PLINK_LOC = '../programs/plink'
MAF = .1
GENO = .02

def plink_run(input_fn, output_fn, maf=MAF, geno=GENO, plink_loc=PLINK_LOC):
    command = "{} --bfile {} --autosome --maf {} --geno {} --make-bed --out {}"
    command = command.format(plink_loc, input_fn, maf, geno, output_fn)
    print command
    subprocess.call(command, shell=True)


if __name__ == "__main__":
    expanded_hmdp_fn = '../hmdp/expanded/hmdp'
    expanded_hmdp_output_fn = '../hmdp/clean_plink/hmdp'

    plink_run(input_fn=expanded_hmdp_fn, output_fn=expanded_hmdp_output_fn)

    haec_fn = '../eqtl/raw_plink/HAEC'
    haec_output_fn = '../eqtl/clean_plink/HAEC'
    plink_run(input_fn=haec_fn, output_fn=haec_output_fn)
