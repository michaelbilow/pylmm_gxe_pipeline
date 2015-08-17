import subprocess
from os.path import join
import sys

MAF = .1
GENO = .02
MAC_PLINK_LOC = 'plink1.07_mac'
LINUX_PLINK_LOC = 'plink1.07linux'


def plink_run(input_fn, output_fn, maf=MAF, geno=GENO):
    if sys.platform == 'darwin':
        plink_loc = MAC_PLINK_LOC
    elif sys.platform == 'linux2':
        plink_loc = LINUX_PLINK_LOC
    else:
        raise ValueError

    command = "{} --bfile {} --autosome --maf {} --geno {} --make-bed --out {}"
    command = command.format(plink_loc, input_fn, maf, geno, output_fn)
    print command
    subprocess.call(command, shell=True)


def main(input_folder, output_folder, study_name, kwargs=False):
    print('Running Plink on ')
    input_fn = join(input_folder, study_name)
    output_fn = join(output_folder, study_name)
    if kwargs:
        plink_run(input_folder, output_folder, **kwargs)
    else:
        plink_run(input_fn, output_fn)


if __name__ == "__main__":
    # expanded_hmdp_fn = '../hmdp/expanded/hmdp'
    # expanded_hmdp_output_fn = '../hmdp/clean_plink/hmdp'
    #
    # plink_run(input_fn=expanded_hmdp_fn, output_fn=expanded_hmdp_output_fn)
    #
    # haec_fn = '../eqtl/raw_plink/HAEC'
    # haec_output_fn = '../eqtl/clean_plink/HAEC'
    # plink_run(input_fn=haec_fn, output_fn=haec_output_fn)
    pass