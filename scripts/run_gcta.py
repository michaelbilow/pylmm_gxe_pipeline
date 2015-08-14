import subprocess
from os.path import join, splitext, split, exists
from os import listdir, makedirs
import shutil


GCTA_LOC = '../programs/gcta_mac'


def make_grm(plink_fn):
    make_grm_command = "{gcta_loc} --bfile {plink} --make-grm --autosome --out {plink}"
    this_grm_command = make_grm_command.format(gcta_loc=GCTA_LOC, plink=plink_fn)
    subprocess.call(this_grm_command, shell=True)
    return


def estimate_variance_components_wrapper(args):
    return estimate_variance_components(*args)


def estimate_variance_components(plink_fn, pheno_fn, gxe_fn):
    gxe_command = "{gcta_loc} --reml --reml-maxit 10000 --reml-alg 2 --reml-no-lrt --grm {plink} --pheno {pheno} \
            --gxe {gxe} --reml-lrt 2 --out {pheno}"
    this_gxe_command = gxe_command.format(plink=plink_fn, pheno=pheno_fn, gxe=gxe_fn, gcta_loc=GCTA_LOC)
    subprocess.call(this_gxe_command, shell=True)
    return


def main(project_folder, study_name):
    plink_folder = join(project_folder, 'clean_plink')
    input_plink_fn = join(plink_folder, study_name)
    make_grm(input_plink_fn)

    input_gxe_fn = join(plink_folder, '{}.gxe'.format(study_name))
    pheno_folder = join(project_folder, 'clean_pheno')
    pheno_fns = [join(pheno_folder, x) for x in listdir(pheno_folder)]

    args_list = [(input_plink_fn, fn, input_gxe_fn) for fn in pheno_fns]
    for args in args_list:
        estimate_variance_components_wrapper(args)

    hsq_folder = join(project_folder, 'hsq')
    if not exists(hsq_folder):
        makedirs(hsq_folder)
    for x in listdir(pheno_folder):
        if '.hsq' in x:
            shutil.move(join(pheno_folder, x),
                        join(hsq_folder, x))
    return


if __name__ == "__main__":
    proj_folder = '../eqtl'
    study = 'HAEC'
    main(proj_folder, study)