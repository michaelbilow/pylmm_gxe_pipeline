import subprocess
from os.path import join, split, splitext, exists
from os import listdir, makedirs
import multiprocessing
import itertools
import sys

if sys.platform == 'darwin':
    PYLMM_SCRIPT_LOC = '~/anaconda/bin/pylmmGWAS.py'
elif sys.platform == 'linux2':
    PYLMM_SCRIPT_LOC = '/u/home/m/mbilow/anaconda/bin/pylmmGWAS.py'
else:
    raise ValueError

def pylmm_run_wrapper(kwargs):
    return pylmm_run(**kwargs)


def pylmm_run(plinkfile, phenofile, gxe, kfile, noCorrect, covfile, test_type, output_folder):
    # print locals()
    pylmm_command = 'python {script_loc} -v --bfile {plink} ' \
                    '--kfile {kinship} --phenofile {pheno} --covfile {covfile} ' \
                    '--limit {gxe} {test_type} {noCorrect} {outfile}'

    kfile_name = split(kfile)[1]
    kfile_type = 'K-GxE' if 'gxe' in kfile_name \
        else 'K-Blank' if 'blank' in kfile_name \
        else 'K-2RE' if '2re' in kfile_name \
        else 'K-Id' if 'identity' in kfile_name \
        else 'K-Unknown'
    run_type = 'GxE' if gxe else 'G'
    out_fn = '{}/{}_{}.stats_{}'.format(output_folder,
                                        split(phenofile)[1],
                                        kfile_type,
                                        run_type)

    this_command = pylmm_command.format(script_loc=PYLMM_SCRIPT_LOC,
                                        pheno=phenofile,
                                        plink=plinkfile,
                                        kinship=kfile,
                                        gxe='--GxE' if gxe else '',
                                        noCorrect='--noKCorrection' if noCorrect else '',
                                        outfile=out_fn,
                                        covfile=covfile,
                                        test_type=test_type)
    print this_command
    subprocess.call(this_command, shell=True)
    return


def main(input_plink_folder, input_pheno_folder, results_folder, two_RE_kfile_folder, study_name):
    if not exists(results_folder):
        makedirs(results_folder)
    pheno_files = [join(input_pheno_folder, x) for x in listdir(input_pheno_folder)]
    input_plinkfile = join(input_plink_folder, study_name)
    input_covfile = join(input_plink_folder, '{}.cov'.format(study_name))
    kfile_list = [join(input_plink_folder, x)for x in listdir(input_plink_folder) if splitext(x)[1] == '.kin']
    kwargs_list = [{'plinkfile': input_plinkfile,
                    'kfile': kfile,
                    'phenofile': fn,
                    'gxe': True,
                    'noCorrect': True,
                    'output_folder': results_folder,
                    'covfile': input_covfile,
                    'test_type': "--testOLS" if 'blank' in kfile or 'identity' in kfile else "--testOneRE"}
                   for (kfile, fn) in itertools.product(kfile_list, pheno_files)]

    matched_pheno_k_files = [(join(input_pheno_folder, x),
                              join(two_RE_kfile_folder, x) + '_2re.kin')
                             for x in listdir(input_pheno_folder)
                             if x + '_2re.kin' in listdir(two_RE_kfile_folder)]

    two_RE_kwargs_list = [{'plinkfile': input_plinkfile,
                           'kfile': kfile,
                           'phenofile': fn,
                           'gxe': True,
                           'noCorrect': True,
                           'output_folder': results_folder,
                           'covfile': input_covfile,
                           'test_type': ''}
                          for (fn, kfile) in matched_pheno_k_files]
    # pylmm_wrapper(kwargs_list[0])
    # pylmm_wrapper(two_RE_kwargs_list[0])
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1)
    sorted_execution_order = sorted(two_RE_kwargs_list + kwargs_list, key=lambda x: x['phenofile'])
    pool.map(pylmm_run_wrapper, sorted_execution_order)
    # for u in sorted_execution_order[:12]:
    #     print u
    # import time
    # time.sleep(60)


if __name__ == "__main__":
    project_folder = '../eqtl'
    study = 'HAEC'
    main(input_plink_folder='../eqtl/clean_plink',
         input_pheno_folder='../eqtl/separate_pheno',
         results_folder='../eqtl/pylmm_output',
         two_RE_kfile_folder='../eqtl/2RE_Ks',
         study_name=study)
