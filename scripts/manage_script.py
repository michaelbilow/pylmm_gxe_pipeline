from os.path import join, split, splitext, exists
from os import listdir, makedirs
import sys
import run_plink
import quantile_normalize
import gcta_prepare
import run_gcta
import grm_to_kin
import manage_kinship
import make_2RE_kinship
import run_pylmm
import qsub_preprocess
from shutil import rmtree

PROGRAMS_FOLDER = '../programs'
RAW_DATA_FOLDER_NAME = 'raw_data'
CLEAN_PLINK_FOLDER_NAME = 'clean_plink'
CLEAN_PHENO_FOLDER_NAME = 'clean_pheno'
SEPARATE_PHENO_FOLDER_NAME = 'separate_pheno'
TWO_RE_FOLDER_NAME = '2RE_Ks'
HSQ_FOLDER_NAME = 'hsq'
PYLMM_OUTPUT_FOLDER_NAME = 'pylmm_output'


def main(project_folder, study_name, max_phenos=10, pipeline_step_start=0, run_local=False, limit=True):
    raw_data_folder = join(project_folder, RAW_DATA_FOLDER_NAME)
    clean_plink_folder = join(project_folder, CLEAN_PLINK_FOLDER_NAME)
    clean_pheno_folder = join(project_folder, CLEAN_PHENO_FOLDER_NAME)
    separate_pheno_folder = join(project_folder, SEPARATE_PHENO_FOLDER_NAME)
    hsq_folder = join(project_folder, HSQ_FOLDER_NAME)
    two_re_folder = join(project_folder, TWO_RE_FOLDER_NAME)
    pylmm_output_folder = join(project_folder, PYLMM_OUTPUT_FOLDER_NAME)

    if pipeline_step_start <= 1:
        if not exists(raw_data_folder):
            raise IOError("No Raw Data Folder")

        if not exists(clean_plink_folder):
            makedirs(clean_plink_folder)
        else:
            sys.stderr.write("Clean PLINK folder exists; overwriting\n")
            rmtree(clean_plink_folder)
            makedirs(clean_plink_folder)

        print 'Running PLINK...'
        run_plink.main(input_folder=raw_data_folder, output_folder=clean_plink_folder, study_name=study_name)

    if pipeline_step_start <= 2:
        if not exists(clean_pheno_folder):
            makedirs(clean_pheno_folder)
        else:
            sys.stderr.write("Clean Pheno folder exists; overwriting\n")
            rmtree(clean_pheno_folder)
            makedirs(clean_pheno_folder)
        print 'Quantile Normalizing...'
        quantile_normalize.main(input_folder=raw_data_folder,
                                study_name=study_name,
                                plink_folder=clean_plink_folder,
                                output_pheno_folder=clean_pheno_folder,
                                max_phenos=max_phenos)

    if pipeline_step_start <= 3:
        if not exists(separate_pheno_folder):
            makedirs(separate_pheno_folder)
        else:
            sys.stderr.write("Separate Pheno folder exists; overwriting\n")
            rmtree(separate_pheno_folder)
            makedirs(separate_pheno_folder)
        print 'Preparing for GCTA'
        gcta_prepare.main(pheno_folder=clean_pheno_folder,
                          pheno_output_folder=separate_pheno_folder,
                          max_phenos=max_phenos)

    if pipeline_step_start <= 4:
        if not exists(hsq_folder):
            makedirs(hsq_folder)
        else:
            sys.stderr.write("HSQ folder exists; overwriting\n")
            rmtree(hsq_folder)
            makedirs(hsq_folder)
        run_gcta.main(plink_folder=clean_plink_folder,
                      pheno_folder=separate_pheno_folder,
                      hsq_folder=hsq_folder,
                      study_name=study_name)

    if pipeline_step_start <= 5:
        grm_to_kin.main(input_folder=clean_plink_folder, study_name=study_name)

    if pipeline_step_start <= 6:
        manage_kinship.main(kinship_fn=[join(clean_plink_folder, x)
                                       for x in listdir(clean_plink_folder)
                                       if x.endswith('.kin')][0])
        # Blank is interpreted as the all -1 matrix;

    if pipeline_step_start <= 7:
        if not exists(two_re_folder):
            makedirs(two_re_folder)
        else:
            sys.stderr.write("2RE-Kinship folder exists; overwriting\n")
            rmtree(two_re_folder)
            makedirs(two_re_folder)
        make_2RE_kinship.main(kinship_matrix_folder=clean_plink_folder,
                              hsq_folder=hsq_folder,
                              output_2RE_folder=two_re_folder,
                              study_name=study_name)

    if run_local:
        if not exists(pylmm_output_folder):
            makedirs(pylmm_output_folder)
        else:
            sys.stderr.write("pylmm output folder exists; overwriting;")
            rmtree(pylmm_output_folder)
            makedirs(pylmm_output_folder)
        run_pylmm.main(input_plink_folder=clean_plink_folder,
                       input_pheno_folder=separate_pheno_folder,
                       results_folder=pylmm_output_folder,
                       two_RE_kfile_folder=two_re_folder,
                       study_name=study_name,
                       limit=limit)
        return
    else:
        qsub_preprocess.main(plink_folder=clean_plink_folder,
                             separate_pheno_folder=separate_pheno_folder,
                             two_RE_kin_folder=two_re_folder)


if __name__ == "__main__":
    proj_folder = '../eqtl'
    study = 'HAEC'
    main(project_folder=proj_folder, study_name=study, pipeline_step_start=0,
         max_phenos=100, limit=True, run_local=True)
