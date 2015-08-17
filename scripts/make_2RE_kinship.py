import numpy as np
import pandas as pd
from os.path import join, splitext, split, exists
from os import listdir, makedirs, remove


def get_weights(hsq_fn):
    hsq_df = pd.read_csv(hsq_fn, sep='\t', index_col=0)
    # print hsq_df
    g_weight = hsq_df.loc['V(G)', 'Variance']
    gxe_weight = hsq_df.loc['V(GxE)', 'Variance']
    return g_weight, gxe_weight


def generate_2RE_kinship(hsq_fn, genetic_kinship, gxe_kinship):
    g_weight, gxe_weight = get_weights(hsq_fn)
    kinship_2re = genetic_kinship*g_weight + gxe_kinship*gxe_weight
    return kinship_2re


def main(kinship_matrix_folder, hsq_folder, output_2RE_folder, study_name):
    genetic_kinship_fn = join(kinship_matrix_folder, '{}.grm.kin'.format(study_name))
    gxe_kinship_fn = join(kinship_matrix_folder, '{}.grm_gxe.kin'.format(study_name))

    genetic_kinship_matrix = np.loadtxt(genetic_kinship_fn)
    gxe_kinship_matrix = np.loadtxt(gxe_kinship_fn)

    hsq_fns = [join(hsq_folder, x) for x in listdir(hsq_folder)]
    if not exists(output_2RE_folder):
        makedirs(output_2RE_folder)
    else:  # Clear it out
        for x in listdir(output_2RE_folder):
            remove(join(output_2RE_folder, x))
    for fn in hsq_fns:
        print fn
        twoRE_kinship = generate_2RE_kinship(fn, genetic_kinship_matrix, gxe_kinship_matrix)
        output_fn = join(output_2RE_folder, '{}_2re.kin'.format(splitext(split(fn)[1])[0]))
        np.savetxt(output_fn, twoRE_kinship, delimiter=' ')


if __name__ == "__main__":
    project_folder = '../eqtl'
    main(kinship_matrix_folder='../eqtl/clean_plink',
         hsq_folder='../eqtl/hsq',
         output_2RE_folder='../eqtl/2RE_Ks',
         study_name='HAEC')
