import pandas as pd
import numpy as np
from os.path import join, isfile, splitext
from os import listdir
import itertools
import matplotlib.pyplot as plt
from scipy import stats



def compare_pheno_files(fn1, fn2, output_folder, draw=True):
    df1 = pd.read_csv(fn1, sep='\t', header=None)
    df2 = pd.read_csv(fn2, sep='\t', header=None)

    is_qnormed = lambda x: 'qnormed' in x
    qnormed1, qnormed2 = is_qnormed(fn1), is_qnormed(fn2)
    if all((qnormed1, qnormed2)):
        return False

    type_of_pheno = lambda x: 'QNormed Within Group' if 'within' in x \
        else 'QNormed Between Groups' if 'between' in x else 'Original'
    short_type_of_pheno = lambda x: 'within' if 'Within' in x else 'between' if 'Between' in x else 'original'
    pheno_name = fn1.split('___')[1]

    type1, type2 = type_of_pheno(fn1), type_of_pheno(fn2)
    df1.columns = ['FID', 'IID', type1]
    df2.columns = ['FID', 'IID', type2]

    if qnormed1:
        qnormed_col, not_qnormed_col = type1, type2
    else:
        qnormed_col, not_qnormed_col = type2, type1

    combined_df = pd.merge(df1, df2, on=['FID', 'IID'])
    combined_df.dropna(how='any', inplace=True)
    slope, intercept, r_value, p_value, std_err = stats.linregress(combined_df[qnormed_col],
                                                                   combined_df[not_qnormed_col])
    print slope, intercept
    if draw:
        fig = plt.figure()
        combined_df.plot(kind='scatter', y=not_qnormed_col, x=qnormed_col)
        plt.title('Phenotype: {} -- {} vs. {}'.format(pheno_name.replace('_', ' '),
                                                      qnormed_col, not_qnormed_col))
        best_fit_values = slope*combined_df[qnormed_col] + intercept
        plt.plot(combined_df[qnormed_col], best_fit_values, color='r')
        xlims, ylims = plt.xlim(), plt.ylim()
        bottom_right = (.25*xlims[0] + .75*xlims[1], .75*ylims[0] + .25*ylims[1])
        plt.text(*bottom_right, s=r'$\rho^2$ = {}'.format(r_value**2),
                 horizontalalignment='center')
        plt.xlabel('{} Value'.format(qnormed_col))
        plt.ylabel('{} Value'.format(not_qnormed_col))
        plt.grid()
        plt.savefig(join(output_folder,
                                 '{}_qnorm_{}_compare.png'.format(pheno_name,
                                                                  short_type_of_pheno(qnormed_col))))
        plt.close("all")
    return pheno_name, short_type_of_pheno(qnormed_col), r_value**2


if __name__ == "__main__":
    plt.rc('text', usetex=True)
    input_folder = '../clean_pheno'
    output_folder = join(input_folder, 'compare_qnormed')
    paired_pheno_fns = [(x, y) for x, y in
                        itertools.product(listdir(input_folder), listdir(input_folder))
                        if x < y and isfile(join(input_folder, x)) and isfile(join(input_folder, y))
                        and x.split('___')[1] == y.split('___')[1]]

    good_paired_pheno_fns = [(join(input_folder, x), join(input_folder, y))
                             for (x, y) in paired_pheno_fns]
    output_records = []
    for good_pair in good_paired_pheno_fns:
        this_output = compare_pheno_files(*good_pair, output_folder=output_folder)
        if not this_output:
            continue
        print good_pair
        pheno, qnorm_type, r2 = this_output
        output_dict = {'Phenotype': pheno, 'R^2': r2, 'Normalization': qnorm_type}
        output_records.append(output_dict)
    output_df = pd.DataFrame(output_records)
    output_df = output_df.pivot(index='Phenotype', columns='Normalization', values='R^2')
    output_df.to_csv(join(output_folder, 'QNormed Phenotype Comparisons.txt'), sep='\t')