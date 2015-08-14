from scipy.stats import norm
from scipy.stats import rankdata
import numpy as np
import os
import pandas as pd
from collections import OrderedDict
from os.path import join, split, splitext

pd.options.display.float_format = '{:.3f}'.format


def quantile_normalize(arr, cov):
    """
    Quanitle-normalizes an array with nan values
    :param arr:
    :return:
    """
    assert -9 not in arr
    # -9 is sometimes used as a nan value, but it could be a number. Better to raise
    # this error than allow something unexpected to happen.

    if len(cov):
        cov = cov.flatten()
        unique_cov_vals = np.unique(cov)
        assert len(unique_cov_vals) == 2
        output_array = np.zeros(arr.shape)
        for cov_val in unique_cov_vals:
            cov_val_mask = cov == cov_val
            pheno_vals = arr[cov_val_mask]
            qnormed_vals = qnorm(pheno_vals)
            output_array[cov_val_mask] = qnormed_vals
        return output_array
    else:
        return qnorm(arr)


def qnorm(arr):
    argsort = rankdata(arr, method='average')
    argsort[np.isnan(arr)] = np.nan
    argsort -= .5
    percentage_ranks = argsort / ((1 - np.isnan(arr)).sum())
    try:
        assert np.isclose(np.nanmean(percentage_ranks), .5)  # The mean of the ranks should be .5
    except AssertionError:
        print percentage_ranks
    # Because scipy.stats.norm does not raise an error when the quantile function (norm.ppf)
    # is given a number outside of the open interval (0, 1), we should check.
    assert np.nanmin(percentage_ranks) > 0
    assert np.nanmax(percentage_ranks) < 1
    qnorms = np.array([norm.ppf(x) for x in percentage_ranks])
    qnorms[qnorms == -9] += .001
    # In the highly unlikely case that -9 actually shows up as a qnorm value,
    # we jiggle the value.
    # (The reason, as above, is that -9 is used to code N/A phenotypes in some PLINK files)
    return qnorms


def strip_phenos(fam_fn):
    with open(fam_fn, 'r') as input_f:
        output_fn = os.path.splitext(fam_fn)[0] + '.pheno'
        with open(output_fn, 'w') as pheno_f:
            for line in input_f:
                line = line.strip().split()
                output_line = [line[0], line[1], line[-1]]
                pheno_f.write('{}\n'.format('\t'.join(output_line)))

    return


def qnorm_pheno_file(pheno_fn):
    real_pheno_indivs = []
    real_pheno_vals = []
    output_phenos = []
    with open(pheno_fn, 'r') as pheno_f:
        for line in pheno_f:
            line = line.strip().split()
            if line[-1] in ("-9", "NA", "N/A"):
                output_phenos.append(line)
            else:
                output_phenos.append([])
                real_pheno_indivs.append(line[:-1])
                real_pheno_vals.append(float(line[-1]))

    qnormed_phenos = quantile_normalize(real_pheno_vals)
    real_phenos = [real_pheno_indivs[i] + [qnormed_phenos[i]] for i in range(len(qnormed_phenos))]

    output_ind = 0
    ind = 0
    while output_ind < len(output_phenos):
        if not output_phenos[output_ind]:
            output_phenos[output_ind] = real_phenos[ind]
            ind += 1
        else:
            pass
        output_ind += 1

    output_fn = '{}_qnormed{}'.format(*os.path.splitext(pheno_fn))
    with open(output_fn, 'w') as output_f:
        for row in output_phenos:
            output_str = '{}\t{}\t{}\n'.format(*row)
            output_f.write(output_str)


def make_all_qnormed_phenos(df, exclude_cols, cov_col, start_dict, output_fn):
    cov_values = cov_col if not cov_col else df[cov_col].values
    qnormed_cols = []
    for col in df.columns:
        if col not in exclude_cols:
            # print col, len(cov_values)
            qnormed_vals = quantile_normalize(df[col].values, cov_values)
            qnormed_cols.append((col, qnormed_vals))
    qnormed_dict = OrderedDict(qnormed_cols)
    output_dict = OrderedDict(start_dict.items() + qnormed_dict.items())
    output_df = pd.DataFrame(output_dict)
    output_df.fillna('NA', inplace=True)
    print output_df.head()
    output_df.to_csv(output_fn, index=False, sep='\t')

    return output_df


def qnorm_multiple(input_fn, keep_cols, drop_cols=(), cov_col=()):
    assert len(cov_col) == 1
    df = pd.read_csv(input_fn, sep='\t')
    print df.head()
    start_dict = OrderedDict([(col, df[col]) for col in keep_cols])
    output_fn = '{}_qnormed_between_group'.format(*os.path.splitext(input_fn))
    exclude_cols = list(keep_cols) + list(drop_cols) + list(cov_col)
    make_all_qnormed_phenos(df=df, exclude_cols=exclude_cols, cov_col=(),
                            start_dict=start_dict, output_fn=output_fn)

    if cov_col:
        cov_dict = OrderedDict([(col, df[col]) for col in cov_col])
        cov_output_dict = OrderedDict(start_dict.items() + cov_dict.items())
        cov_df = pd.DataFrame(cov_output_dict)
        cov_fn = '{}.cov'.format(os.path.splitext(input_fn)[0])
        cov_df.to_csv(cov_fn, index=False, sep='\t')
        cov_output_fn = '{}_qnormed_within_group'.format(input_fn)
        make_all_qnormed_phenos(df=df, exclude_cols=exclude_cols,
                                cov_col=cov_col,
                                start_dict=start_dict,
                                output_fn=cov_output_fn)

    return


def haec_preprocess(pheno_fn, cov_fn, bad_rows_fn, fam_fn, output_pheno_fn):
    df = pd.read_csv(pheno_fn, sep='\t', header=None)
    print df.head()
    print df.values.shape
    bad_rows = set([int(x) - 1 for x in open(bad_rows_fn, 'r')])
    cleaned_df = df.iloc[[i for i in range(len(df)) if i not in bad_rows]]
    # print bad_rows
    print cleaned_df.head()
    print cleaned_df.values.shape, len(bad_rows)
    output_df = cleaned_df.T
    output_df.columns = ['E_{:05}'.format(int(x) + 1) for x in output_df.columns]

    fam_df = pd.read_csv(fam_fn, sep=' ', header=None)
    fam_df.columns = ['FID', 'IID', 'FatherID', 'MotherID', 'Sex', 'Pheno']
    fam_df = fam_df[['FID', 'IID', 'Sex']]
    print fam_df.values.shape
    output_df = pd.merge(output_df, fam_df, left_index=True, right_index=True)

    output_df['Env'] = [int(x) for x in open(cov_fn, 'r')]
    output_df.to_csv(output_pheno_fn, sep='\t', index=False)
    print output_df.head()
    return


if __name__ == "__main__":
    # z = np.repeat(np.arange(20), 2)
    # z = np.concatenate((np.arange(10), [np.nan] * 10, np.arange(10), np.arange(30)))
    # # z = np.arange(40)
    # output = quantile_normalize(z)
    # print output

    # strip_phenos('./mousedata/hmdp1.fam')
    # qnorm_pheno_file('./mousedata/hmdp1.pheno')
    # print rankdata(z)
    # print np.mean(z[~np.isnan(z)])
    # print z[~np.isnan(z)]
    # print quantile_normalize(z)

    # input_fn = '../hmdp/raw_pheno/hmdp.pheno'
    # keep_these = ["Strain", "Animal_Id"]
    # drop_these = ["Sex", "heart_wt"]
    # covariates = ["covariate_thioglycolate_treated"]
    # qnorm_multiple(input_fn, keep_these, drop_these, covariates)

    eqtl_folder = '../eqtl/raw_pheno'
    input_fn = join(eqtl_folder, 'expr.merged.R.txt')
    exclude_rows = join(eqtl_folder, 'normal_test_0_05_sd_test_bad_expr.index.txt')
    covariate_fn = join(eqtl_folder, 'env_0_1.txt')
    output_fn = join(eqtl_folder, 'HAEC.pheno')

    plink_folder = '../eqtl/raw_plink'

    haec_preprocess(input_fn, covariate_fn, exclude_rows,
                    fam_fn=join(plink_folder, 'HAEC.fam'), output_pheno_fn=output_fn)

    qnorm_multiple(output_fn, keep_cols=['FID', 'IID'], drop_cols=['Sex'], cov_col=['Env'])
