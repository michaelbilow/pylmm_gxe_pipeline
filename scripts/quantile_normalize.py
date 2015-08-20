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


def make_all_qnormed_phenos(df, keep_cols, cov_col, start_dict, output_fn, no_qnorm=False, max_phenos=None):
    cov_values = cov_col if not cov_col else df[cov_col].values
    exclude_cols = keep_cols + list(cov_col) + [x for x in df.columns if x.startswith('Cov')]
    qnormed_cols = []
    count = 0
    if no_qnorm:
        output_df = df[keep_cols + [x for x in df.columns if x not in exclude_cols][:max_phenos]]
    else:
        for col in df.columns:
            if col in exclude_cols:
                continue
            count += 1
            if max_phenos and count > max_phenos:
                break
            qnormed_vals = quantile_normalize(df[col].values, cov_values)
            qnormed_cols.append((col, qnormed_vals))
        qnormed_dict = OrderedDict(qnormed_cols)
        output_dict = OrderedDict(start_dict.items() + qnormed_dict.items())
        output_df = pd.DataFrame(output_dict)
    output_df.fillna('NA', inplace=True)
    print output_df.head()
    output_df.to_csv(output_fn, index=False, sep='\t')
    return output_df


def qnorm_multiple(input_fn, keep_cols, cov_col, output_folder, max_phenos):
    assert len(cov_col) == 1
    df = pd.read_csv(input_fn, sep='\t')
    # print df.head()
    start_dict = OrderedDict([(col, df[col]) for col in keep_cols])
    pheno_file_name = split(input_fn)[1]
    original_output_fn = join(output_folder,
                              '{}_original'.format(pheno_file_name))
    print 'Copying original phenotypes...'
    make_all_qnormed_phenos(df=df,
                            keep_cols=keep_cols,
                            cov_col=cov_col,
                            start_dict=start_dict,
                            output_fn=original_output_fn,
                            no_qnorm=True,
                            max_phenos=max_phenos)

    print 'Quantile Normalizing whole dataset'
    between_output_fn = join(output_folder,
                     '{}_qnormed_between_group'.format(pheno_file_name))
    make_all_qnormed_phenos(df=df, keep_cols=keep_cols, cov_col=(),
                            start_dict=start_dict, output_fn=between_output_fn,
                            max_phenos=max_phenos)

    if cov_col:
        print 'Quantile Normalizing within groups'
        within_output_fn = join(output_folder,
                             '{}_qnormed_within_group'.format(pheno_file_name))
        make_all_qnormed_phenos(df=df, keep_cols=keep_cols,
                                cov_col=cov_col,
                                start_dict=start_dict,
                                output_fn=within_output_fn,
                                max_phenos=max_phenos)

    return


def main(input_folder, study_name, plink_folder, output_pheno_folder, max_phenos=None):
    input_fn = join(input_folder, '{}.pheno'.format(study_name))
    keep_cols = ['FID', 'IID']
    env_col = ['Env']
    input_df = pd.read_csv(input_fn, sep='\t')
    all_cov_cols = keep_cols + [x for x in input_df.columns if x.startswith('Cov')] + env_col
    cov_df = input_df[all_cov_cols]
    cov_output_fn = join(plink_folder, '{}.cov'.format(study_name))
    cov_df.to_csv(cov_output_fn,
                  index=False, sep='\t', header=False)
    gxe_df = input_df[keep_cols + env_col]
    gxe_output_fn = join(plink_folder, '{}.gxe'.format(study_name))
    gxe_df.to_csv(gxe_output_fn,
                  index=False, sep='\t', header=False)

    qnorm_multiple(input_fn, keep_cols=keep_cols, cov_col=env_col, output_folder=output_pheno_folder,
                   max_phenos=max_phenos)
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

    # eqtl_folder = '../eqtl/raw_pheno'
    # input_fn = join(eqtl_folder, 'expr.merged.R.txt')
    # exclude_rows = join(eqtl_folder, 'normal_test_0_05_sd_test_bad_expr.index.txt')
    # covariate_fn = join(eqtl_folder, 'env_0_1.txt')
    # output_fn = join(eqtl_folder, 'HAEC.pheno')
    #
    # plink_folder = '../eqtl/raw_plink'
    #
    # haec_preprocess(input_fn, covariate_fn, exclude_rows,
    #                 fam_fn=join(plink_folder, 'HAEC.fam'), output_pheno_fn=output_fn)
    #
    # qnorm_multiple(output_fn, keep_cols=['FID', 'IID'], drop_cols=['Sex'], cov_col=['Env'])

    main('../eqtl/raw_data', 'HAEC', '../eqtl/clean_plink/', '../eqtl/clean_pheno', max_phenos=10)
