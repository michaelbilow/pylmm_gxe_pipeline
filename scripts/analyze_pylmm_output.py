import pandas as pd
import matplotlib.pyplot as plt
from os.path import join, split, splitext, exists
from os import listdir, makedirs
import numpy as np


def translate_title(title):
    K_type = title.split('_')[-1].split('.')[0]
    pretty_K_type = '2RE' if K_type == 'K-2RE' \
        else '1RE' if K_type == 'K-G' \
        else 'OLS' if K_type == 'K-Id' \
        else 'K-Unknown'
    if 'original' in title:
        phen = title.split('_original_')[0]
        qnorm = 'Original*Phenotype'
    else:
        phen = title.split('_qnormed_')[0]
        raw_qnorm = ' '.join(title.split('_qnormed_')[1].split('_')[:-1])
        good_qnorm = ' '.join(x.capitalize() for x in raw_qnorm.split(' '))
        qnorm = 'Quanitle Normalize*{}'.format(good_qnorm)
    return phen, qnorm, pretty_K_type


def read_input(fn, output_folder, graph=False):
    try:
        df = pd.read_csv(fn, sep='\t')
    except ValueError:
        return np.nan
    # print df.head()
    p_vals = df.P_VALUE.values
    p_vals = p_vals[~np.isnan(p_vals)]
    p_vals.sort()
    sorted_p_vals = p_vals
    sorted_logp = -np.log10(sorted_p_vals)
    try:
        lambda_val = sorted_logp[-len(sorted_logp)/2]/(-np.log10(.5))
    except IndexError:
        return np.nan
    if graph:
        hypothesized_p_vals = np.linspace(0, 1, len(sorted_p_vals) + 1)[1:]
        hypothesized_logp = -np.log10(hypothesized_p_vals)

        plt.figure()
        df = pd.DataFrame({"-LogP": sorted_logp,
                          "Expected -LogP": hypothesized_logp})
        df.plot(kind='scatter', y='-LogP', x="Expected -LogP")

        plt.plot(hypothesized_logp, hypothesized_logp, 'black')
        title = split(fn)[-1].split('___')[1].split('.')[0]
        print title
        pheno_name, K_name = translate_title(title)
        good_title = 'Phenotype: {}\t K: {}'.format(pheno_name, K_name)
        plt.title(good_title)
        plt.plot(hypothesized_logp, lambda_val*hypothesized_logp, 'r')
        middle_point = (max(hypothesized_logp)/2.0, lambda_val * max(hypothesized_logp)/2.0)
        text_point = (middle_point[0] - 1, middle_point[1] + 1)
        plt.xlabel('Predicted $-\log_{10}(p)$')
        plt.ylabel('Experimental $-\log_{10}(p)$')
        plt.annotate(xy=middle_point, xytext=text_point, s="$\lambda=$" + '{:.3}'.format(lambda_val))
        plt.savefig('{}_qq_plot.png'.format(join(output_folder, split(splitext(fn)[0]))[1]))
    return lambda_val


def main(input_folder, reuse=False):
    output_folder = join(split(input_folder)[0], 'pylmm_stats')

    if not reuse:
        input_fns = [join(input_folder, x) for x in listdir(input_folder) if 'stats_GxE' in x]
        print input_fns
        if not exists(output_folder):
            makedirs(output_folder)
        output_records = []
        for fn in input_fns:
            output_lambda = read_input(fn, output_folder)
            title = split(fn)[-1].split('___')[1].split('.')[0]
            pheno_name, qnorm, K_name = translate_title(title)
            output_record = {'lambda': output_lambda,
                             'QNormed': qnorm,
                             'Phenotype': pheno_name,
                             'K_type': K_name}
            if 'Unknown' in K_name:
                continue
            output_records.append(output_record)
        df = pd.DataFrame(output_records)
        df.to_csv(join(output_folder, 'pylmm_inflation_raw.csv'), index=True)
        df.set_index(['Phenotype', 'QNormed', 'K_type'], inplace=True)
        df = df.unstack((-2, -1))
        # df.sortlevel(level=2, axis=1, inplace=True)
        col_tuples = [x[1:] for x in df.columns]
        df.columns = pd.MultiIndex.from_tuples(col_tuples)
        df.to_csv(join(output_folder, 'pylmm_inflation.csv'), index=True)
    else:
        df = pd.read_csv(join(output_folder, 'pylmm_inflation.csv'), header=[0, 1],
                         skipinitialspace=True, tupleize_cols=True)
        df.columns = pd.MultiIndex.from_tuples(df.columns)


    plt.rc('text', usetex=True)
    df.columns = ['\n'.join(x).replace('*', '\n') for x in df.columns]
    for i in range(3):
        plt.figure(figsize=(10, 10))
        tiny_df = df[df.columns[3*i: 3*i+3]]
        print tiny_df.columns
        correct_col_order = sorted(tiny_df.columns, key=lambda col: 0 if col.endswith('OLS')
                                                        else 1 if col.endswith('1RE') else 2)
        tiny_df = tiny_df[correct_col_order]
        norm_type = ' '.join(tiny_df.columns[0].split('\n')[:2])
        tiny_df.boxplot()
        plt.title('Distribution of Inflation Factors for HAEC eQTL Studies ({})'.format(norm_type))
        plt.savefig(join(output_folder, 'boxplot_{}.png'.format(''.join((x[0] for x in norm_type.split())))))
    return


if __name__ == "__main__":
    input_folder = '../eqtl/pylmm_output (2)'
    main(input_folder, reuse=False)


