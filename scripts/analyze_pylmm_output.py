import pandas as pd
import matplotlib.pyplot as plt
from os.path import join, split, splitext, exists
from os import listdir, makedirs
import numpy as np


def translate_title(title):
    K_type = title.split('_')[-1]
    phen = ' '.join(title.split('_')[:-1])
    if any(u in K_type for u in ('GxE',)):
        pretty_K_type = K_type
    elif 'BLANK' in K_type.upper():
        pretty_K_type = 'OLS'
    elif '2RE' in K_type.upper():
        pretty_K_type = 'Two RE'
    else:
        pretty_K_type = 'One RE'
    return phen, pretty_K_type


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


def main(input_folder):
    plt.rc('text', usetex=True)
    input_fns = [join(input_folder, x) for x in listdir(input_folder) if 'stats_GxE' in x]
    print input_fns
    output_folder = join(split(input_folder)[0], 'pylmm_stats')
    if not exists(output_folder):
        makedirs(output_folder)
    output_records = []
    for fn in input_fns:
        output_lambda = read_input(fn, output_folder)
        title = split(fn)[-1].split('___')[1].split('.')[0]
        pheno_name, K_name = translate_title(title)
        output_record = {'lambda': output_lambda,
                         'QNormed': 'QNormed' if 'qnormed' in fn else 'Not QNormed',
                         'Phenotype': pheno_name,
                         'K_type': K_name}
        output_records.append(output_record)
    import pandas as pd
    df = pd.DataFrame(output_records)
    df.to_csv(join(output_folder, 'pylmm_inflation_raw.csv'))
    df.set_index(['Phenotype', 'QNormed', 'K_type'], inplace=True)
    df = df.unstack((-1, -2))
    df.to_csv(join(output_folder, 'pylmm_inflation.csv'))
    plt.figure()
    df.columns = ['\n'.join(x) for x in df.columns]
    df.boxplot()
    plt.savefig(join(output_folder, 'boxplot.png'))
    return


if __name__ == "__main__":
    input_folder = '../eqtl/pylmm_output'
    main(input_folder)


