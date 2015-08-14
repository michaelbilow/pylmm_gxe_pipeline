import pandas as pd
from os.path import join, split, splitext, exists
from os import listdir, makedirs
import shutil

BAD_COLS = ['FID', 'IID', 'Sex', 'covariate_thioglycolate_treated']


def make_trivial_gxe(cov_fn):
    shutil.copy(src=cov_fn, dest='{}.gxe'.format(splitext(cov_fn)[0]))
    return


def make_gxe(cov_fn, fam_fn):
    fam_df = pd.read_csv(fam_fn, sep=' ', header=None)
    fam_df = fam_df[fam_df.columns[:2]]
    fam_df.columns = ['FID', 'IID']
    cov_df = pd.read_csv(cov_fn, sep='\t')
    gxe_df = pd.merge(fam_df, cov_df, on='IID')
    gxe_fn = '{}.gxe'.format(splitext(cov_fn)[0])
    gxe_df.to_csv(gxe_fn, sep='\t', index=False, header=False)
    return


def preprocess_hmdp(hmdp_fn, hmdp_output_fn):
    df = pd.read_csv(hmdp_fn, sep='\t')
    iids = df['Strain'] + '__' + df['Animal_Id'].astype(str)
    df['IID'] = iids
    df.drop(['Strain', 'Animal_Id'], axis=1, inplace=True)
    df.to_csv(hmdp_output_fn, index=False, sep='\t')
    return


def make_phenos(pheno_fn, max_phenos):
    pheno_df = pd.read_csv(pheno_fn, sep='\t')
    pheno_df.fillna('NA', inplace=True)
    # print pheno_df.head()
    count = 0
    pheno_fn_start = join(split(pheno_fn)[0], split(pheno_fn)[1].split('_')[0])
    print pheno_fn_start
    ext = '_'.join(split(pheno_fn)[1].split('_')[1:])
    print ext
    for col in pheno_df.columns:
        count += 1
        if count > max_phenos:
            break
        if col in BAD_COLS:
            continue
        output_fn = '{}___{}_{}'.format(pheno_fn_start, col, ext)
        output_df = pheno_df[['FID', 'IID', col]]
        output_df.to_csv(output_fn, sep='\t', index=False, header=False)
    return


if __name__ == "__main__":
    # input_pheno_folder = '../raw_pheno'
    # pheno_fns = [os.path.join(input_pheno_folder, x) for x in os.listdir(input_pheno_folder)
    #              if os.path.splitext(x)[-1] == '.pheno']
    # output_plink_folder = '../clean_plink/'
    # filtered_expanded_fam_fn = os.path.join(output_plink_folder, 'hmdp.fam')
    # for fn in pheno_fns:
    #     make_phenos(fn, filtered_expanded_fam_fn)
    # covariate_fn = '../raw_pheno/hmdp.cov'
    # make_gxe(covariate_fn, filtered_expanded_fam_fn)
    #
    # output_pheno_folder = '../clean_pheno'
    # for x in os.listdir(input_pheno_folder):
    #     if '___' in x:
    #         shutil.move(os.path.join(input_pheno_folder, x), os.path.join(output_pheno_folder, x))
    #
    # if 'hmdp.gxe' in os.listdir(input_pheno_folder):
    #     shutil.move(os.path.join(input_pheno_folder, 'hmdp.gxe'),
    #                 os.path.join(output_plink_folder, 'hmdp.gxe'))

    input_pheno_folder = '../eqtl/raw_pheno'
    output_pheno_folder = '../eqtl/clean_pheno'

    if not exists(output_pheno_folder):
        makedirs(output_pheno_folder)
    max_phenos = 100

    pheno_fns = [join(input_pheno_folder, x) for x in listdir(input_pheno_folder)
                 if splitext(x)[-1].startswith('.pheno')]
    output_plink_folder = '../eqtl/clean_plink/'

    for fn in pheno_fns:
        make_phenos(fn, max_phenos)

    for x in listdir(input_pheno_folder):
        if '___' in x:
            shutil.move(join(input_pheno_folder, x), join(output_pheno_folder, x))
