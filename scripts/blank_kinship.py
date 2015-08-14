import pandas as pd
import os
import numpy as np


def blank_kinship(input_fn):
    df = pd.read_csv(input_fn, sep=' ', header=None)
    df.loc[:, :] = 'NA'
    df.to_csv('{}_blank{}'.format(*os.path.splitext(input_fn)), index=False, sep=' ', header=False)


def identity_kinship(input_fn):
    df = pd.read_csv(input_fn, sep=' ', header=None)
    output_df = pd.DataFrame(np.eye(len(df)))
    output_df.to_csv('{}_identity{}'.format(*os.path.splitext(input_fn)), index=False, sep=' ', header=False)


def main(kinship_fn):
    blank_kinship(kinship_fn)
    # identity_kinship(kinship_fn)
    return

if __name__ == "__main__":
    input_kinship_fn = '../hmdp/clean_plink/hmdp.grm.kin'

    input_kinship_fn = '../eqtl/clean_plink/HAEC.grm.kin'
    main(input_kinship_fn)

