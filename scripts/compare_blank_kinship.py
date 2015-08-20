import pandas as pd
from os.path import join, split, splitext
from os import listdir
import numpy as np

if __name__ == "__main__":
    input_folder = '../eqtl/pylmm_output'
    fns = [join(input_folder, x) for x in listdir(input_folder) if 'Blank' in x or 'Id' in x]
    print fns
    paired_fns = [(x, y) for x in fns for y in fns if x < y and x.split('K')[0] == y.split('K')[0]]
    for x in paired_fns:
        print x

    for pair in paired_fns:
        blank_df, id_df = pd.read_csv(pair[0], sep='\t'), pd.read_csv(pair[1], sep='\t')
        # print blank_df.head()
        # print id_df.head()
        print np.allclose(blank_df['P_VALUE'], id_df['P_VALUE'])