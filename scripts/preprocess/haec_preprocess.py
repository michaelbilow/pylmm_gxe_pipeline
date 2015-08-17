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