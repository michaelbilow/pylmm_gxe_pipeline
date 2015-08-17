from os.path import join, split, splitext, exists
from os import listdir, rename


def main(plink_folder, two_RE_kin_folder, separate_pheno_folder):
    output_fn = join(plink_folder, 'pheno_aliases.txt')

    pheno_fns = listdir(separate_pheno_folder)
    two_re_fns = listdir(two_RE_kin_folder)

    assert len(pheno_fns) == len(two_re_fns)
    for x in zip(pheno_fns, two_re_fns):
        print x

    assert all(two_re_fns[i].startswith(pheno_fns[i]) for i in range(len(pheno_fns)))

    with open(output_fn, 'w') as f:
        for i in range(len(pheno_fns)):
            pheno_fn = pheno_fns[i]
            k_fn = two_re_fns[i]
            new_pheno_fn = 'pheno.{}'.format(i)
            new_k_fn = k_fn.replace(pheno_fn, new_pheno_fn)
            line = '{}\t{}\n'.format(i, pheno_fn)
            f.write(line)
            rename(join(separate_pheno_folder, pheno_fn),
                   join(separate_pheno_folder, new_pheno_fn))
            rename(join(two_RE_kin_folder, k_fn),
                   join(two_RE_kin_folder, new_k_fn))


if __name__ == "__main__":
    main('../eqtl/clean_plink', '../eqtl/2RE_Ks', '../eqtl/separate_pheno')
