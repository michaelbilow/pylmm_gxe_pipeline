from plinkio import plinkfile
import os
import pandas as pd
from collections import defaultdict

def create_output_sample_list(sample_list, iids):
    sample_info_dict = {x.iid: (x.fid, x.father_iid, x.mother_iid, x.sex, x.affection, x.phenotype)
                        for x in sample_list}
    output_sample_list = []
    bad_strains = []
    good_strains = []
    for iid in iids:
        strain = iid.split('__')[0].replace('-', '')
        try:
            fid, father_iid, mother_iid, sex, affection, phenotype = sample_info_dict[strain]
            output_sample = plinkfile.Sample(fid, iid, father_iid, mother_iid, sex, affection, phenotype)
            output_sample_list.append(output_sample)
            good_strains.append(strain)
        except KeyError:
            bad_strains.append(strain)
    print 'Have Pheno Only:', sorted(set(bad_strains))
    print 'Have Pheno + Geno:', sorted(set(good_strains))
    print sorted(sample_info_dict.keys())
    return output_sample_list


def recode_for_multiple_indivs(plink_fn, pheno_fn, output_plink_fn):
    plink_file = plinkfile.open(plink_fn)
    if not plink_file.one_locus_per_row():
        print("This script requires that snps are rows and samples columns.")
        exit(1)
    sample_list = plink_file.get_samples()
    for x in sample_list:
        print x.iid
    pheno_df = pd.read_csv(pheno_fn, sep="\t")
    iids = pheno_df['Strain'] + '__' + pheno_df['Animal_Id'].astype(str)
    # print iids
    strain_count = pheno_df.groupby(['Strain']).count()['Animal_Id'].to_dict()
    # print strain_count
    output_sample_list = create_output_sample_list(sample_list, iids)
    # for x in output_sample_list:
    #     print x.iid
    out_plink = plinkfile.create(output_plink_fn, output_sample_list)
    locus_list = plink_file.get_loci()
    count = 0
    for locus, row in zip(locus_list, plink_file):
        count += 1
        if count % 1000 == 0:
            print 'At SNP {}'.format(count)
        strain_genotype_dict = {x[0].iid: x[1] for x in zip(sample_list, row)}
        sample_strains = [sample.iid.split('__')[0].replace('-', '') for sample in output_sample_list]
        output_row = [strain_genotype_dict[x] for x in sample_strains]
        out_plink.write_row(locus, output_row)
    return


def check_plinkfile(plink_fn):
    out_plink = plinkfile.open(plink_fn)
    samples = out_plink.get_samples()
    locuses = out_plink.get_loci()
    count = 0
    for locus, row in zip(locuses, out_plink):
        for sample, genotype in zip(samples, row):
            print("Individual {0} has genotype {1} for snp {2}.".format(sample.iid, genotype, locus.name))
        count += 1
        if count >= 2:
            break
    print len(locuses)

if __name__ == "__main__":
    plink_folder = '../raw_plink'
    pheno_folder = '../raw_pheno'
    output_fn = '../expanded/hmdp'
    recode_for_multiple_indivs(plink_folder + '/hmdp', pheno_folder + '/hmdp.pheno', output_fn)
    check_plinkfile(output_fn)
