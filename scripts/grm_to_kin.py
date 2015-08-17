import os
import numpy as np
import gzip


def get_n_indivs(grm_fn):
    with open(grm_fn, "rb") as f:
        first = f.readline()     # Read the first line.
        f.seek(-2, 2)            # Jump to the second last byte.
        while f.read(1) != "\n": # Until EOL is found...
            f.seek(-2, 1)        # ...jump back the read byte plus one more.
        last = f.readline()      # Read last line.
    n_indivs = int(last.split()[0])
    return n_indivs


def populate_kinship(n_indivs, grm_fn):
    output_kinship = np.zeros((n_indivs, n_indivs))
    with open(grm_fn, 'r') as f:
        for line in f:
            vals = line.strip().split()
            indiv1, indiv2 = int(vals[0]) - 1, int(vals[1]) - 1
            kinship_val = float(vals[-1])
            output_kinship[indiv1, indiv2] += kinship_val
            if indiv1 == indiv2:
                continue
            else:
                output_kinship[indiv2, indiv1] += kinship_val
    return output_kinship


def read_gzipped_grm(grm_gz_fn):
    output_fn = os.path.splitext(grm_gz_fn)[0]
    with open(output_fn, 'w') as output_f:
        with gzip.open(grm_gz_fn, 'rb') as input_f:
            file_content = input_f.read()
            output_f.write(file_content)
    return


def populate_gxe_kinship(kinship, gxe_fn):
    vals = []
    with open(gxe_fn, 'r') as gxe_f:
        for row in gxe_f:
            val = int(row.strip().split()[-1])
            vals.append(val)
    vals = np.array(vals)
    mask = np.array([vals == vals[i] for i in range(len(vals))])
    mask = mask.astype(int)
    gxe_kinship = kinship*mask  # This is an element-wise product
    return gxe_kinship


def main(input_folder, study_name):
    input_grm_gz = os.path.join(input_folder, '{}.grm.gz'.format(study_name))
    read_gzipped_grm(input_grm_gz)
    input_grm_fn = os.path.join(input_folder, '{}.grm'.format(study_name))
    study_size = get_n_indivs(input_grm_fn)
    kinship_matrix = populate_kinship(study_size, input_grm_fn)

    input_gxe_fn = os.path.join(input_folder, '{}.gxe'.format(study_name))
    gxe_kinship_matrix = populate_gxe_kinship(kinship_matrix, input_gxe_fn)

    np.savetxt(input_grm_fn + '.kin', kinship_matrix, delimiter=' ')
    np.savetxt(input_grm_fn + '_gxe.kin', gxe_kinship_matrix, delimiter=' ')


if __name__ == "__main__":
    input_folder = '../hmdp/clean_plink'
    study = 'hmdp'
    # main(input_folder, study)

    input_folder = '../eqtl/clean_plink'
    study = 'HAEC'
    main(input_folder, study)
