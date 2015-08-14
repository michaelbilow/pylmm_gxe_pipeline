import subprocess
import os
import numpy as np
PYLMM_SCRIPT_LOC = '~/anaconda/bin/pylmmKinship.py'


def kinship_run(plink):
    pylmm_command = 'python {script_loc} -v --bfile {plink} {outfile}'

    out_fn = '{}.pylmm.kin'.format(os.path.splitext(plink)[0])
    this_command = pylmm_command.format(plink=plink, script_loc=PYLMM_SCRIPT_LOC,
                                        outfile=out_fn)
    subprocess.call(this_command, shell=True)
    return out_fn


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


if __name__ == "__main__":
    input_folder = '../clean_plink/'
    input_plink_fn = os.path.join(input_folder, 'hmdp')
    kinship_fn = kinship_run(input_plink_fn)

    kinship_matrix = np.loadtxt(kinship_fn)

    gxe_filename = os.path.join(input_folder, 'hmdp.gxe')
    gxe_kinship_matrix = populate_gxe_kinship(kinship_matrix, gxe_filename)
    np.savetxt('{}_gxe{}'.format(*os.path.splitext(kinship_fn), delimiter=' '), gxe_kinship_matrix)