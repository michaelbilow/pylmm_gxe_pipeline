import numpy as np
import os

def compare_matrices(X, Y):
    """
    Checks to what extent the matrices X and Y are
    :param X:
    :param Y:
    :return:
    """
    assert X.shape == Y.shape  # Should have the same shape
    assert not np.allclose(X, np.zeros(X.shape)) and not np.allclose(Y, np.zeros(Y.shape))  # Should not be zero matrices

    X_zeros = np.isclose(X, np.zeros(X.shape))
    Y_zeros = np.isclose(Y, np.zeros(Y.shape))

    same_zeroes = np.logical_and(X_zeros, Y_zeros)
    same_not_zeroes = np.logical_and(~X_zeros, ~Y_zeros)

    assert np.ndarray.all(np.logical_or(same_zeroes, same_not_zeroes)) # Should have the same set of zeroes
    X[same_zeroes] = 10**-6
    Y[same_zeroes] = 10**-6

    output = X/Y
    assert np.allclose(output, output.T)  # Should be symmetrical
    median_factor = np.median(output)
    print 'Median difference in factor: {}'.format(median_factor)
    output /= median_factor
    return output

def find_differences(ratio_matrix):
    """
    Find anything out of +/- 5% from equal
    :param ratio_matrix:
    :return: list of indices where the two differ
    """

    mask = np.logical_or(ratio_matrix < .95, ratio_matrix > 1.05)
    indices = np.array([[(ind1 + 1, ind2 + 1) for ind2 in range(ratio_matrix.shape[1])]
                         for ind1 in range(ratio_matrix.shape[0])])
    output = indices[mask]
    return output


if __name__ == "__main__":
    input_folder = '../clean_plink'
    grm_fn = os.path.join(input_folder, 'hmdp.grm.kin')
    pylmm_kinship_fn = os.path.join(input_folder, 'hmdp.pylmm.kin')
    grm_matrix = np.loadtxt(grm_fn)
    pylmm_matrix = np.loadtxt(pylmm_kinship_fn)
    pylmm_grm_ratio_matrix = compare_matrices(grm_matrix, pylmm_matrix)
    pylmm_vs_grm_fn = os.path.join(input_folder, 'grm_vs_pylmm.txt')
    np.savetxt(pylmm_vs_grm_fn, pylmm_grm_ratio_matrix, fmt='%.4f')
    different_indices = find_differences(pylmm_grm_ratio_matrix)

    for x in different_indices:
        print x