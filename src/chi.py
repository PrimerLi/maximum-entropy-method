def chi(G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A, omega, LambdaInverse):
    import numpy as np

    Niom = len(G_real_rotated)
    Nomega = len(omega)

    vector = G_real_rotated - K_real_rotated.dot(A)
    result = np.transpose(vector).dot(LambdaInverse).dot(vector)

    vector = G_imag_rotated - K_imag_rotated.dot(A)
    result = result + np.transpose(vector).dot(LambdaInverse).dot(vector)
    return result
