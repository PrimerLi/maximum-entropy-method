def f(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A, D, omega, LambdaInverse, lambd):
    import numpy as np

    domega = omega[1] - omega[0]
    Nomega = len(omega)
    result = np.zeros(Nomega)
    for nw in range(Nomega):
        off_diagonal_terms = 0 
        if (nw + 2 < Nomega):
            off_diagonal_terms += 0.25*lambd/domega*A[nw+2]
        if (nw - 2 >= 0):
            off_diagonal_terms += 0.25*lambd/domega*A[nw-2]
        result[nw] = -alpha*domega*(1 + np.log(A[nw]/D[nw])) - lambd/domega*0.5*A[nw] + off_diagonal_terms

    vector = G_real_rotated - K_real_rotated.dot(A)
    result = result + np.transpose(K_real_rotated).dot(LambdaInverse).dot(vector)

    vector = G_imag_rotated - K_imag_rotated.dot(A)
    result = result + np.transpose(K_imag_rotated).dot(LambdaInverse).dot(vector)

    return result
