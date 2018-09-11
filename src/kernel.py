import numpy as np

def K_real(omega_n, omega):
    return -omega/(omega_n**2 + omega**2)

def K_imag(omega_n, omega):
    return -omega_n/(omega_n**2 + omega**2)

def K_matrix_real(omega_n, omega):
    Niom = len(omega_n)
    Nomega = len(omega)
    matrix = np.zeros((Niom, Nomega))
    domega = omega[1] - omega[0]
    for i in range(Niom):
        for j in range(Nomega):
            matrix[i, j] = K_real(omega_n[i], omega[j])*domega
    return matrix

def K_matrix_imag(omega_n, omega): 
    Niom = len(omega_n)
    Nomega = len(omega)
    domega = omega[1] - omega[0]
    matrix = np.zeros((Niom, Nomega))
    for i in range(Niom):
        for j in range(Nomega):
            matrix[i, j] = K_imag(omega_n[i], omega[j])*domega
    return matrix
