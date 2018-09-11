import kernel
import numpy as np

def readA(fileName):
    omega = []
    A = []
    ifile = open(fileName, "r")
    for i, string in enumerate(ifile):
        string = string.replace(",", " ")
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return omega, A

def K(tau, omega, beta):
    import numpy as np
    return -np.exp(-tau*omega)/(1 + np.exp(-beta*omega))

def printFile(x, y, output):
    ofile = open(output, "w")
    for i in range(len(x)):
        ofile.write(str(x[i]) + "   " + str(y[i]) + "\n")
    ofile.close()

def main():
    import os
    import sys

    if (len(sys.argv) != 3):
        print "beta = sys.argv[1], fileName = sys.argv[2]. "
        return -1 

    omega_n = []
    ifile = open("G_cc_real.txt", "r")
    for index, string in enumerate(ifile):
        omega_n.append(float(string.split()[0]))
    ifile.close()

    beta = float(sys.argv[1])
    fileName = sys.argv[2]

    omega, A = readA(fileName)
    domega = omega[1] - omega[0]
    tau = []
    dtau = 0.1
    ntau = int(beta/dtau) + 1
    for i in range(ntau):
        tau.append(i*dtau)

    A = np.asarray(A)
    G = []
    for i in range(len(tau)):
        KA = 0
        for nw in range(len(omega)):
            KA = KA + domega*K(tau[i], omega[nw], beta)*A[nw]
        G.append(KA)

    Giom_real = kernel.K_matrix_real(omega_n, omega).dot(A)
    Giom_imag = kernel.K_matrix_imag(omega_n, omega).dot(A)
    printFile(omega_n, Giom_real, "Giom_real_test.txt")
    printFile(omega_n, Giom_imag, "Giom_imag_test.txt")
    printFile(tau, G, "G_tau.txt")
    return 0

main()
