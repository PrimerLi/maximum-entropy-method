import matrix
import kernel
import model
import printFile

def main():
    import sys
    import os 
    import numpy as np
    
    omega_n = []
    ifile = open("G_cc_imag.txt", "r")
    for index, string in enumerate(ifile):
        omega_n.append(float(string.split()[0]))
    ifile.close()
    Niom = len(omega_n)

    G_real = np.zeros(Niom)
    G_imag = np.zeros(Niom)
    ifile = open("G_cc_real.txt", "r")
    for index, string in enumerate(ifile):
        G_real[index] = float(string.split()[1])
    ifile.close()
    ifile = open("G_cc_imag.txt", "r")
    for index, string in enumerate(ifile):
        G_imag[index] = float(string.split()[1])
    ifile.close()

    C_real = np.zeros((Niom, Niom))
    ifile = open("CM_cc_real.txt")
    for i in range(Niom):
        for j in range(Niom):
            string = ifile.readline()
            C_real[i,j] = float(string.split()[2])
    ifile.close()

    u, s, ut = np.linalg.svd(C_real)#C_real = u.dot(s).dot(ut)
    ofile = open("eigenvalues.txt", "w")
    for i in range(len(s)):
        ofile.write(str(s[i]) + "\n")
    ofile.close()
    ofile = open("eig_log.txt", "w")
    for i in range(len(s)):
        ofile.write(str(i+1) + "    " + str(np.log(s[i])) + "\n")
    ofile.close()

    ratio = 1.0e-7
    counter = 0
    for i in range(Niom):
        if (s[i] >= ratio*s[0]):
            counter = counter + 1

    ofile = open("s.txt", "w")
    for i in range(counter):
        ofile.write(str(s[i]) + "\n")
    ofile.close()

    u_truncated = np.zeros((counter, counter))
    for i in range(counter):
        for j in range(counter):
            u_truncated[i,j] = u[i,j]
    G_real_truncated = np.zeros(counter)
    G_imag_truncated = np.zeros(counter)
    for i in range(counter):
        G_real_truncated[i] = G_real[i]
        G_imag_truncated[i] = G_imag[i]
    G_real_rotated = np.transpose(u_truncated).dot(G_real_truncated)
    G_imag_rotated = np.transpose(u_truncated).dot(G_imag_truncated)

    ofile = open("G_real_rotated.txt", "w")
    for i in range(len(G_real_rotated)):
        ofile.write(str(omega_n[i]) + "   " + str(G_real_rotated[i]) + "\n")
    ofile.close()
    ofile = open("G_imag_rotated.txt", "w")
    for i in range(len(G_imag_rotated)):
        ofile.write(str(omega_n[i]) + "    " + str(G_imag_rotated[i]) + "\n")
    ofile.close()

    if (not os.path.exists("model.txt")):
        model.generateModel(0, 1)
    omega = []
    model_function = []
    ifile = open("model.txt", "r")
    for i, string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        model_function.append(float(a[1]))
    ifile.close()
    Nomega = len(omega)
   
    K_matrix_real = kernel.K_matrix_real(omega_n, omega)
    K_matrix_imag = kernel.K_matrix_imag(omega_n, omega)
    K_real_truncated = np.zeros((counter, Nomega))
    K_imag_truncated = np.zeros((counter, Nomega))
    for i in range(counter):
        for j in range(Nomega):
            K_real_truncated[i,j] = K_matrix_real[i,j]
            K_imag_truncated[i,j] = K_matrix_imag[i,j]
    K_real_rotated = np.transpose(u_truncated).dot(K_real_truncated)
    K_imag_rotated = np.transpose(u_truncated).dot(K_imag_truncated)
    printFile.printMatrix(K_real_rotated, "K_real_rotated.txt")
    printFile.printMatrix(K_imag_rotated, "K_imag_rotated.txt")       
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
