import kernel
import printFile
import numpy as np

def gauss(x, mu, sigma):
    return 1.0/(sigma*np.sqrt(2*np.pi))*np.exp(-(x - mu)**2/(2*sigma**2))

def spectral(x):
    return 0.36*gauss(x, -2.5, 0.45) + 0.36*gauss(x, 2.5, 0.45) + 0.28*gauss(x, 0, 0.65)

def readA(spectralFileName):
    import os
    assert(os.path.exists(spectralFileName))
    omega = []
    A = []
    ifile = open(spectralFileName, "r")
    for (index, string) in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return np.asarray(omega), np.asarray(A)

def generateGreenFunction(omega_n, omega, A):
    KReal = kernel.K_matrix_real(omega_n, omega)
    KImag = kernel.K_matrix_imag(omega_n, omega)
    return KReal.dot(A), KImag.dot(A)

def main():
    import sys
    import os

    if (len(sys.argv) != 3):
        print "Niom = sys.argv[1], Nomega = sys.argv[2]. "
        return -1

    Niom = int(sys.argv[1])
    Nomega = int(sys.argv[2])

    s = []
    if (os.path.exists("s.txt")):
        ifile = open("s.txt", "r")
        for (index, string) in enumerate(ifile):
            s.append(float(string))
        ifile.close()
    
    if (len(s) != 0):
        Niom = len(s)

    omega_n = np.zeros(Niom)
    omega = np.zeros(Nomega)
    A = np.zeros(Nomega)
    omegaLower = -8.0
    omegaUpper = 8.0
    domega = (omegaUpper - omegaLower)/float(Nomega - 1)
    for i in range(Nomega):
        omega[i] = omegaLower + i*domega
        A[i] = spectral(omega[i])
    
    printFile.printFile(omega, A, "A.txt")

    beta = 10.0
    for i in range(Niom):
        omega_n[i] = (2*i + 1)*np.pi/beta

    GReal, GImag = generateGreenFunction(omega_n, omega, A)

    sigma = 0.001
    if (len(s) != 0):
        noise_real = np.asarray(map(lambda ele: np.random.normal(0, ele), s))
        noise_imag = np.asarray(map(lambda ele: np.random.normal(0, ele), s))
    else:
        noise_real = np.random.normal(0, sigma, Niom)
        noise_imag = np.random.normal(0, sigma, Niom)
    GReal += noise_real
    GImag += noise_imag
    ofile = open("noise.txt", "w")
    ofile.write(str(sigma) + "\n")
    ofile.close()

    printFile.printFile(omega_n, GReal, "G_real.txt")
    printFile.printFile(omega_n, GImag, "G_imag.txt")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
