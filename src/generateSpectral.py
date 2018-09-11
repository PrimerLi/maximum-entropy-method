import printFile
import numpy as np

def gauss(x, mu, sigma):
    return 1.0/(np.sqrt(2*np.pi)*sigma)*np.exp(-(x-mu)**2/(2*sigma**2))

def spectral(x):
    return 0.49*gauss(x, -1.1, 0.3) + 0.51*gauss(x, 1.1, 0.29)

def main():
    if (len(sys.argv) != 2):
        print "Nomega = sys.argv[1]. "
        return -1

    Nomega = int(sys.argv[1])
    omega = np.zeros(Nomega)
    A = np.zeros(Nomega)
    omegaLower = -5
    omegaUpper = 5
    domega = (omegaUpper - omegaLower)/float(Nomega-1)
    for i in range(Nomega):
        omega[i] = omegaLower + i*domega
        A[i] = spectral(omega[i])
    printFile.printFile(omega, A, "A.txt")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
