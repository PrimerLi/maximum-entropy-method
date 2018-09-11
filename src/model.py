import numpy as np
import printFile

def gauss(x, mu, sigma):
    return 1.0/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def generateModel(mu, sigma):
    import sys
    omega = []
    model = []
    Nomega = 501
    omega_lower = -20
    omega_upper = -omega_lower
    domega = (omega_upper - omega_lower)/float(Nomega-1)
    for i in range(Nomega):
        omega.append(omega_lower + i*domega)
        model.append(gauss(omega[i], mu, sigma))

    printFile.printFile(omega, model, "model.txt")
    
def main():
    import sys
    generateModel(0, 1)
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
