def printFile(x, y, fileName):
    ofile = open(fileName, "w")
    for i in range(len(x)):
        ofile.write(str(x[i]) + "    " + str(y[i]) + "\n")
    ofile.close()

def printMatrix(matrix, fileName):
    shape = matrix.shape
    ofile = open(fileName, "w")
    for i in range(shape[0]):
        for j in range(shape[1]):
            ofile.write(str(i) + "    " + str(j) + "    " + str(matrix[i,j]) + "\n")
    ofile.close()
