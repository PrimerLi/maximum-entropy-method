import numpy as np

def generate_matrix(dimension):
    assert(dimension >= 3)
    result = np.zeros((dimension, dimension))
    for i in range(dimension):
        for j in range(dimension):
            if (i == j - 2):
                result[i, j] = -1
            elif (i == j + 2):
                result[i, j] = -1
            elif(i == j):
                result[i, j] = 2
            else:
                pass
    return result

def print_matrix(matrix):
    (row, col) = matrix.shape
    for i in range(row):
        for j in range(col):
            print matrix[i,j], 
        print 

def main():
    import sys
    if (len(sys.argv) != 2):
        print "dimension = sys.argv[1]. "
        return -1

    dimension = int(sys.argv[1])
    matrix = generate_matrix(dimension)
    print_matrix(matrix)
    eigenvalues = np.linalg.eigvals(matrix)
    print "\n".join(map(str, sorted(eigenvalues)))
    print "Det = " + str(np.linalg.det(matrix))
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
