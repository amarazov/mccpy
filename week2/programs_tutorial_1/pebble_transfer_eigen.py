import numpy

neighbor = [[1, 3, 0, 0], [2, 4, 0, 1], [2, 5, 1, 2],
            [4, 6, 3, 0], [5, 7, 3, 1], [5, 8, 4, 2],
            [7, 6, 6, 3], [8, 7, 6, 4], [8, 8, 7, 5]]
transfer = numpy.zeros((9, 9))
for k in range(9):
    for neigh in range(4): transfer[neighbor[k][neigh], k] += 0.25
eigenvalues, eigenvectors = numpy.linalg.eig(transfer)
print(eigenvalues)

for iter in range(len(eigenvalues)):
    print("eigenvalue:", eigenvalues[iter])
    print(eigenvectors[:, iter] / sum(eigenvectors[:, iter]))
