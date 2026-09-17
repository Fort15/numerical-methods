from utils.io_helpers import read_matrix, print_matrix
from utils.linalg_helper import matrix_transposition, matmul, matrix_equals
from math import pi, atan, sin, cos


def find_max_not_in_diag(A):
    n = len(A)
    a_ij_max = abs(A[1][0])
    i_max = 1
    j_max = 0
    for i in range(1, n):
        for j in range(i):
            if abs(A[i][j]) > a_ij_max:
                a_ij_max = abs(A[i][j])
                i_max = i
                j_max = j

    return a_ij_max, min(i_max, j_max), max(i_max, j_max)


def matrix_rotation(A, i, j):
    n = len(A)
    U = [[1.0 if i1 == j1 else 0.0 for j1 in range(n)] for i1 in range(n)]
    if A[i][i] == A[j][j]:
        phi = pi / 4
    else:
        phi = 1 / 2 * atan(2 * A[i][j] / (A[i][i] - A[j][j]))

    sinphi = sin(phi)
    cosphi = cos(phi)

    U[i][i] = cosphi
    U[i][j] = -sinphi
    U[j][i] = sinphi
    U[j][j] = cosphi

    return U


def t(A):
    s = 0
    for i in range(1, len(A)):
        for j in range(i):
            s += A[i][j] ** 2
    return s ** (1 / 2)


def jakobi(A, eps=0.3):
    max_iter = 10000
    n = len(A)
    Ak = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for k in range(max_iter):
        aij_max, i, j = find_max_not_in_diag(Ak)
        U = matrix_rotation(Ak, i, j)
        V = matmul(V, U)
        U_T = matrix_transposition(U)
        Ak = matmul(matmul(U_T, Ak), U)
        if t(Ak) < eps:
            break

    return V, Ak


def main():
    A = read_matrix('data/1.4.txt')
    eigenvectors, eigenvalues = jakobi(A, eps=0.01)

    print_matrix(eigenvectors, title='Матрица собственных векторов:')
    print_matrix(eigenvalues, title='Матрица собственных значений:')

    if matrix_equals(matmul(A, eigenvectors), matmul(eigenvectors, eigenvalues)):
        print('AV=VΛ')
    else:
        print('AV!=VΛ')


if __name__ == '__main__':
    main()