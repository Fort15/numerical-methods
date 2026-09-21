from utils.io_helpers import read_matrix, print_matrix
from utils.linalg_helper import (mat_sub, vecvec, mat_scale,
                                 scalar_product, identity, matmul,
                                 matrix_equals, matrix_transposition)
import numpy as np


def sign(x):
    if x >= 0:
        return 1
    return -1


def householder_matrix(v):
    n = len(v)
    vvT = vecvec(v, v)
    vTv = scalar_product(v, v)
    E = identity(n)
    scaled = mat_scale(vvT, 2 / vTv)
    H = mat_sub(E, scaled)
    return H


def qr_decompose(A):
    n = len(A)
    Ak = [row[:] for row in A]
    Q = identity(n)
    for i in range(n - 1):
        v = [0.0] * n
        for j in range(n):
            if j < i:
                v[j] = 0
            elif i == j:
                v[j] = Ak[i][i] + sign(Ak[i][i]) * sum([Ak[z][i] ** 2 for z in range(i, n)]) ** 0.5
            else:
                v[j] = Ak[j][i]
        H = householder_matrix(v)
        Ak = matmul(H, Ak)
        Q = matmul(Q, H)

    return Q, Ak


def qr_eigenvalues(A, eps=1e-6):
    max_iter = 10000
    n = len(A)
    Ak = [row[:] for row in A]
    for k in range(max_iter):
        Q, R = qr_decompose(Ak)
        Ak = matmul(R, Q)
        if sum(Ak[i][j] ** 2 for i in range(n) for j in range(i)) ** 0.5 < eps:
            break
    return [Ak[i][i] for i in range(n)]


def main():
    A = read_matrix('data/1.5.txt')
    n = len(A)

    print_matrix(A, title='Исходная матрица A:')
    print()

    Q, R = qr_decompose(A)

    print_matrix(Q, title='Матрица Q:')
    print_matrix(R, title='Матрица R:')

    QR = matmul(Q, R)
    print_matrix(QR, title='Q * R:')
    if matrix_equals(QR, A):
        print("Q * R = A")
    else:
        print("Q * R != A")

    Qt = matrix_transposition(Q)
    QtQ = matmul(Qt, Q)
    E = identity(n)
    if matrix_equals(QtQ, E):
        print("Q^T*Q = E")
    else:
        print("Q^T*Q != E")
    print()

    eigvals = qr_eigenvalues(A)
    print(f'Собственные значения матрицы A:')
    for i in range(n):
        lam = eigvals[i]
        print(f'lambda_{i + 1} = {lam:.6f}')

    print()

    eigvals_np = np.linalg.eigvals(np.array(A, dtype=float))
    print('Проверка (NumPy):')
    for i in range(n):
        lam = eigvals_np[i]
        if abs(lam.imag) < 1e-10:
            print(f'lambda_{i + 1} = {lam.real:.6f}')
        else:
            print(f'lambda_{i + 1} = {lam.real:.6f} + {lam.imag:.6f}i')

    print()
    print('=' * 60)
    print('Пример матрицы с комплексно-сопряжёнными lambda:')

    B = [
        [1, 3, 1],
        [1, 1, 4],
        [4, 3, 1]
    ]
    print_matrix(B, title='Матрица B:')

    eigvals_B = np.linalg.eigvals(np.array(B, dtype=float))
    print('Собственные значения (NumPy):')
    for i in range(n):
        lam = eigvals_B[i]
        print(f'lambda_{i + 1} = {lam:.6f}')


if __name__ == '__main__':
    main()