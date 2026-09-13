from utils.io_helpers import read_matrix, read_vector, print_matrix, print_vector
from utils.linalg_helper import matmul, matvec, matrix_equals


def lu_decompose(A):
    n = len(A)
    M = [row[:] for row in A]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    P = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    swaps = 0

    for k in range(n - 1):
        pivot_row = k
        max_abs = abs(M[k][k])
        for i in range(k + 1, n):
            if abs(M[i][k]) > max_abs:
                max_abs = abs(M[i][k])
                pivot_row = i

        if max_abs <= 1e-15:
            raise ValueError('Матрица вырождена: нулевой столбец')

        if pivot_row != k:
            M[k], M[pivot_row] = M[pivot_row], M[k]
            P[k], P[pivot_row] = P[pivot_row], P[k]

            for j in range(k):
                L[k][j], L[pivot_row][j] = L[pivot_row][j], L[k][j]
            swaps += 1

        pivot = M[k][k]
        for i in range(k + 1, n):
            mu = M[i][k] / pivot
            L[i][k] = mu
            for j in range(k, n):
                M[i][j] -= mu * M[k][j]

    return L, M, P, swaps


def solve_lu(L, U, P, b):
    """
    LU = PA -> Ax = b -> LUx = Pb
    Lz = Pb -> Ux = z
    """
    n = len(L)

    Pb = matvec(P, b)

    # 1) Lz = Pb
    z = [0.0] * n
    for i in range(n):
        s = Pb[i]
        for j in range(i):
            s -= L[i][j] * z[j]
        z[i] = s

    # 2) Ux = z
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = z[i]
        for j in range(i + 1, n):
            s -= U[i][j] * x[j]
        x[i] = s / U[i][i]

    return x


def inverse_lu(L, U, P, n):
    A_inv = [[0.0] * n for _ in range(n)]

    for i in range(n):
        e_i = [0.0] * n
        e_i[i] = 1.0

        x_i = solve_lu(L, U, P, e_i)
        for j in range(n):
            A_inv[j][i] = x_i[j]
    return A_inv


def det_lu(U, swaps):
    n = len(U)
    det = 1.0
    for i in range(n):
        det *= U[i][i]
    if swaps % 2 == 1:
        det = -det
    return det


def main():
    A = read_matrix('data/1.1_A.txt')
    b = read_vector('data/1.1_b.txt')
    n = len(A)

    print_matrix(A, title='Исходная матрица A:')
    print_vector(b, title='Вектор b:')

    L, U, P, swaps = lu_decompose(A)

    print_matrix(L, title='Матрица L:')
    print_matrix(U, title='Матрица U:')
    print_matrix(P, title='Матрица перестановок P:')

    LU = matmul(L, U)
    PA = matmul(P, A)
    print_matrix(LU, title='L * U:')

    x = solve_lu(L, U, P, b)
    print_vector(x, title='Решение системы x:')

    A_inv = inverse_lu(L, U, P, n)
    print_matrix(A_inv, title='A^(-1):')

    detA = det_lu(U, swaps)
    print(f'detA = {detA:.6f}')

    A_Ainv = matmul(A, A_inv)
    E = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    if matrix_equals(A_Ainv, E):
        print(f'A * A^(-1) = E')
    else:
        print(f'Обратная матрица найдена неверно')

    if matrix_equals(LU, PA):
        print(f'L * U = P * A')
    else:
        print(f' L * U != P * A')


if __name__ == '__main__':
    main()