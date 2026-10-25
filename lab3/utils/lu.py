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
    LU = PA →  Lz = Pb → Ux = z
    Возвращает x.
    """
    n = len(L)

    Pb = [sum(P[i][j] * b[j] for j in range(n)) for i in range(n)]

    z = [0.0] * n
    for i in range(n):
        s = Pb[i]
        for j in range(i):
            s -= L[i][j] * z[j]
        z[i] = s

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = z[i]
        for j in range(i + 1, n):
            s -= U[i][j] * x[j]
        x[i] = s / U[i][i]

    return x