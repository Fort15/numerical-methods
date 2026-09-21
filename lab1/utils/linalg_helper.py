def matmul(A, B):
    """Умножение матриц"""
    m = len(A)
    n = len(A[0])
    p = len(B[0])

    if len(B) != n:
        raise ValueError(f'Несовместные размерности {m} * {n} и {len(B)} * {p}')
    C = [[0.0] * p for _ in range(m)]

    for i in range(m):
        for j in range(p):
            s = 0.0
            for k in range(n):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C


def matvec(A, x):
    """Умножение матрицы на вектор"""
    m = len(A)
    n = len(A[0])

    if len(x) != n:
        raise ValueError(f'Несовместные размерности')

    result = [0.0] * m
    for i in range(m):
        s = 0.0
        for j in range(n):
            s += A[i][j] * x[j]
        result[i] = s
    return result


def vecvec(u, v):
    """Умножение вектора на вектор"""
    m = len(u)
    n = len(v)
    C = [[0.0] * n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            C[i][j] = u[i] * v[j]
    return C


def matrix_equals(A, B, eps=1e-10):
    m = len(A)
    n = len(A[0])
    for i in range(m):
        for j in range(n):
            if abs(A[i][j] - B[i][j]) > eps:
                return False
    return True


def matrix_transposition(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def mat_sub(A, B):
    """Вычитание матриц"""
    m = len(A)
    n = len(A[0])
    if len(B) != m or len(B[0]) != n:
        raise ValueError('Размеры матриц не совпадают')
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(m)]


def mat_scale(A, c):
    """Умножение матрицы на число"""
    m = len(A)
    n = len(A[0])
    return [[c * A[i][j] for j in range(n)] for i in range(m)]


def scalar_product(u, v):
    """Скалярное произведение векторов u и v."""
    if len(u) != len(v):
        raise ValueError('Разные длины векторов')
    s = 0.0
    for i in range(len(u)):
        s += u[i] * v[i]
    return s


def identity(n):
    """Единичная матрица n x n."""
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

