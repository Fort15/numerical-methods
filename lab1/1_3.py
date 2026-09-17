from utils.io_helpers import read_matrix, read_vector, print_vector
import numpy as np


def norm_matrix(A):
    """Норма матрицы по строкам"""
    n = len(A)
    return max([sum([abs(A[i][j]) for j in range(n)]) for i in range(n)])


def simple_iterations(alpha, beta, eps=0.01):
    max_iter = 10000
    n = len(beta)
    x = beta[:]

    for k in range(1, max_iter + 1):
        x_new = [0.0] * n
        for i in range(n):
            s = beta[i]
            for j in range(n):
                if i != j:
                    s += alpha[i][j] * x[j]
            x_new[i] = s

        diff = max([abs(x_new[i] - x[i]) for i in range(n)])
        x = x_new

        if diff < eps:
            return x, k

    return x, max_iter


def zeidel(alpha, beta, eps=0.01):
    max_iter = 10000
    n = len(beta)
    x = beta[:]

    for k in range(1, max_iter + 1):
        x_old = x[:]
        for i in range(n):
            s = beta[i]
            for j in range(n):
                if i != j:
                    s += alpha[i][j] * x[j]
            x[i] = s

        diff = max([abs(x[i] - x_old[i]) for i in range(n)])

        if diff < eps:
            return x, k

    return x, max_iter


def create_alpha_beta(A, b):
    """
    Ax=b -> x = alpha * x + beta
    """
    n = len(A)

    for i in range(n):
        # нулевой диагональный элемент -> перестановка строк
        if abs(A[i][i]) < 1e-15:
            for k in range(i + 1, n):
                if abs(A[k][i]) > 1e-15:
                    A[i], A[k] = A[k], A[i]
                    b[i], b[k] = b[k], b[i]
                    break
            else:
                raise ValueError(f'Нулевой столбец {i + 1}')

    alpha = [[0.0] * n for _ in range(n)]
    beta = [0.0] * n
    for i in range(n):
        beta[i] = b[i] / A[i][i]
        for j in range(n):
            if i != j:
                alpha[i][j] = -A[i][j] / A[i][i]

    return alpha, beta


def main():
    A = read_matrix('data/1.3_A.txt')
    b = read_vector('data/1.3_b.txt')
    n = len(A)

    eps = 0.01
    print(f'Заданная точность eps = {eps}')
    print()

    alpha, beta = create_alpha_beta(A, b)

    if norm_matrix(alpha) < 1:
        print(f'||alpha|| < 1, условие сходимости выполнено')
    else:
        print(f'||alpha|| >= 1, сходимость не гаранитрована')
    print()

    x_simp_it, count_simp_it = simple_iterations(alpha, beta, eps=eps)
    print(f'Метод простой итерации:')
    print(f'Число итераций: {count_simp_it}')
    print_vector(x_simp_it, title='Решение x:')
    print()

    x_zeid, count_zeid_it = zeidel(alpha, beta, eps=eps)
    print(f'Метод Зейделя:')
    print(f'Число итераций: {count_zeid_it}')
    print_vector(x_zeid, title='Решение x:')
    print()

    print('Сравнение методов:')
    if count_zeid_it < count_simp_it:
        print(f'Зейдель быстрее в {count_simp_it / count_zeid_it:.2f} раз')
    elif count_zeid_it > count_simp_it:
        print(f'Простые итерации быстрее в {count_zeid_it / count_simp_it:.2f} раз')
    else:
        print('Одинаковое число итераций')

    x_expected = np.linalg.solve(np.array(A), np.array(b))
    x_expected = [float(v) for v in x_expected]
    print_vector(x_expected, title='Эталонное решение:')
    diff_simp = max([abs(x_simp_it[i] - x_expected[i]) for i in range(n)])
    diff_zeid = max([abs(x_zeid[i] - x_expected[i]) for i in range(n)])
    print(f'Отклонение от эталона:')
    print(f'Простые итерации: {diff_simp:.3f}')
    print(f'Метод Зейделя: {diff_zeid:.3f}')


if __name__ == '__main__':
    main()