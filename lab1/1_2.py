from utils.io_helpers import print_vector


def read_tridiagonal_matrix(path):
    """
    Читает трёхдиагональную матрицу из файла
    Формат
    n
    a1 ... an (поддиагональ, a1 = 0)
    b1 ... bn (главная диагональ)
    c1 ... cn (наддиагональ, cn = 0)
    d1 ... dn (правая часть)
    Возвращает (a, b, c, d, n)
    """
    with open(path, 'r', encoding='utf-8') as f:
        tokens = f.read().split()

    n = int(tokens[0])
    if len(tokens) != 1 + 4 * n:
        raise ValueError(f'Введено неверное количество чисел: {len(tokens)}')

    a = [float(x) for x in tokens[1:n + 1]]
    b = [float(x) for x in tokens[n + 1:2 * n + 1]]
    c = [float(x) for x in tokens[2 * n + 1:3 * n + 1]]
    d = [float(x) for x in tokens[3 * n + 1:4 * n + 1]]

    return a, b, c, d, n


def runoff_coef(a, b, c, d, n):
    P = [-c[0] / b[0]]
    Q = [d[0] / b[0]]

    for i in range(1, n):
        Pi = -c[i] / (b[i] + a[i] * P[i - 1])
        Qi = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1])

        P.append(Pi)
        Q.append(Qi)

    return P, Q


def reverse_move(P, Q, n):
    x = [0.0] * n
    x[-1] = Q[-1]
    for i in range(n-2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]

    return x


def main():
    a, b, c, d, n = read_tridiagonal_matrix('data/1.2.txt')

    P, Q = runoff_coef(a, b, c, d, n)

    for i in range(n):
        print(f'P{i + 1} = {P[i]:.4f}\nQ{i + 1} = {Q[i]:.4f}')
    print()

    x = reverse_move(P, Q, n)
    print_vector(x, title='Итоговое решение СЛАУ:')


if __name__ == '__main__':
    main()