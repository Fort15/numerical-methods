import math

X_star = 1.0
X = [-1.0, 0.0, 1.0, 2.0, 3.0]
Y = [-0.5, 0.0, 0.5, 0.86603, 1.0]


def find_interval(X, x_star):
    n = len(X)
    for i in range(n - 1):
        if X[i] <= x_star <= X[i + 1]:
            return i
    return None


def choose_nodes(X, x_star):
    n = len(X)
    idx = find_interval(X, x_star)
    if idx is None:
        raise ValueError('x* не лежит в таблице')

    if idx + 2 > n - 1:
        idx = n - 3
    return idx


# f[x_i, x_{i+1}]
def divided_diff_1(X, Y, i):
    return (Y[i+1] - Y[i]) / (X[i+1] - X[i])


def divided_diff_2(X, Y, i):
    d1 = divided_diff_1(X, Y, i)
    d2 = divided_diff_1(X, Y, i + 1)
    return (d2 - d1) / (X[i + 2] - X[i])


def first_derivative(X, Y, x_star, i):
    d1 = divided_diff_1(X, Y, i)
    d2 = divided_diff_2(X, Y, i)
    return d1 + d2 * (2 * x_star - X[i] - X[i + 1])


def second_derivative(X, Y, i):
    return 2 * divided_diff_2(X, Y, i)


def main():
    print('=====ИСХОДНА ТАБЛИЦА ЗНАЧЕНИЙ=====')
    print(f"{'i':>3} | {'x_i':>10} | {'y_i':>12}")
    print('-' * 32)
    for i in range(len(X)):
        print(f"{i:>3} | {X[i]:>10.5f} | {Y[i]:>12.5f}")
    print()

    print(f'Точка в которой считаем производную: x* = {X_star}')
    print()

    i = choose_nodes(X, X_star)

    print('=====УЗЛЫ ИСПОЛЬЗОВАННЫЕ ДЛ АППРОКСИМАЦИИ=====')
    print(f'Используются 3 узла: x_{i}, x_{i + 1}, x_{i + 2}')
    print(f'x_{i} = {X[i]:.4f}, y_{i} = {Y[i]:.6f}')
    print(f'x_{i + 1} = {X[i + 1]:.4f}, y_{i + 1} = {Y[i + 1]:.6f}')
    print(f'x_{i + 2} = {X[i + 2]:.4f}, y_{i + 2} = {Y[i + 2]:.6f}')
    print()

    y1 = first_derivative(X, Y, X_star, i)
    y2 = second_derivative(X, Y, i)

    print('=====ЗНАЧЕНИЯ ПРОИЗВОДНЫХ=====')
    print(f"y'(x*) = {y1:.6f}")
    print(f"y'(x*) = {y2:.6f}")
    print()


if __name__ == '__main__':
    main()