import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from utils.lu import lu_decompose, solve_lu


X = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
Y = [0.86603, 1.0, 0.86603, 0.50, 0.0, -0.50]


def build_mnk_system(X, Y, degree):
    N = len(X)

    S = [sum(X[j] ** p for j in range(N)) for p in range(2 * degree + 1)]

    A = [[S[i + k] for i in range(degree + 1)] for k in range(degree + 1)]
    b = [sum(Y[j] * (X[j] ** k) for j in range(N)) for k in range(degree + 1)]

    return A, b


def poly_value(coeffs, x):
    result = 0.0
    for c in coeffs[::-1]:
        result = result * x + c
    return result


# Ф
def sum_squared(coeffs, X, Y):
    N = len(X)
    return sum((poly_value(coeffs, X[j]) - Y[j]) ** 2 for j in range(N))


# MSE
def mse(coeffs, X, Y):
    return sum_squared(coeffs, X, Y) / len(X)


def print_polynom(coeffs, name):
    parts = []
    for i in range(len(coeffs)):
        a = coeffs[i]
        if i == 0:
            parts.append(f'{a:.6f}')
        else:
            sign = '+' if a >= 0 else '-'
            power = '*x' if i == 1 else f'*x^{i}'
            parts.append(f'{sign} {abs(a):.6f}{power}')
    print(f'{name}(x) = ' + ' '.join(parts))


def main():
    N = len(X)

    print('=====ИСХОДНЫЕ ДАННЫЕ=====')
    for i in range(len(X)):
        print(f'x_{i} = {X[i]:.2f}, y_{i} = {Y[i]:.5f}')
    print()

    degrees = [1, 2, 3]
    polys = {}
    mses = {}
    msds = {}

    for deg in degrees:
        print(f'=====МНОГОЧЛЕН СТЕПЕНИ {deg}=====')

        A, b = build_mnk_system(X, Y, deg)

        L, U, P, swaps = lu_decompose(A)
        coeffs = solve_lu(L, U, P, b)
        polys[deg] = coeffs

        print('Коэффициенты:')
        for i in range(len(coeffs)):
            print(f'a_{i} = {coeffs[i]:.6f}')
        print()

        print('Многочлен в явном виде:')
        print_polynom(coeffs, f'P_{deg}')
        print()

        f = sum_squared(coeffs, X, Y)
        m = mse(coeffs, X, Y)
        mses[deg] = m
        msds[deg] = f
        print(f'Ф = {f:.10f}    MSE_{deg} = {m:.10f}')
        print()

    print('=====СРАВНЕНИЕ MSE=====')
    mses_list = []
    for deg in degrees:
        MSE_deg = mses[deg]
        mses_list.append(MSE_deg)
        print(f'MSE_{deg} = {MSE_deg:.10f} Ф = {msds[deg]:.10f}')
    print()
    best_degree = min(mses, key=mses.get)
    best = mses[best_degree]
    print(f'Лучший приближающий многочлен: P_{best_degree}(x), MSE = {best:.10f} Ф = {msds[best_degree]:.10f}')

    plt.figure(figsize=(11, 7))
    plt.plot(X, Y, 'ko', markersize=10, zorder=5, label='Исходные данные')

    x_min, x_max = min(X), max(X)
    xs = [x_min + (x_max - x_min) * t / 200 for t in range(201)]

    colors = {1: 'tab:blue', 2: 'tab:orange', 3: 'tab:green'}
    for deg in [1, 2, 3]:
        ys = [poly_value(polys[deg], x) for x in xs]
        plt.plot(xs, ys, color=colors[deg], linewidth=2, label=f'P_{deg}(x), MSE = {mses[deg]:.6f}')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Метод наименьших квадратов')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig('img/3_3.png', dpi=150)
    plt.show()


if __name__ == '__main__':
    main()