import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


f = math.cos
X = [0, math.pi / 6, 2 * math.pi / 6, 3 * math.pi / 6]
# X = [0, math.pi / 6, 5 * math.pi / 12, math.pi / 2]
X_star = math.pi / 4


# Значение многочлена Лагранжа в x
def Ln(X, Y, x):
    n = len(X)
    result = 0.0
    for i in range(n):
        tmp = Y[i]
        for j in range(n):
            if i != j:
                tmp *= (x - X[j]) / (X[i] - X[j])
        result += tmp
    return result


# Многочлен Лагранжа
def print_lagrange(X, Y):
    n = len(X)
    print(f'L(x) = ', end='')
    first = True
    for i in range(n):
        denom = 1
        for j in range(n):
            if i != j:
                denom *= (X[i] - X[j])
        coef = Y[i] / denom

        if abs(coef) < 1e-12:
            continue

        if not first:
            print(' + ' if coef >= 0 else ' - ', end='')
            print(f'{abs(coef):.3f}', end='')
        else:
            print(f'{coef:.3f}', end='')
            first = False

        for j in range(n):
            if i != j:
                print(f'(x - {X[j]:.3f})', end='')
        print()
        if i < n - 1:
            print('      ', end='')
    print()


def divided_differences(X, Y):
    n = len(X)
    table = [[0.0] * n for _ in range(n)]
    for i in range(n):
        table[i][0] = Y[i]
    for k in range(1, n):
        for i in range(n - k):
            table[i][k] = (table[i + 1][k - 1] - table[i][k - 1]) / (X[i + k] - X[i])
    return table


def newton_value(X, Y, x):
    table = divided_differences(X, Y)
    coeffs = [table[0][k] for k in range(len(X))]
    result = coeffs[-1]
    for k in range(len(X) - 2, -1, -1):
        result = result * (x - X[k]) + coeffs[k]
    return result


def print_divided_differences(X, Y):
    n = len(X)
    table = divided_differences(X, Y)
    print('Таблица разделённых разностей')
    header = '  i |   x_i   |' + ''.join(f'  пор.{k:<2}   |' for k in range(n))
    print(header)
    print('-' * len(header))

    for i in range(n):
        row = f'  {i} | {X[i]:>9.6f} |'
        for k in range(n):
            if k <= n - 1 - i:
                row += f' {table[i][k]:>10.6f} |'
            else:
                row += ' ' * 12 + '|'
        print(row)
    print()

    return table


def plot_graph(X, Y, X_star, f, a=None, b=None):
    if a is None:
        a = min(X) - 1
    if b is None:
        b = max(X) + 1

    N = 500
    xs = [a + (b - a) * i / N for i in range(N + 1)]

    ys_f = [f(x) for x in xs]
    ys_L = [Ln(X, Y, x) for x in xs]
    ys_P = [newton_value(X, Y, x) for x in xs]

    plt.figure(figsize=(10,6))
    plt.plot(xs, ys_f, 'k-', linewidth=2, label='f(x) - точная функция')
    plt.plot(xs, ys_L, 'r--', linewidth=1.5, label='L(x) - Лагранж')
    plt.plot(xs, ys_P, 'b:', linewidth=1.5, label='P(x) - Ньютон')

    plt.plot(X, Y, 'go', markersize=8, label='узлы интерполяции')

    f_star = f(X_star)
    plt.plot(X_star, f_star, 'm*', markersize=15, label=f'X* = {X_star:.4f}')

    plt.axvline(X_star, color='gray', linestyle='--', alpha=0.5)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Интерполяционные многочлены Лагранжа и Ньютона')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('img/3_1.png', dpi=150)
    plt.show()


def print_newton(X, Y, table):
    n = len(X)
    coeffs = [table[0][k] for k in range(n)]
    print('P(x) =', end='')
    print(f'{coeffs[0]:.5f}', end='')
    for k in range(1, n):
        print(f" {'+' if coeffs[k] >= 0 else '-'} {abs(coeffs[k]):.3f}", end='')
        for j in range(k):
            print(f'(x - {X[j]:.3f})', end='')
        if k < n - 1:
            print()
            print('      ', end='')
    print()
    print()
    return coeffs


def main():
    Y = [f(x) for x in X]

    print('=====ИСХОДНЫЕ ДАННЫЕ=====')
    for i in range(len(X)):
        print(f'x_{i} = {X[i]:.3f}, f_{i} = {Y[i]:.3f}')
    print(f'X* = {X_star:.3f}')
    print()

    print('=====МНОГОЧЛЕН ЛАГРАНЖА=====')
    print_lagrange(X, Y)

    print('=====МНОГОЧЛЕН НЬЮТОНА=====')
    table = print_divided_differences(X, Y)
    print_newton(X, Y, table)

    L_star = Ln(X, Y, X_star)
    f_star = f(X_star)
    P_star = newton_value(X, Y, X_star)
    print('=====СРАВНЕНИЕ МЕТОДОВ=====')
    print(f'L(X*) = {L_star:.10f}')
    print(f'P(X*) = {P_star:.10f}')
    print(f'f(X*) = {f_star:.10f}')
    print()
    print(f'|f(X*) - L(X*)| = {abs(f_star - L_star):.6f}')
    print(f'|f(X*) - P(X*)| = {abs(f_star - P_star):.6f}')
    print(f'|P(X*) - L(X*)| = {abs(P_star - L_star):.6f}')

    # ГРАФИК
    plot_graph(X, Y, X_star, f)


if __name__ == '__main__':
    main()