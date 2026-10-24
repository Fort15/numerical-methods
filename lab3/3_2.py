import matplotlib.pyplot as plt
from utils.progonka import runoff_coef, reverse_move

X = [0.0, 1.0, 2.0, 3.0, 4.0]
F = [1.0, 0.86603, 0.5, 0.0, -0.5]
X_star = 1.5


def build_spline(X, F):
    n = len(X) - 1
    h = [X[i + 1] - X[i] for i in range(n)]

    a = F[:-1][:]
    b = [0.0] * n
    c = [0.0] * (n + 1)
    d = [0.0] * n

    if n == 1:
        b[0] = (F[1] - F[0]) / h[0]
        return a, b, c, d, h

    m = n - 1

    A = [0.0] * m  # поддиагональ
    B = [0.0] * m  # главная диагональ
    C = [0.0] * m  # наддиагональ
    D = [0.0] * m  # правая часть

    if m == 1:
        B[0] = 2.0 * (h[0] + h[1])
        D[0] = 3.0 * ((F[2] - F[1]) / h[1] - (F[1] - F[0]) / h[0])
    else:
        B[0] = 2.0 * (h[0] + h[1])
        C[0] = h[1]
        D[0] = 3.0 * ((F[2] - F[1]) / h[1] - (F[1] - F[0]) / h[0])

        for j in range(1, m - 1):
            i = j + 2
            A[j] = h[i - 2]
            B[j] = 2.0 * (h[i - 2] + h[i - 1])
            C[j] = h[i - 1]
            D[j] = 3.0 * ((F[i] - F[i - 1]) / h[i - 1] - (F[i - 1] - F[i - 2]) / h[i - 2])

        A[m - 1] = h[n - 2]
        B[m - 1] = 2.0 * (h[n - 2] + h[n - 1])
        D[m - 1] = 3.0 * ((F[n] - F[n - 1]) / h[n - 1] - (F[n - 1] - F[n - 2]) / h[n - 2])

    P, Q = runoff_coef(A, B, C, D, m)
    c_inner = reverse_move(P, Q, m)

    for j in range(m):
        c[j + 1] = c_inner[j]

    for i in range(n):
        b[i] = (F[i + 1] - F[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i])

    return a, b, c, d, h


def find_interval(x, X):
    n = len(X) - 1
    if x <= X[0]:
        return 0
    if x >= X[n]:
        return n - 1
    for i in range(n):
        if X[i] <= x <= X[i + 1]:
            return i
    return n - 1


def spline_value(x, X, a, b, c, d):
    i = find_interval(x, X)
    dx = x - X[i]
    return a[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3, i


def main():
    n = len(X) - 1

    print('=====ИСХОДНЫЕ ДАННЫЕ=====')
    for i in range(len(X)):
        print(f'x_{i} = {X[i]:.3f}, f_{i} = {F[i]:.5f}')
    print(f'X* = {X_star:.3f}')
    print()

    a, b, c, d, h = build_spline(X, F)
    print('=====КОЭФФИЦИЕНТЫ СПЛАЙНА=====')
    print(f"{'i':>3} | {'[x_{i-1}, x_i]':>18} | {'a_i':>10} | {'b_i':>10} | {'c_i':>10} | {'d_i':>10}")
    for i in range(n):
        interval = f"[{X[i]:.2f}, {X[i + 1]:.2f}]"
        print(f"{i + 1:>3} | {interval:>18} | {a[i]:>10.6f} | {b[i]:>10.6f} | "
              f"{c[i]:>10.6f} | {d[i]:>10.6f}")
    print()

    print('=====ИНТЕРВАЛЫ ДЕЙСТВИЯ КАЖДОГО МНОГОЧЛЕНА=====')
    for i in range(n):
        print(f"  Участок {i + 1}: x ∈ [{X[i]:.4f}, {X[i + 1]:.4f}]")
        print(f"    S{i + 1}(x) = {a[i]:.6f} + {b[i]:.6f}·(x-{X[i]:.4f}) "
              f"+ {c[i]:.6f}·(x-{X[i]:.4f})^2 + {d[i]:.6f}·(x-{X[i]:.4f})^3")
    print()

    print('=====ЗНАЧЕНИЕ СПЛАЙНА В ТОЧКЕ x*=====')
    S_star, idx_star = spline_value(X_star, X, a, b, c, d)
    print(f"  x* = {X_star}")
    print(f"  Номер интервала, содержащего x*: {idx_star + 1}")
    print(f"  Интервал: [{X[idx_star]:.4f}, {X[idx_star + 1]:.4f}]")
    print(f"  S(x*) = {S_star:.10f}")
    print()

    print('=====ПРОВЕРКА ПРОХОЖДЕНИЯ СПЛАЙНА ЧЕРЕЗ УЗЛЫ=====')
    print(f"{'x_i':>10} | {'f_i':>12} | {'S(x_i)':>14} | {'|f - S|':>12}")
    print("-" * 56)
    for i in range(len(X)):
        S_i, _ = spline_value(X[i], X, a, b, c, d)
        print(f"{X[i]:>10.4f} | {F[i]:>12.6f} | {S_i:>14.10f} | {abs(F[i] - S_i):>12.2e}")
    print()

    print('=====ВТОРАЯ ПРОИЗВОДНАЯ НА КОНЦАХ=====')
    S2_left = 2.0 * c[0]
    S2_right = 2.0 * c[n - 1] + 6.0 * d[n - 1] * h[n - 1]
    print(f"  Левый конец:   S''(x_0) = 2·c_1 = {S2_left:.6e}")
    print(f"  Правый конец:  S''(x_n) = 2·c_n + 6·d_n·h_n = {S2_right:.6e}")
    ok = abs(S2_left) < 1e-9 and abs(S2_right) < 1e-9
    print(f"  Граничные условия {('ВЫПОЛНЕНЫ' if ok else 'НАРУШЕНЫ')}")
    print()

    plt.figure(figsize=(11, 7))
    colors = plt.cm.tab10.colors
    for i in range(n):
        xs = [X[i] + (X[i + 1] - X[i]) * t / 100 for t in range(101)]
        ys = [a[i] + b[i] * (x - X[i]) + c[i] * (x - X[i]) ** 2 + d[i] * (x - X[i]) ** 3
              for x in xs]
        plt.plot(xs, ys, color=colors[i % len(colors)], linewidth=2,
                 label=f'S{i + 1}(x) на [{X[i]:.2f}, {X[i + 1]:.2f}]')

    plt.plot(X, F, 'ko', markersize=10, zorder=5, label='Узлы')
    plt.plot(X_star, S_star, 'r*', markersize=20, zorder=6,
             label=f'S(x*)={S_star:.4f} в x*={X_star}')

    plt.xlabel('x')
    plt.ylabel('S(x)')
    plt.title('Кубический сплайн')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()