import matplotlib.pyplot as plt
import math

X0 = 0.0
X_END = 1
H = 0.1

# y(0) = 1, y'(0) = 0
Y_0 = [1.0, 0.0]

# y' = f(x, y, z) = z
# z' = g(x, y, z) = 2cosx - y
def f(x, y, z):
    return z


def g(x, y, z):
    return 2 * math.cos(x) - y


def euler_method(h):
    n = round((X_END - X0) / h)

    x = [X0]
    y = [Y_0[0]]
    z = [Y_0[1]]

    for i in range(n):
        x_next = x[i] + h
        y_next = y[i] + h * f(x[i], y[i], z[i])
        z_next = z[i] + h * g(x[i], y[i], z[i])

        x.append(x_next)
        y.append(y_next)
        z.append(z_next)

    return x, y, z


def runge_kutta_step(x, y, z, h):
    k1 = h * f(x, y, z)
    l1 = h * g(x, y, z)

    k2 = h * f(x + h / 2, y + k1 / 2, z + l1 / 2)
    l2 = h * g(x + h / 2, y + k1 / 2, z + l1 / 2)

    k3 = h * f(x + h / 2, y + k2 / 2, z + l2 / 2)
    l3 = h * g(x + h / 2, y + k2 / 2, z + l2 / 2)

    k4 = h * f(x + h, y + k3, z + l3)
    l4 = h * g(x + h, y + k3, z + l3)

    y_next = y + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    z_next = z + 1 / 6 * (l1 + 2 * l2 + 2 * l3 + l4)

    return x + h, y_next, z_next


def runge_kutta_method(h):
    n = round((X_END - X0) / h)

    x = [X0]
    y = [Y_0[0]]
    z = [Y_0[1]]

    for i in range(n):
        x_next, y_next, z_next = runge_kutta_step(x[i], y[i], z[i], h)

        x.append(x_next)
        y.append(y_next)
        z.append(z_next)
    return x, y, z


def adams_method(h):
    n = round((X_END - X0) / h)

    x = [X0 + i * h for i in range(n + 1)]
    y = [0.0] * (n + 1)
    z = [0.0] * (n + 1)

    y[0] = Y_0[0]
    z[0] = Y_0[1]

    for i in range(3):
        _, y[i + 1], z[i + 1] = runge_kutta_step(x[i], y[i], z[i], h)

    for k in range(3, n):
        g0 = g(x[k], y[k], z[k])
        g1 = g(x[k - 1], y[k - 1], z[k - 1])
        g2 = g(x[k - 2], y[k - 2], z[k - 2])
        g3 = g(x[k - 3], y[k - 3], z[k - 3])

        f0 = f(x[k], y[k], z[k])
        f1 = f(x[k - 1], y[k - 1], z[k - 1])
        f2 = f(x[k - 2], y[k - 2], z[k - 2])
        f3 = f(x[k - 3], y[k - 3], z[k - 3])

        y[k + 1] = y[k] + h / 24 * (55 * f0 - 59 * f1 + 37 * f2 - 9 * f3)
        z[k + 1] = z[k] + h / 24 * (55 * g0 - 59 * g1 + 37 * g2 - 9 * g3)
    return x, y, z


def exact_solution(x):
    return x * math.sin(x) + math.cos(x)


# def exact_derivative(x):
#     return x * math.cos(x)


def runge_romberg(coarse, fine, p):
    errors = []

    for i in range(len(coarse)):
        error = abs(fine[2 * i] - coarse[i]) / (2 ** p - 1)
        errors.append(error)

    return errors


def print_results(name, x, y, y_half, p):
    print()
    print('=' * 100)
    print(name)
    print('=' * 100)

    print(f'{"x":>6} | {"y числ.":>15} | {"y точн.":>15} | {"Погрешность y":>17}')

    errors = []

    for i in range(len(x)):
        y_exact = exact_solution(x[i])
        error = abs(y[i] - y_exact)
        errors.append(error)

        print(f'{x[i]:6.2f} | {y[i]:15.9f} | {y_exact:15.9f} | {error:17.9e}')

    rr = runge_romberg(y, y_half, p)
    print()

    print(f'Максимальная погрешность y по точному решению = {max(errors)}')
    print(f'Максимальная погрешность по Рунге-Ромбергу = {max(rr)}')


def plot_results(x, y_euler, y_rk, y_adams):
    plt.figure(figsize=(11, 7))

    plt.plot(x, y_euler, color='blue', marker='o', markersize=10, label='Метод Эйлера')
    plt.plot(x, y_rk, color='black', marker='o', markersize=10, label='Метод Рунге-Кутты 4ого порядка')
    plt.plot(x, y_adams, color='green', marker='o', markersize=10, label='Метод Адамса 4ого порядка')

    x_exact = [X0 + i * (X_END - X0) / 300 for i in range(301)]
    y_exact = [exact_solution(x1) for x1 in x_exact]

    plt.plot(x_exact, y_exact, color='red', linewidth=2, label='Точное решение')
    plt.xlabel('x')
    plt.ylabel('y(x)')

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.show()


def main():
    print("y'' + y - 2*cos(x) = 0")
    print('y(0) = 1')
    print("y'(0) = 0")
    print('x принадлежит [0; 1]')
    print('h = 0.1')

    x_euler, y_euler, z_euler = euler_method(H)
    x_rk, y_rk, z_rk = runge_kutta_method(H)
    x_adams, y_adams, z_adams = adams_method(H)

    x_euler2, y_euler2, z_euler2 = euler_method(H / 2)
    x_rk2, y_rk2, z_rk2 = runge_kutta_method(H / 2)
    x_adams2, y_adams2, z_adams2 = adams_method(H / 2)

    print_results('Метод Эйлера', x_euler, y_euler, y_euler2, 1)
    print_results('Метод Рунге-Кутты 4ого порядка', x_rk, y_rk, y_rk2, 4)
    print_results('Метод Адамса 4ого порядка', x_adams, y_adams, y_adams2, 4)

    plot_results(x_euler, y_euler, y_rk, y_adams)


if __name__ == '__main__':
    main()