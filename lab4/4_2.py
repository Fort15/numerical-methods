import math
import matplotlib.pyplot as plt

from utils.progonka import runoff_coef, reverse_move


X0 = 1.0
X_END = 2.0
Y_X0 = math.exp(-1)
Y_XEND = 0.5 * math.exp(-2)


def p(x):
    return 2.0 / x


def q(x):
    return -1.0


def rhs_f(x):
    return 0.0


def exact_solution(x):
    return math.exp(-x) / x


def f(x, y, z):
    return z


def g(x, y, z):
    return y - 2 * z / x


def rk4_step(x, y, z, h):
    k1 = h * f(x, y, z)
    l1 = h * g(x, y, z)

    k2 = h * f(x + h / 2, y + k1 / 2, z + l1 / 2)
    l2 = h * g(x + h / 2, y + k1 / 2, z + l1 / 2)

    k3 = h * f(x + h / 2, y + k2 / 2, z + l2 / 2)
    l3 = h * g(x + h / 2, y + k2 / 2, z + l2 / 2)

    k4 = h * f(x + h, y + k3, z + l3)
    l4 = h * g(x + h, y + k3, z + l3)

    y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    z_next = z + (l1 + 2 * l2 + 2 * l3 + l4) / 6.0
    return y_next, z_next


def solve_inp(eta, h):
    n = round((X_END - X0) / h)
    x = X0
    y = Y_X0
    z = eta

    for _ in range(n):
        y, z = rk4_step(x, y, z, h)
        x += h
    return y


def shooting_method(h):
    n = round((X_END - X0) / h)
    x_vals = [X0 + i * h for i in range(n + 1)]
    y_vals = []

    def F(eta):
        return solve_inp(eta, h) - Y_XEND

    eta0 = -0.5
    eta1 = -1.0
    eps = 1e-8

    for _ in range(20):
        f0 = F(eta0)
        f1 = F(eta1)

        if abs(f1) < eps:
            break

        eta2 = eta1 - f1 * (eta1 - eta0) / (f1 - f0)
        eta0, eta1 = eta1, eta2

    best_eta = eta1

    y = Y_X0
    z = best_eta
    y_vals.append(y)
    x_current = X0

    for _ in range(n):
        y, z = rk4_step(x_current, y, z, h)
        x_current += h
        y_vals.append(y)

    return x_vals, y_vals


def finite_difference_method(h):
    n = round((X_END - X0) / h)
    x_vals = [X0 + i * h for i in range(n + 1)]

    num_unknowns = n - 1

    A = [0.0] * num_unknowns
    B = [0.0] * num_unknowns
    C = [0.0] * num_unknowns
    D = [0.0] * num_unknowns

    for k in range(num_unknowns):
        i = k + 1
        xi = x_vals[i]

        A[k] = 1.0 / (h ** 2) - p(xi) / (2.0 * h)
        B[k] = -2.0 / (h ** 2) + q(xi)
        C[k] = 1.0 / (h ** 2) + p(xi) / (2.0 * h)
        D[k] = rhs_f(xi)

    D[0] -= A[0] * Y_X0
    D[-1] -= C[-1] * Y_XEND

    P, Q = runoff_coef(A, B, C, D, num_unknowns)
    y_internal = reverse_move(P, Q, num_unknowns)

    y_vals = [Y_X0] + y_internal + [Y_XEND]
    return x_vals, y_vals


def runge_romberg(coarse, fine, p):
    errors = []
    for i in range(len(coarse)):
        error = abs(fine[2 * i] - coarse[i]) / (2 ** p - 1)
        errors.append(error)
    return errors


def print_table_and_errors(name, x, y, y_half, p_order):
    print("\n" + "=" * 85)
    print(f"--- {name} ---")
    print("=" * 85)
    print(f"{'x':>6} | {'y числ.':>15} | {'y точн.':>15} | {'Погрешность y':>20}")
    print("-" * 85)

    max_exact_error = 0.0
    for i in range(len(x)):
        y_exact = exact_solution(x[i])
        err = abs(y[i] - y_exact)
        if err > max_exact_error:
            max_exact_error = err
        print(f"{x[i]:6.2f} | {y[i]:15.8f} | {y_exact:15.8f} | {err:20.2e}")

    rr_errors = runge_romberg(y, y_half, p_order)
    max_rr_error = max(rr_errors)

    print("-" * 85)
    print(f"Максимальная погрешность y по точному решению:  {max_exact_error:.2e}")
    print(f"Максимальная оценка погрешности (Рунге-Ромберг, p={p_order}): {max_rr_error:.2e}")


def plot_results(x1, y1_sh, y1_fd, x2, y2_sh, y2_fd):
    plt.figure(figsize=(11, 7))

    plt.scatter(x1, y1_sh, color='tab:blue', s=40, label='shooting method, step=0.1', zorder=3)
    plt.scatter(x1, y1_fd, color='tab:green', s=40, label='finite difference method, step=0.1', zorder=3)

    plt.scatter(x2, y2_sh, color='tab:orange', s=15, label='shooting method, step=0.05', zorder=2)
    plt.scatter(x2, y2_fd, color='tab:red', s=15, label='finite difference method, step=0.05', zorder=2)

    x_ex = [X0 + i * (X_END - X0) / 200 for i in range(201)]
    y_ex = [exact_solution(xi) for xi in x_ex]
    plt.plot(x_ex, y_ex, color='black', linewidth=1.5, label='exact solution', zorder=1)

    plt.xlabel('x')
    plt.ylabel('y(x)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    H = 0.1

    x1, y1_shoot = shooting_method(H)
    _, y1_fd = finite_difference_method(H)

    x2, y2_shoot = shooting_method(H / 2)
    _, y2_fd = finite_difference_method(H / 2)

    print_table_and_errors("МЕТОД СТРЕЛЬБЫ (h=0.1)", x1, y1_shoot, y2_shoot, p_order=4)

    print_table_and_errors("КОНЕЧНО-РАЗНОСТНЫЙ МЕТОД (h=0.1)", x1, y1_fd, y2_fd, p_order=2)

    plot_results(x1, y1_shoot, y1_fd, x2, y2_shoot, y2_fd)


if __name__ == '__main__':
    main()
