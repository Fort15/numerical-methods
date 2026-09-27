from math import log, sqrt


def f(x):
    return 2 ** x - x ** 2 - 0.5


def df_dx(x):
    return 2 ** x * log(2) - 2 * x


def d2f_dx2(x):
    return 2 ** x * log(2) ** 2 - 2


def dphi_dx(x):
    return 2 ** x * log(2) / (2 * sqrt(2 ** x - 0.5))


def phi(x):
    return sqrt(2 ** x - 0.5)


def newton_method(x0, eps=0.01):
    max_iter = 10000
    xk = x0
    for k in range(max_iter):
        f_xk = f(xk)
        df_xk = df_dx(xk)
        xk_1 = xk - f_xk / df_xk
        if abs(xk_1 - xk) < eps:
            return xk_1
        xk = xk_1

    print('too small count iter')
    return xk


def simple_iteration_method(x0, q, eps=0.01):
    max_iter = 10000
    xk = x0
    for k in range(max_iter):
        xk_1 = phi(xk)
        if (q / (1 - q)) * abs(xk_1 - xk) <= eps:
            return xk_1
        xk = xk_1

    print('too small count iter')
    return xk


# =0.47, 1.57
def main():
    eps = 0.001

    a_newton, b_newton = -0.6, 2

    print("     МЕТОД НЬЮТОНА")

    if f(b_newton) * d2f_dx2(b_newton) > 0:
        print(f"f({b_newton}) * f''({b_newton}) > 0")
    else:
        print(f"f({b_newton}) * f''({b_newton}) < 0")
        print('Теорема не выполняется, выберите другое начальное приближение')
        return

    x = newton_method(b_newton, eps)
    print(f"x = {x:.5f}")
    print()

    print("     МЕТОД ПРОСТОЙ ИТЕРАЦИИ")

    a_simp, b_simp = 1.4, 1.7
    q = dphi_dx(1.7) # возрастает на отрезке и < 1
    x0 = (a_simp + b_simp) / 2
    x = simple_iteration_method(x0, q, eps=eps)
    print(f"x = {x:.5f}")


if __name__ == '__main__':
    main()