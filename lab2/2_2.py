from math import sqrt


def f1(x1, x2, a):
    return (x1 ** 2 + a ** 2) * x2 - a ** 3


def f2(x1, x2, a):
    return (x1 - a / 2) ** 2 + (x2 - a / 2) ** 2 - a ** 2


# Якобиан
def j11(x1, x2, a):
    return 2 * x1 * x2


def j12(x1, x2, a):
    return x1 ** 2 + a ** 2


def j21(x1, x2, a):
    return 2 * (x1 - a / 2)


def j22(x1, x2, a):
    return 2 * (x2 - a / 2)


# phi для простой итерации
def phi1(x2, a):
    return a / 2 + sqrt(a ** 2 - (x2 - a / 2) ** 2)


def phi2(x1, a):
    return a ** 3 / (x1 ** 2 + a ** 2)


def newton_method(x1, x2, a, eps=0.01):
    max_iter = 10000

    for k in range(max_iter):
        F1 = f1(x1, x2, a)
        F2 = f2(x1, x2, a)

        a11 = j11(x1, x2, a)
        a12 = j12(x1, x2, a)
        a21 = j21(x1, x2, a)
        a22 = j22(x1, x2, a)

        detJk = a11 * a22 - a21 * a12
        if abs(detJk) < 1e-15:
            print('detJ = 0')
            break

        detA1k = F1 * a22 - F2 * a12
        detA2k = a11 * F2 - a21 * F1

        x1_new = x1 - detA1k / detJk
        x2_new = x2 - detA2k / detJk

        if max(abs(x1_new - x1), abs(x2_new - x2)) < eps:
            return x1_new, x2_new

        x1, x2 = x1_new, x2_new

    print('too small count iter')
    return x1, x2


def simple_iterations(x1, x2, a, q, eps=1e-4):
    max_iter = 10000
    
    for k in range(max_iter):
        x1_new = phi1(x2, a)
        x2_new = phi2(x1, a)

        err = max(abs(x1_new - x1), abs(x2_new - x2))
        if (q / (1 - q)) * err <= eps:
            return x1_new, x2_new

        x1, x2 = x1_new, x2_new

    print('too small count iter')
    return x1, x2


def start_x(a):
    if a == 2:
        return 3.0, 0.6
    elif a == 3:
        return 4.5, 0.95
    elif a == 4:
        return 6.0, 1.25


# max||phi'(x)|| <= q < 1
def q_needed(a):
    if a == 2:
        return 0.39
    elif a == 3:
        return 0.35
    elif a == 4:
        return 0.33


def main():
    eps = 1e-4

    for a in [2, 3, 4]:
        print(f"---- a = {a} ----")
        print()
        x1_start, x2_start = start_x(a)

        print("     МЕТОД НЬЮТОНА")
        x_newton = newton_method(x1_start, x2_start, a, eps)
        print(f"x1 = {x_newton[0]:.5f}, x2 = {x_newton[1]:.5f}")
        print()

        print("     МЕТОД ПРОСТЫХ ИТЕРАЦИЙ")
        q = q_needed(a)
        x_simp = simple_iterations(x1_start, x2_start, a, q, eps)
        print(f'x1 = {x_simp[0]:.5f}, x2 = {x_simp[1]:.5f}, q = {q}')
        print()


if __name__ == '__main__':
    main()