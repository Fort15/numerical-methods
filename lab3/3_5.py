import math


# аналитическое решение интеграла
I_exact = 1 / 9 * math.log(4) - 1 / 12

A = 0.0
B = 4.0
H1 = 1.0
H2 = 0.5


def f(x):
    return x / (3 * x + 4) ** 2


def left_rect(a, b, h):
    n = round((b - a) / h)
    s = 0.0
    for i in range(n):
        x = a + i * h
        s += f(x)
    return s * h


def right_rect(a, b, h):
    n = round((b - a) / h)
    s = 0.0
    for i in range(1, n + 1):
        x = a + i * h
        s += f(x)
    return s * h


def mid_rect(a, b, h):
    n = round((b - a) / h)
    s = 0.0
    for i in range(n):
        x = a + i * h + h / 2.0
        s += f(x)
    return s * h


def trapezoid(a, b, h):
    n = round((b - a) / h)
    s = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        x = a + i * h
        s += f(x)
    return s * h


def simpson(a, b, h):
    n = round((b - a) / h)
    if n % 2 == 1:
        raise ValueError('Для метода Симпсона число шагов должно быть чётным')

    s = f(a) + f(b)
    for i in range(1, n):
        x = a + i * h
        if i % 2 == 0:
            s += 2 * f(x)
        else:
            s += 4 * f(x)
    return s * h / 3.0


def runge_method(I1, I2, p):
    correction = (I2 - I1) / (2 ** p - 1)
    I_refined = I2 + correction

    error = abs(correction)
    return I_refined, error


def main():
    print(f'Точное значение интеграла: {I_exact:.6f}\n')

    methods = [
        ("Левые прямоуг.", left_rect, 1),
        ("Правые прямоуг.", right_rect, 1),
        ("Средние прямоуг.", mid_rect, 2),
        ("Трапеции", trapezoid, 2),
        ("Симпсон", simpson, 4)
    ]

    results = []

    print(
        f"{'Метод':<18} | {'I(h1)':>10} | {'I(h2)':>10} | {'p':^2} | {'Уточн. Рунге':>13} | {'Оценка погр.':>13} | {'Абс. погр.':>13}")
    print("-" * 95)

    for name, func, p in methods:
        I1 = func(A, B, H1)
        I2 = func(A, B, H2)

        I_ref, err = runge_method(I1, I2, p)
        abs_err = abs(I2 - I_exact)

        results.append((name, abs_err))

        print(f"{name:<18} | {I1:>10.6f} | {I2:>10.6f} | {p:^2} | {I_ref:>13.8f} | {err:>13.8f} | {abs_err:>13.8f}")
    print("-" * 95)

    best_name = ""
    min_error = float('inf')

    for name, err in results:
        if err < min_error:
            min_error = err
            best_name = name

    print("\n" + "=" * 50)
    print(f"Самый точный метод (по h2): {best_name}")
    print(f"Его абсолютная погрешность: {min_error:.8f}")


if __name__ == '__main__':
    main()