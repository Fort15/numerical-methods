import math


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


def main():
    Y = [f(x) for x in X]

    print('=====ИСХОДНЫЕ ДАННЫЕ=====')
    for i in range(len(X)):
        print(f'x_{i} = {X[i]:.3f}, f_{i} = {Y[i]:.3f}')
    print(f'X* = {X_star:.3f}')
    print()

    print('=====МНОГОЧЛЕН ЛАГРАНЖА=====')
    print_lagrange(X, Y)

    L_star = Ln(X, Y, X_star)
    f_star = f(X_star)
    print(f'L(X*) = {L_star:.10f}')
    print(f'f(X*) = {f_star:.10f}')
    print(f'Абсолютная погрешность интерполяции: |f(X*) - L(X*)| = {abs(f_star - L_star):.6f}')


if __name__ == '__main__':
    main()