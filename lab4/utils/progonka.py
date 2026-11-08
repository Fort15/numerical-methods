def runoff_coef(a, b, c, d, n):
    P = [-c[0] / b[0]]
    Q = [d[0] / b[0]]

    for i in range(1, n):
        Pi = -c[i] / (b[i] + a[i] * P[i - 1])
        Qi = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1])

        P.append(Pi)
        Q.append(Qi)

    return P, Q


def reverse_move(P, Q, n):
    x = [0.0] * n
    x[-1] = Q[-1]
    for i in range(n-2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]

    return x