import os


def read_matrix(path):
    '''
    Читает матрицу из файла по пути path
    Формат:
    n
    a11 a12 ... a1n
    a21 a22 ... a2n
    ...
    an1 an2 ... ann
    Возвращает двумерный список типа float
    '''

    if not os.path.isfile(path):
        raise FileNotFoundError(f'Файл не найден: {path}')

    with open(path, 'r', encoding='utf-8') as f:
        tokens = f.read().split()

    if not tokens:
        raise ValueError(f'Файл пустой: {path}')

    n = int(tokens[0])
    values = [float(x) for x in tokens[1:]]

    if len(values) != n * n:
        raise ValueError(f'Неверное количество чисел для матрицы {n} * {n}')

    matrix = [values[i * n:(i + 1) * n] for i in range(n)]
    return matrix


def read_vector(path):
    """
    Читает вектор из файла по пути path
    Формат:
    n
    b1
    ...
    bn
    Возвращает список чисел типа float
    """

    if not os.path.isfile(path):
        raise FileNotFoundError(f'Файл не найден: {path}')

    with open(path, 'r', encoding='utf-8') as f:
        tokens = f.read().split()

    if not tokens:
        raise ValueError(f'Файл пустой: {path}')

    n = int(tokens[0])
    values = [float(x) for x in tokens[1:]]

    if len(values) != n:
        raise ValueError(f'Неверное количество чисел для вектора длины {n}')

    return values


def print_matrix(matrix, title="Матрица"):
    if title:
        print(title)
    for row in matrix:
        for x in row:
            print(f"{x:10.4f}", end="")
        print()
    print()


def print_vector(v, title="Вектор"):
    if title:
        print(title)
    for x in v:
        print(f"{x:10.4f}")
    print()
