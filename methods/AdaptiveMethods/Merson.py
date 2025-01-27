
import numpy as np
import matplotlib.pyplot as plt

def R1(x, y, m):
  return pow(pow(x + m, 2) + pow(y,2), 3/2)

def R2(x, y, M):
  return pow(pow(x - M, 2) + pow(y,2), 3/2)

def system(t, vars, M, m):
    """
    Правая часть системы уравнений.
    vars = [x, y, u, v], где u = x', v = y'.
    """
    x, y, u, v = vars

    # Уравнения
    dx_dt = u
    dy_dt = v
    du_dt = x + 2 * v - M * (x + m) / R1(x,y,m) - m * (x - M) / R2(x,y,M)
    dv_dt = y - 2 * u - M * y / R1(x,y,m) - m * y / R2(x,y,M)

    return np.array([dx_dt, dy_dt, du_dt, dv_dt])

def runge_kutta_4(system, t0, tf, vars0, h, params):
    """
    Решение системы ОДУ методом Рунге-Кутты 4-го порядка.
    """
    t_values = [t0]
    vars_values = [vars0]

    t = t0
    vars = np.array(vars0)

    while t < tf:
        if t + h > tf:  # Последний шаг
            h = tf - t

        k1 = system(t, vars, *params)
        k2 = system(t + h / 2, vars + h * k1 / 2, *params)
        k3 = system(t + h / 2, vars + h * k2 / 2, *params)
        k4 = system(t + h, vars + h * k3, *params)

        vars += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t += h

        t_values.append(t)
        vars_values.append(vars)

    return np.array(t_values), np.array(vars_values)

# Параметры задачи
M = 0.9
m = 0.1

# Начальные условия
x0 = 0.994
y0 = 0.0
u0 = 0.0  # u = x'
v0 = -2.031732629557337  # v = y'
vars0 = [x0, y0, u0, v0]

# Интервал времени и шаг
t0 = 0.0
tf = 10.0  # Время интегрирования
h = 0.0001  # Шаг

# Решение
params = (M, m)
t_values, vars_values = runge_kutta_4(system, t0, tf, vars0, h, params)

# Извлечение результатов
x_values = vars_values[:, 0]
y_values = vars_values[:, 1]

# Построение графиков
plt.figure(figsize=(10, 5))
plt.plot(x_values, y_values, label="x(y)")
plt.grid()
plt.savefig("pythsol.png")
