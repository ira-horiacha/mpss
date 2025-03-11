import numpy as np
import matplotlib.pyplot as plt


### 1. ПЕРША СИСТЕМА: ХИЖАКИ-ЖЕРТВИ

# Вхідні параметри для моделі
a11 = 0.04
a12 = 0.0004
a21 = 0.0004
a22 = 0.16

x0 = 960
y0 = 660

h = 0.1
T = 150

t_values = []
t = 0
i = 0
while i <= T:
    t_values.append(t)
    t += h
    i += 1

#масиви для збереження обчислених значень
x_values = [0] * len(t_values)
y_values = [0] * len(t_values)

# рівняння системи
def f1(x, y):
    return a11 * x - a12 * x * y

def f2(x, y):
    return a22 * y - a21 * x * y

#метод Рунге-Кутта 4-го порядку
def lotka_volterra_rk4():

    x_values[0] = x0
    y_values[0] = y0

    for i in range(1, len(t_values)):
        x = x_values[i - 1]
        y = y_values[i - 1]

        # Обчислюємо коефіцієнти
        k1x = h * f1(x, y)
        k1y = h * f2(x, y)

        k2x = h * f1(x + k1x / 2, y + k1y / 2)
        k2y = h * f2(x + k1x / 2, y + k1y / 2)

        k3x = h * f1(x + k2x / 2, y + k2y / 2)
        k3y = h * f2(x + k2x / 2, y + k2y / 2)

        k4x = h * f1(x + k3x, y + k3y)
        k4y = h * f2(x + k3x, y + k3y)

        # Оновлюємо значення
        x_values[i] = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 6
        y_values[i] = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6

def plot():
    # Побудова графіків результатів
    plt.figure(figsize=(12, 5))

    # Графік зміни популяції жертв у часі
    plt.subplot(1, 3, 1)
    plt.plot(t_values, x_values, label="Жертви (x)", color='blue')
    plt.xlabel("Час (дні)")
    plt.ylabel("Популяція")
    plt.legend()
    plt.grid()

    # Графік зміни популяції хижаків у часі
    plt.subplot(1, 3, 2)
    plt.plot(t_values, y_values, label="Хижаки (y)", color='pink')
    plt.xlabel("Час (дні)")
    plt.legend()
    plt.grid()

    # Фазова траєкторія (залежність хижаків від жертв)
    plt.subplot(1, 3, 3)
    plt.plot(x_values, y_values, label="y(x)", color='purple')
    plt.xlabel("Жертви (x)")
    plt.ylabel("Хижаки (y)")
    plt.legend()
    plt.grid()

    # Відображаємо всі графіки
    plt.tight_layout()
    plt.show()

def print_results():
    print(f"{'t':<10}{'x':<10}{'y':<10}")
    print("-" * 30)
    for i in range(len(t_values)):
        print(f"{t_values[i]:<10.2f}{x_values[i]:<10.2f}{y_values[i]:<10.2f}")


### 2. ДРУГА СИСТЕМА: ЕПІДЕМІЯ

# Вхідні параметри
H = 996 #кількість людей
beta = 21 #здорові
gamma = 4 #кількість днів потрібних для одужання

x0_2 = 896
y0_2 = 86
z0_2 = H - x0_2 - y0_2

h2 = 0.1
T2 = 40

t_values_2 = []
t = 0
i = 0
while i <= T2:
    t_values_2.append(t)
    t += h
    i += 1

# Масиви для збереження результатів
x_values_2 = np.zeros(len(t_values_2))
y_values_2 = np.zeros(len(t_values_2))
z_values_2 = np.zeros(len(t_values_2))

# Рівняння системи
def g1(x, y):
    return - (beta / H) * x * y

def g2(x, y):
    return (beta / H) * x * y - (1 / gamma) * y

def g3(y):
    return (1 / gamma) * y

# Метод Рунге-Кутта 4-го порядку
def infection_model_rk4():
    x_values_2[0] = x0_2
    y_values_2[0] = y0_2
    z_values_2[0] = z0_2

    for i in range(1, len(t_values_2)):
        x = x_values_2[i - 1]
        y = y_values_2[i - 1]
        z = z_values_2[i - 1]

        k1x = h2 * g1(x, y)
        k1y = h2 * g2(x, y)
        k1z = h2 * g3(y)

        k2x = h2 * g1(x + k1x / 2, y + k1y / 2)
        k2y = h2 * g2(x + k1x / 2, y + k1y / 2)
        k2z = h2 * g3(y + k1y / 2)

        k3x = h2 * g1(x + k2x / 2, y + k2y / 2)
        k3y = h2 * g2(x + k2x / 2, y + k2y / 2)
        k3z = h2 * g3(y + k2y / 2)

        k4x = h2 * g1(x + k3x, y + k3y)
        k4y = h2 * g2(x + k3x, y + k3y)
        k4z = h2 * g3(y + k3y)

        x_values_2[i] = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 6
        y_values_2[i] = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6
        z_values_2[i] = z + (k1z + 2 * k2z + 2 * k3z + k4z) / 6


### ПОБУДОВА ГРАФІКІВ

def plot_2():
    plt.figure(figsize=(12, 6))

    # Друга система
    plt.subplot(2, 3, 4)
    plt.plot(t_values_2, x_values_2, label="x(t)", color='green')
    plt.xlabel("Час")
    plt.legend()
    plt.grid()

    plt.subplot(2, 3, 5)
    plt.plot(t_values_2, y_values_2, label="y(t)", color='orange')
    plt.xlabel("Час")
    plt.legend()
    plt.grid()

    plt.subplot(2, 3, 6)
    plt.plot(t_values_2, z_values_2, label="z(t)", color='black')
    plt.xlabel("Час")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()

def print_results_2():
    print(f"{'t':<10}{'x':<10}{'y':<10}{'z':<10}")
    print("-" * 40)
    for i in range(len(t_values_2)):
        print(f"{t_values_2[i]:<10.2f}{x_values_2[i]:<10.2f}{y_values_2[i]:<10.2f}{z_values_2[i]:<10.2f}")



def main():

    lotka_volterra_rk4()
    print_results()
    plot()

    print("\n\n\n")

    infection_model_rk4()
    print_results_2()
    plot_2()


if __name__ == "__main__":
    main()