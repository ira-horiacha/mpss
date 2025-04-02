import numpy as np
import matplotlib.pyplot as plt


# Функція для розв’язання рівняння методом Рунге-Кутта
def rk4(a, L, T, N, h, phi1, phi2, phi_y):
    dy = L / N  # Крок за глибиною
    dt = h  # Крок за часом
    M = int(T / dt)  # Кількість часових кроків

    # сітка для температури
    temperature = np.zeros((M + 1, N + 1))

    # Початкові та граничні умови
    for i in range(N + 1):
        temperature[0][i] = phi_y  # Початковий розподіл температури

    for t in range(M + 1):
        temperature[t][0] = phi1  # Ліва гранична умова
        temperature[t][-1] = phi2  # Права гранична умова

    # Функція для похідної температури
    def temperature_derivative(temp_array):
        d_temp = np.zeros_like(temp_array)  # Масив для збереження похідних
        for i in range(1, N):
            d_temp[i] = a * (temp_array[i + 1] - 2 * temp_array[i] + temp_array[i - 1]) / (dy ** 2)
        return d_temp

    # Чисельне розв’язання методом Рунге-Кутта
    for t in range(M):
        k1 = dt * temperature_derivative(temperature[t])
        k2 = dt * temperature_derivative(temperature[t] + 0.5 * k1)
        k3 = dt * temperature_derivative(temperature[t] + 0.5 * k2)
        k4 = dt * temperature_derivative(temperature[t] + k3)

        # Оновлення значень температури
        for i in range(1, N):
            temperature[t + 1][i] = temperature[t][i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6

    return temperature


# Функція для знаходження аналітичного розв’язку
def analytical_solution(a, L, T, N, h, alpha, beta, terms=30):
    t_steps = int(T / h) + 1  # Кількість кроків за часом
    y_steps = N + 1  # Кількість точок за глибиною

    t_values = np.linspace(0, T, t_steps)  # Масив часу
    y_values = np.linspace(0, L, y_steps)  # Масив глибини

    U_analytical = np.zeros((t_steps, y_steps))  # Матриця аналітичного розв’язку

    # Обчислення значень в кожній точці
    for t_idx in range(len(t_values)):
        for y_idx in range(len(y_values)):
            sum_series = 0
            for k in range(1, terms + 1):
                exp_term = -a * (np.pi * k / L) ** 2 * t_values[t_idx]
                exp_term = np.exp(max(exp_term, -100))  # Уникаємо проблем з експонентою

                # Додаємо повний множник (β (-1)^k - α)
                term = (1 / k) * (beta * (-1) ** k - alpha) * exp_term * np.sin(np.pi * k * y_values[y_idx] / L)

                sum_series += term

            U_analytical[t_idx][y_idx] = alpha + (beta - alpha) * y_values[y_idx] / L + 2 / np.pi * sum_series

    return U_analytical

# Функція для обчислення похибки
def calculate_errors(M, N, U_numerical, U_analytical):
    absolute_error = np.abs(U_numerical - U_analytical)
    relative_error = (1 / (M * N)) * np.sum((U_numerical - U_analytical) ** 2)
    return absolute_error, relative_error

# побудова графіка
def plot_surface(Y, T_vals, U, title, cmap):
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(Y, T_vals, U, cmap=cmap, alpha=0.8)
    ax.set_xlabel('Глибина (м)')
    ax.set_ylabel('Час (с)')
    ax.set_zlabel('Температура (°C)')
    ax.set_title(title)
    plt.show()


def main():
    # дані
    a = 111.0e-6  # Коефіцієнт теплопровідності
    L = 0.07  # Глибина стінки
    T = 1  # Час симуляції
    N = 100  # Кількість вузлів по глибині
    h = 0.15  # Часовий крок
    M = (T/h)-1
    alpha = 4  # Ліва гранична умова
    beta = 23  # Права гранична умова
    phi_y = 0  # Початкова температура

    # Обчислення чисельного та аналітичного розв’язку
    U_numerical = rk4(a, L, T, N, h, alpha, beta, phi_y)
    U_analytical = analytical_solution(a, L, T, N, h, alpha, beta, terms=30)

    # Обчислення похибок
    absolute_error, relative_error = calculate_errors(M, N, U_numerical, U_analytical)

    print("Максимальна абсолютна похибка:", np.max(absolute_error))
    print("Максимальна відносна похибка:", np.max(relative_error))

    # Побудова графіків
    y_positions = np.linspace(0, L, N + 1)
    t_positions = np.linspace(0, T, U_numerical.shape[0])
    Y, T_vals = np.meshgrid(y_positions, t_positions)

    plot_surface(Y, T_vals, U_numerical, 'Чисельний розв’язок', 'plasma')
    plot_surface(Y, T_vals, U_analytical, 'Аналітичний розв’язок', 'coolwarm')


if __name__ == "__main__":
    main()
