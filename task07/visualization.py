import matplotlib.pyplot as plt

P_values = [1, 2, 4, 8, 16, 32]
T_P = [63.2158, 32.606, 15.876, 9.62968, 4.37121, 2.56795]  # время выполнения для каждого P при фиксированном N
S_P = [T_P[0] / t for t in T_P]    # ускорение при фиксированном N
E_P = [s / p for s, p in zip(S_P, P_values)]  # эффективность

# Построение графиков для T(P), S(P), E(P)
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(P_values, T_P, marker='o', label="T(P)")
plt.xlabel("Число потоков")
plt.ylabel("Время в секундах")
plt.title("Зависимость Т от P")
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(P_values, S_P, marker='o', color="green", label="S(P)")
plt.xlabel("Число потоков")
plt.ylabel("Ускорение")
plt.title("Зависимость S от P")
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(P_values, E_P, marker='o', color="red", label="E(P)")
plt.xlabel("Число потоков")
plt.ylabel("Эффективность распараллеливания")
plt.title("Зависимость E от P")
plt.grid(True)

plt.tight_layout()
plt.show()
