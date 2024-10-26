import matplotlib.pyplot as plt

# Примерные данные (замените на реальные данные, полученные при измерениях)
N_values = [25000, 75000, 125000, 175000, 225000, 275000]
T_N = [0.894272, 2.53002, 4.13716, 5.72486, 7.43626, 8.99389]  # время выполнения для каждого N при фиксированном P
S_N = [0.139024, 0.132142, 0.128293, 0.128184, 0.128900, 0.126942]  # ускорение при фиксированном P
P = 8  # число потоков для графиков по N
E_N = [s / P for s in S_N]       # эффективность при фиксированном P

P_values = [1, 2, 4, 8, 16]
T_P = [44.7306, 22.5052, 11.2483, 5.80179, 3.05503]  # время выполнения для каждого P при фиксированном N
S_P = [T_P[0] / t for t in T_P]    # ускорение при фиксированном N
E_P = [s / p for s, p in zip(S_P, P_values)]  # эффективность

# Построение графиков для T(N), S(N), E(N)
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(N_values, T_N, marker='o', label="T(N)")
plt.xlabel("Число точек")
plt.ylabel("Время в секундах")
plt.title("Зависимость Т от N")
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(N_values, S_N, marker='o', color="green", label="S(N)")
plt.xlabel("Число точек")
plt.ylabel("Ускорение")
plt.title("Зависимость S от N")
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(N_values, E_N, marker='o', color="red", label="E(N)")
plt.xlabel("Число точек")
plt.ylabel("Эффективность распараллеливания")
plt.title("Зависимость E от N")
plt.grid(True)

plt.tight_layout()
plt.show()

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
