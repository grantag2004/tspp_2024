import matplotlib.pyplot as plt

P_values = [i for i in range(1, 17)]
T_P = [19.1178, 11.2088, 6.87324, 5.16196, 4.47874, 4.45119, 4.11927, 4.05589, 3.67259, 3.11033, 3.56383, 3.84294, 3.2908, 2.8619, 3.38021, 2.77889]  # время выполнения для каждого P при фиксированном N
S_P = [T_P[0] / t for t in T_P]    # ускорение при фиксированном N
E_P = [s / p for s, p in zip(S_P, P_values)]  # эффективность

# Построение графиков для T(P), S(P), E(P)
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(P_values, T_P, marker='o', label="T(P)")
plt.axhline(y=17.035, color='blue', linestyle='--', label="y = 17.035s")
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
