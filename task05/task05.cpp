#include <iostream>
#include <omp.h>
#include <random>

#define MAX_STEPS 10000

int main() 
{
    int a, b, x, N, P, steps, position;
    double p;
    std::cout << "Enter a and b" << std::endl;
    std::cin >> a >> b;
    std::cout << "Enter p" << std::endl;
    std::cin >> p;
    std::cout << "Enter x" << std::endl;
    std::cin >> x;
    std::cout << "Enter N" << std::endl;
    std::cin >> N;
    std::cout << "Enter P" << std::endl;
    std::cin >> P;

    int get_b = 0;         
    double total_lifetime = 0;

    double start_time = omp_get_wtime();

    #pragma omp parallel num_threads(P) default(none) private(steps, position) shared(a, b, x, N, p, P) reduction(+:get_b, total_lifetime)
    {
        std::mt19937 generator(omp_get_thread_num());  
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) 
        {
            position = x;
            steps = 0;

            while (position != a && position != b && steps < MAX_STEPS) {
                steps++;
                if (distribution(generator) < p)
                    position += 1;  
                else
                    position -= 1;  
            }

            if (position == b) 
                ++get_b;

            total_lifetime += steps;
        }
    }

    double end_time = omp_get_wtime();
    double time = end_time - start_time;

    // Вывод результатов
    std::cout << "Вероятность достижения b: " << static_cast<double>(get_b) / N << std::endl;
    std::cout << "Среднее время жизни одной частицы: " << total_lifetime / N << " шагов" << std::endl;
    std::cout << "Время выполнения основного цикла: " << time << " секунд" << std::endl;

    return 0;
}
