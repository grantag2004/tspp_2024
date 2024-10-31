#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <random>
#include <cstdlib>

const int MAX_LEN = 50000;

void merge(std::vector<int>& mas, int left, int mid, int right) 
{
    int len1 = mid - left + 1;
    int len2 = right - mid;
    std::vector<int> mas1(len1);
    std::vector<int> mas2(len2);

    for (int i = 0; i < len1; ++i)
        mas1[i] = mas[left + i];
    for (int i = 0; i < len2; ++i)
        mas2[i] = mas[mid + i + 1];

    int i = 0, j = 0, k = left;
    while (i < len1 && j < len2) 
        if (mas1[i] <= mas2[j]) 
            mas[k++] = mas1[i++];
        else 
            mas[k++] = mas2[j++];

    while (i < len1) 
        mas[k++] = mas1[i++];
    while (j < len2) 
        mas[k++] = mas2[j++];
}

int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

void parallel_sort(std::vector<int>& mas, int left, int right, int k = 0) 
{
    if (right - left + 1 <= MAX_LEN) 
        std::qsort(&mas[left], right - left + 1, sizeof(int), [](const void* a, const void* b) {
            return (*(int*)a - *(int*)b);
        });
    else 
    {
        int mid = left + (right - left) / 2;

        if (k < omp_get_max_threads()) 
        {
            #pragma omp task shared(mas)
            parallel_sort(mas, left, mid, k + 1);
            #pragma omp task shared(mas)
            parallel_sort(mas, mid + 1, right, k + 1);
            #pragma omp taskwait
        } 
        else 
        {
            parallel_sort(mas, left, mid, k + 1);
            parallel_sort(mas, mid + 1, right, k + 1);
        }

        merge(mas, left, mid, right);
    }
}


int main(int argc, char* argv[]) 
{
    if (argc < 2) 
        return 1;
    
    int N = std::atoi(argv[1]); 

    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<> dis(0, 1000);

    //initialization of mas
    std::vector<int> mas1(N);
    for (int i = 0; i < N; i++)
        mas1[i] = dis(gen);

    std::vector<int> mas2 = mas1;
    double start = omp_get_wtime();
    qsort(&mas2[0], mas2.size(), sizeof(int), compare);
    double end = omp_get_wtime();
    double time1 = end - start;
    std::cout << "Qsort time: " << time1 << " seconds\n";

    std::vector<int> mas3 = mas1;
    double time2;
    for (int p = 1; p <= 16; ++p) 
    {
        omp_set_num_threads(p);
        mas3 = mas1;
        start = omp_get_wtime();
        #pragma omp parallel
        {
            #pragma omp single
            parallel_sort(mas3, 0, N - 1);
        }
        end = omp_get_wtime();
        time2 = end - start;
        std::cout << p << " pthreads" << std::endl;
        
        std::cout << "Parallel merge sort time: " << time2 << " seconds\n";

        std::cout << "Parallel time / std::sort time: " << (time2 / time1) * 100 << "% of std::sort\n";

        //Match check
        bool is_correct = true;
        for (int i = 0; i < N; ++i)
            if (mas3[i] != mas2[i]) 
            {
                is_correct = false;
                break;
            }

        if (is_correct) 
            std::cout << "The parallel merge sort matches qsort.\n";
        else 
            std::cout << "The parallel merge sort does not match qsort.\n";
    }

    return 0;
}
