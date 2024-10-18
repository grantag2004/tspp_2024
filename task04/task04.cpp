#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <arm_neon.h>

//MULTIPLICATION WITH A ROW VIEW(Row-Major)

//Serial version
void serial_multiply(float* A, float* B, float* C, int N) 
{
    for (int i = 0; i < N; ++i) 
    {
        for (int j = 0; j < N; ++j) 
        {
            float curr_sum = 0.0f;
            for (int k = 0; k < N; ++k) 
                curr_sum += A[i * N + k] * B[k * N + j];
            C[i * N + j] = curr_sum;
        }
    }
}

//Vectorize version
void vectorize_multiply(float* A, float* B, float* C, int N) 
{
    for (int i = 0; i < N; ++i) 
    {
        for (int j = 0; j < N; j += 4) 
        { 
            float32x4_t c_vec = vdupq_n_f32(0.0f); 
            
            for (int k = 0; k < N; ++k) 
            {
                float32x4_t b_vec = vld1q_f32(&B[k * N + j]);  
                float32x4_t a_val = vdupq_n_f32(A[i * N + k]);  
                
                c_vec = vmlaq_f32(c_vec, a_val, b_vec);        
            }

            vst1q_f32(&C[i * N + j], c_vec); 
        }
    }
}

void initialize(float* matrix, int N) 
{
    for (int i = 0; i < N * N; i++) 
        matrix[i] = static_cast<float>(rand()) / RAND_MAX; 
}

void test(int N) 
{
    float* A = new float[N * N];
    float* B = new float[N * N];
    float* C_serial = new float[N * N];
    float* C_vectorize = new float[N * N];

    
    initialize(A, N);
    initialize(B, N);

    auto start = std::chrono::high_resolution_clock::now();
    serial_multiply(A, B, C_serial, N);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = end - start;
    std::cout << "Serial_multiply time: " << time1.count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    vectorize_multiply(A, B, C_vectorize, N);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time2 = end - start;
    std::cout << "Vectorize_multiply time: " << time2.count() << " seconds" << std::endl;

    delete[] A;
    delete[] B;
    delete[] C_serial;
    delete[] C_vectorize;
}


int main() {
    int sizes[] = {512, 1024, 2048};

    for (int N: sizes) 
    {
        std::cout << "Matrix size: " << N << "x" << N << std::endl;
        test(N);
    }

    return 0;
}