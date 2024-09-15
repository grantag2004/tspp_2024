#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

unsigned long int number_of_segments, number_of_pthreads;
double answer;
pthread_t *pthreads_vector;

typedef struct
{
    unsigned long int start;
    unsigned long int end;
    double step;
    double part_of_pi;
} info;

info *pthreads_info;


void *F(void *arg)
{
    info *inf = (info*) arg;
    inf -> part_of_pi = 0.0;
    for (unsigned long int i = inf -> start;i < inf -> end;++i)
    {
        double x = i * inf -> step + 0.5 * inf -> step;
        inf -> part_of_pi += 4.0 / (1.0 + x * x);    
    }

    return NULL;
}

int main(int argc, char **argv)
{
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    number_of_segments = atoi(argv[1]);
    number_of_pthreads = atoi(argv[2]);

    answer = 0.0;
    pthreads_info = (info *) calloc(number_of_pthreads, sizeof(info));
    pthreads_vector = (pthread_t *) calloc(number_of_pthreads, sizeof(pthread_t));
    unsigned long int steps_for_each_pthread = number_of_segments / number_of_pthreads;
    unsigned long int remaining_steps = number_of_segments % number_of_pthreads;

    for (unsigned long int i = 0;i < number_of_pthreads;++i)
    {
        pthreads_info[i].start = i * steps_for_each_pthread;
        pthreads_info[i].end = (i + 1) * steps_for_each_pthread;
        if (i == number_of_pthreads - 1)
            pthreads_info[i].end += remaining_steps;
            
        pthreads_info[i].step = 1.0 / number_of_segments;
        
        if (pthread_create(&pthreads_vector[i], NULL, F, &pthreads_info[i]) != 0)
        {
            printf("Error during creating %lu pthread\n", i);   
            return 1;         
        }
    }
    
    for (unsigned long int i = 0;i < number_of_pthreads;++i)
    {
        if (pthread_join(pthreads_vector[i], NULL) != 0)
        {
            printf("Error during joining %lu pthread\n", i); 
            return 1;
        }
        answer += pthreads_info[i].part_of_pi;
    }
    answer *= (1.0 / number_of_segments);

    clock_gettime(CLOCK_REALTIME, &end);

    double start_sec = start.tv_sec + start.tv_nsec / 1e9;
    double end_sec = end.tv_sec + end.tv_nsec / 1e9;
    double time = (double) (end_sec - start_sec);

    printf("%lf\n", answer);   
    printf("Elapsed time: %lf s\n", time);
    return 0;  
}