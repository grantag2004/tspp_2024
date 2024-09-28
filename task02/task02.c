#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

int SIZE;
int WRITERS;
int WRITERS_ACTIONS;
int READERS;
int READERS_ACTIONS;

//POD-type = int 
typedef struct 
{
    int *queue;
    size_t max_size;
    size_t curr_size;
    size_t head;
    size_t tail;
    pthread_cond_t is_empty;
    pthread_cond_t is_full;
    pthread_mutex_t mutex;
    pthread_mutex_t print_mutex; 
    int writers_done; 
    int readers_done; 
} MyConcurrentQueue;

void create_queue(MyConcurrentQueue *q, size_t size)
{
    q -> queue = (int *) malloc(size * sizeof(int));
    q -> max_size = size;
    q -> curr_size = 0;
    q -> head = 0;
    q -> tail = 0;
    q -> writers_done = 0;
    q -> readers_done = 0;
    pthread_cond_init(&(q -> is_empty), NULL);
    pthread_cond_init(&(q -> is_full), NULL);
    pthread_mutex_init(&(q -> mutex), NULL);
    pthread_mutex_init(&(q -> print_mutex), NULL); 
}

int put(MyConcurrentQueue *q, int element) 
{
    pthread_mutex_lock(&(q -> mutex));
    //Waiting until the queue is no longer the maximum size, the mutex release
    while ((q -> curr_size) == (q -> max_size) && !(q -> readers_done)) 
        pthread_cond_wait(&(q -> is_full), &(q -> mutex));

    if (q -> curr_size == q -> max_size && q->readers_done) 
    {
        pthread_mutex_unlock(&(q->mutex));
        return -1; 
    }

    q -> queue[q -> tail] = element;
    q -> tail = (q -> tail + 1) % q -> max_size;//cyclically tranfer to the beginning
    ++(q -> curr_size);

    //Signal for a thread waiting for an element to appear
    pthread_cond_signal(&(q -> is_empty));
    pthread_mutex_unlock(&(q -> mutex));
    return 0;
}

int get(MyConcurrentQueue *q) 
{
    pthread_mutex_lock(&(q -> mutex));
    while (q -> curr_size == 0 && !q->writers_done) 
        pthread_cond_wait(&(q -> is_empty), &(q -> mutex));

    if (q -> curr_size == 0 && q->writers_done) 
    {
        pthread_mutex_unlock(&(q->mutex));
        return -1; 
    }

    int element = q -> queue[q -> head];
    q -> head = (q -> head + 1) % q -> max_size;
    --(q -> curr_size);

    //Signal for a thread waiting for an element to be removed
    pthread_cond_signal(&(q -> is_full));
    pthread_mutex_unlock(&(q -> mutex));
    return element;
}

void delete_queue(MyConcurrentQueue *q)
{
    free(q -> queue);
    pthread_cond_destroy(&(q -> is_empty));
    pthread_cond_destroy(&(q -> is_full));
    pthread_mutex_destroy(&(q -> mutex));
    pthread_mutex_destroy(&(q->print_mutex)); 
}

void* writer(void* arg) 
{
    MyConcurrentQueue* queue = (MyConcurrentQueue*)arg;
    for (int i = 0; i < WRITERS_ACTIONS; ++i) 
    {
        int r = rand(); // Returns a pseudo-random integer between 0 and RAND_MAX.
        int flag = put(queue, r);
        if (flag == -1)
            break;
        pthread_mutex_lock(&(queue->print_mutex));
        printf("Write: %d\n", r);
        pthread_mutex_unlock(&(queue->print_mutex));
        usleep(100000);
    }

    pthread_mutex_lock(&(queue->mutex));
    queue-> writers_done = 1; 
    pthread_cond_broadcast(&(queue->is_empty));
    pthread_mutex_unlock(&(queue->mutex));
    return NULL;
}

void* reader(void* arg) 
{
    MyConcurrentQueue* queue = (MyConcurrentQueue*)arg;
    for (int i = 0; i < READERS_ACTIONS; ++i) {
        int element = get(queue);
        if (element == -1) 
            break;
        pthread_mutex_lock(&(queue->print_mutex));
        printf("Read: %d\n", element);
        pthread_mutex_unlock(&(queue->print_mutex));
        usleep(150000);
    }
    pthread_mutex_lock(&(queue->mutex));
    queue-> readers_done = 1; 
    pthread_cond_broadcast(&(queue->is_full));
    pthread_mutex_unlock(&(queue->mutex));
    return NULL;
}

int test(int a, int b, int c, int d, int e)
{
    SIZE = a;
    WRITERS = b;
    WRITERS_ACTIONS = c;
    READERS = d;
    READERS_ACTIONS = e;
    MyConcurrentQueue queue;
    create_queue(&queue, SIZE);
    pthread_t writers[WRITERS];
    pthread_t readers[READERS];
    for (int i = 0; i < WRITERS; ++i) {
        if (pthread_create(&writers[i], NULL, writer, (void*)&queue)) {
            fprintf(stderr, "Error creating writer thread\n");
            return 1;
        }
    }
    for (int i = 0; i < READERS; ++i) {
        if (pthread_create(&readers[i], NULL, reader, (void*)&queue)) {
            fprintf(stderr, "Error creating reader thread\n");
            return 1;
        }
    }
    for (int i = 0; i < WRITERS; ++i) {
        pthread_join(writers[i], NULL);
    }
    for (int i = 0; i < READERS; ++i) {
        pthread_join(readers[i], NULL);
    }
    delete_queue(&queue); 
    return 0;
}

int main()
{
    srand(time(NULL));   // Initialization, should only be called once.

    test(3, 1, 5, 1, 5); // {1, 1}
    putchar('\n');

    test(3, 5, 2, 1, 1); // {5, 1} 
    test(3, 5, 2, 1, 10); // {5, 1} 
    putchar('\n'); 

    test(3, 1, 1, 5, 2); // {1, 5} 
    test(3, 1, 10, 5, 2); // {1, 5} 
    putchar('\n'); 

    test(3, 5, 2, 3, 2); // {5, 3} 
    test(3, 5, 2, 2, 5); // {5, 2} 
    putchar('\n'); 

    return 0;
}