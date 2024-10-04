#include <iostream>
#include <vector>
#include <limits>
#include "papi.h"
#include <unordered_map>

class CSR_graph {
    int row_count; //number of vertices in graph
    unsigned int col_count; //number of edges in graph
    
    std::vector<unsigned int> row_ptr;
    std::vector<int> col_ids;
    std::vector<double> vals;
    
    std::unordered_map<int, double> calculated_weights;
public:

    void read(const char* filename) {
        FILE *graph_file = fopen(filename, "rb");
        fread(reinterpret_cast<char*>(&row_count), sizeof(int), 1, graph_file);
        fread(reinterpret_cast<char*>(&col_count), sizeof(unsigned int), 1, graph_file);

        std::cout << "Row_count = " << row_count << ", col_count = " << col_count << std::endl;
        
        row_ptr.resize(row_count + 1);
        col_ids.resize(col_count);
        vals.resize(col_count);
         
        fread(reinterpret_cast<char*>(row_ptr.data()), sizeof(unsigned int), row_count + 1, graph_file);
        fread(reinterpret_cast<char*>(col_ids.data()), sizeof(int), col_count, graph_file);
        fread(reinterpret_cast<char*>(vals.data()), sizeof(double), col_count, graph_file);
        fclose(graph_file);
    }

    void print_vertex(int idx) {
        for (int col = row_ptr[idx]; col < row_ptr[idx + 1]; col++) {
            std::cout << col_ids[col] << " " << vals[col] <<std::endl;
        }
        std::cout << std::endl;
    }

    void reset() {
        row_count = 0;
        col_count = 0;
        row_ptr.clear();
        col_ids.clear();
        vals.clear();
    }

    void find_vertex_max()
    {
        double max_sum = -std::numeric_limits<double>::min();
        int max_vertex = INT_MIN;
        for (int i = 0;i < row_count;++i)
        {
            double curr_sum = 0.0;
            for (int j = row_ptr[i];j < row_ptr[i + 1];++j)
            {
                if (col_ids[j] % 2 == 0)
                    curr_sum += vals[j];
            }
            if (curr_sum > max_sum)
            {
                max_sum = curr_sum;
                max_vertex = i;
            }
        }

        std::cout << "Vertex with the largest total weight of incident edges " << max_vertex << ", total weight: " << max_sum << std::endl;
    }

    void find_max_rank_vertex()
    {
        double max_rank = -std::numeric_limits<double>::min();
        int max_vertex = INT_MIN;

        for (int i = 0;i < row_count;++i) 
        {
            double curr_rank = 0.0;
            for (int j = row_ptr[i];j < row_ptr[i + 1];++j) 
            {
                int neighbor = col_ids[j];
                double w_edge = vals[j];

                double neighbor_weight = 0.0;
                for (int k = row_ptr[neighbor];k < row_ptr[neighbor + 1];++k) 
                {
                    int n = col_ids[k];
                    neighbor_weight += vals[k] * (row_ptr[n + 1] - row_ptr[n]); 
                }
                curr_rank += w_edge * neighbor_weight;
            }
            if (curr_rank > max_rank) {
                max_rank = curr_rank;
                max_vertex = i;
            }
        }  
        std::cout << "Vertex with max rank " << max_vertex << ", Rank: " << max_rank << std::endl;
    }
};


#define N_TESTS 5

int main () {
    const char* filenames[N_TESTS];
    filenames[0] = "synt";
    filenames[1] = "road_graph";
    filenames[2] = "stanford";
    filenames[3] = "youtube";
    filenames[4] = "syn_rmat";

    /* https://drive.google.com/file/d/183OMIj56zhqN12Aui1kxv76_Ia0vTPIF/view?usp=sharing архив с тестами, 
        распаковать командой tar -xzf 
    */

    for (int n_test = 0; n_test < N_TESTS; n_test++) {
        CSR_graph a;
        a.read(filenames[n_test]);
        
        PAPI_library_init(PAPI_VER_CURRENT);

        int event = PAPI_NULL;
        PAPI_create_eventset(&event);
        int event_code;
        int flag1 = PAPI_event_name_to_code("PAPI_L1_TCM", &event_code);
        if (!flag1)
            PAPI_add_event(event, event_code); 
        
        int flag2 = PAPI_event_name_to_code("PAPI_L2_TCM", &event_code);
        if (!flag2)
            PAPI_add_event(event, event_code); 

        int flag3 = PAPI_event_name_to_code("PERF_COUNT_HW_CPU_CYCLES", &event_code);
        if (!flag3)
            PAPI_add_event(event, event_code); 

        long long values[3];

        PAPI_start(event); 
        a.find_vertex_max();
        PAPI_stop(event, values); 
        if (!flag1)
            std::cout << "Test " << n_test << " algo1 L1 cache misses: " << values[0] << std::endl;
        if (!flag2)
            std::cout << "Test " << n_test << " algo1 L2 cache misses: " << values[1] << std::endl;
        if (!flag3)
            std::cout << "Test " << n_test << " algo1 PERF_COUNT_HW_CPU_CYCLES " << values[2] << std::endl;

        values[0] = 0;
        values[1] = 0;
        values[2] = 0;
        PAPI_start(event); 
        a.find_max_rank_vertex();
        PAPI_stop(event, values); 
        if (!flag1)
             std::cout << "Test " << n_test << " algo2 L1 cache misses: " << values[0] << std::endl;
        if (!flag2)
            std::cout << "Test " << n_test << " algo2 L2 cache misses: " << values[1] << std::endl;
        if (!flag3)
            std::cout << "Test " << n_test << " algo2 PERF_COUNT_HW_CPU_CYCLES " << values[2] << std::endl;

        PAPI_cleanup_eventset(event);
        PAPI_destroy_eventset(&event);
        PAPI_destroy_eventset(&event);
        PAPI_shutdown();
    }
}