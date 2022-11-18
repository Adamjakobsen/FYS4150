#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, char const *argv[])
{

    #pragma omp parallel for 
    for(int i = 0; i < 10; i++){
        #pragma omp critical 
        {
            cout << omp_get_num_threads() << " " << omp_get_thread_num() << endl;
        }
    }

    return 0;
}