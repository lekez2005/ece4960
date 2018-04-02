#include <cstdlib>
#include <cstdio>
#include <cmath>

double monte_pi(long N){
    long long sum = 0;
    for (long i = 0; i < N; i++){
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        if ((x*x + y*y) < 1)
            sum++;
    }
    return 4*((double) sum)/N;
}

int main(){
    
    srand(1);
    printf("%12s \t %10s\n", "N", "Pi");
    for (int i = 2; i < 10; i++){
        long N = std::pow(10, i);
        double pi_est = monte_pi(N);
        printf("%12ld \t %10f\n", N, pi_est);
    }

    return 0;
}

