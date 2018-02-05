#include <stdio.h>

int main() {
    int i;
    int int_ans = 1;
    long long_ans = 1;
    float float_ans = 1;
    double double_ans = 1;

    for (i = 1; i <= 100; i++) {
        int_ans *= i;
        long_ans *= i;
        float_ans *= i;
        double_ans *= i;
    }

    printf("Int ans = %d\n", int_ans);
    printf("Long ans = %ld\n", long_ans);
    printf("Float ans = %f\n", float_ans);
    printf("Double ans = %f\n", double_ans);

    return 0;
}

