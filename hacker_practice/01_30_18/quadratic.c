#include <stdio.h>
#include <math.h>

//#define decimal float
#define decimal double

int main() {
    decimal a = 1e-20;
    decimal b = 1e3;
    decimal c = 1e3;

    decimal root = sqrt(pow(b, 2) - 4*a*c);
    
    decimal res1_1, res1_2, res2_1, res2_2, res3_1, res3_2;

    res1_1 = (-b + root)/(2*a);
    res1_2 = (-b - root)/(2*a);

    res2_1 = (-2*c)/(-b + root);
    res2_2 = (-2*c)/(-b - root);

    res3_1 = -c/b;
    res3_2 = -b/a + c/b;

    printf("Res 1 = %f, %f \n", res1_1, res1_2);
    printf("Res 2 = %f, %f \n", res2_1, res2_2);
    printf("Res 3 = %f, %f \n", res3_1, res3_2);

    return 0;
}

