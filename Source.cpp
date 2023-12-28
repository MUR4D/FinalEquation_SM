#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include<iostream>
#include<string>
#include "dopri5.h"
#define delta 1.e-6
#define N 2
#define T 1.0
using namespace std;



int main(void) {

    int i;
    double z[N], tend = 1.0, eps, x1, y1, glob1;
    eps = 1.e-7;

    double t0 = 0, x0 = 0, h, tn = 1.0, xn = 0, z0, m1, m2, m3, b, b1, b2, b3, e;
    b = xn;
    z[1] = -1;
    z[0] = -4.766669009;//-4.7666987;

    b1 = dopri5(z, t0, T, eps);
    printf("\nB1 is %.8lf", b1);

    return 0;
}
