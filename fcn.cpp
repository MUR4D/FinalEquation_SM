#include<iostream>
#include<cmath>
#define N 2

void fcn(double t, double z[N], double f[N])
{

    f[0] = z[1];
    f[1] = -sinh(t) * z[1] - t * z[0] - t;

}
