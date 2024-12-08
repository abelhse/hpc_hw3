#import <iostream>

using std::endl;
using std::cout;
using std::swap;


double u_true(
    double x,
    double t,
    double k,
    double l,
    double u0,
    int M
) 
{
    int u = 0;
    for (int m = 0; m <= M; m++) {
        float c0 = pi * (2*m + 1) / l;
        float c1 = exp(-k * t * c0 * c0) / (2*m + 1);
        u += (4 * u0 / pi) * c1 * sin(c0 * x);
    }
    return u;
}


int main(int argc, char **argv)
{
    // h = 0.02
    double l = 1;
    int Nx = 50 + 1;
    double h = l / (Nx - 1);

    // tau = 0.0002
    double T = 0.1
    double Nt = 500 + 1;
    double tau = T / (Nt - 1);

    double k = 1;

    vector<double> u(1 + Nx + 1);
    u[0] = 0;
    u[Nx + 1] = 0;
    vector<double> u_new(1 + Nx + 1);
    u_new[0] = 0;
    u_new[Nx + 1] = 0;

    for (int j = 0; j < Nt; j++) {
        for (int i = 1; i <= Nx; i++) {
            u_new[i] = u[i] + (k * tau) / (h * h) * (u[i+1] - 2*u[i] + u[i-1]);
        }
        swap(u, u_new);
    }

    double maxAbsErr = 0;
    for (int i = 1; i <= Nx; i++) {
        double t = T;
        double x = 0 + (i-1)*h;
        u_true_i = u_true(x, t, k, l, u0, 5);
        cout << u[i] << "\t" << u_true_i << endl;
        maxAbsErr = max(maxAbsErr, abs(u[i] - u_true_i));
    }   


}