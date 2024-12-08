#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cassert>

#include <mpi.h>

using std::endl;
using std::cout;
using std::swap;
using std::vector;
using std::max;
using std::stoi;
using std::stod;
using std::ofstream;


double u_true(
    double x,
    double t,
    double k,
    double l,
    double u0,
    int M
) 
{
    double u = 0;
    for (int m = 0; m <= M; m++) {
        float c0 = M_PI * (2*m + 1) / l;
        float c1 = exp(-k * t * c0 * c0) / (2*m + 1);
        u += (4 * u0 / M_PI) * c1 * sin(c0 * x);
    }
    return u;
}


int main(int argc, char **argv)
{
    int size, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    double u0 = 1;

    // // h = 0.02
    // double l = 1;
    // int Nx_global = 50 + 1;
    // double h = l / (Nx_global - 1);

    // // tau = 0.0002
    // double T = 0.1;
    // double Nt = 500 + 1;
    // double tau = T / (Nt - 1);

    double l = 1;
    int Nx_global = stoi(argv[3]); 
    double h = l / (Nx_global - 1);

    double T = stod(argv[1]);
    int Nt = stoi(argv[2]);
    double tau = T / (Nt - 1);

    double k = 1;

    // размеры блоков
    int chuck_size = (Nx_global + size - 1) / size;
    int Nx = chuck_size;
    if (myrank == size-1) { // последний кусочек может быть меньше
        Nx = Nx_global - (size-1) * chuck_size;
    }
    cout << "rank=" << myrank << " Nx=" << Nx << endl; // TODO barrier

    // проверим, что правильно рассчитали размеры блоков
    int Nx_global_check;
    MPI_Reduce(&Nx, &Nx_global_check, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
        assert(Nx_global_check == Nx_global);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    vector<double> u(1 + Nx + 1, 1);
    vector<double> u_new(1 + Nx + 1, 1);
    if (myrank == 0) {
        u[0] = 0;
        u_new[0] = 0;
    }

    if (myrank == size - 1) {
        u[u.size() - 1] = 0;
        u_new[u_new.size() - 1] = 0;
    }
    
    // 1-процессная версия
    vector<double> ug(1 + Nx_global + 1, 1);
    vector<double> ug_new(1 + Nx_global + 1, 1);
    ug[0] = ug[ug.size() - 1] = 0;
    ug_new[0] = ug_new[ug_new.size() - 1] = 0;


    for (int j = 0; j < Nt; j++) {
        MPI_Status status;

        // >
        MPI_Sendrecv(
            &u[u.size() - 2], 1, MPI_DOUBLE, (myrank+1 + size)%size, 0,
            &u[0], 1, MPI_DOUBLE, (myrank-1 + size)%size, 0,
            MPI_COMM_WORLD, &status
        );

        // <
        MPI_Sendrecv(
            &u[1], 1, MPI_DOUBLE, (myrank-1 + size)%size, 1,
            &u[u.size() - 1], 1, MPI_DOUBLE, (myrank+1 + size)%size, 1,
            MPI_COMM_WORLD, &status
        );

        if (myrank == 0) {
            u[0] = 0;
            u_new[0] = 0;
        }

        if (myrank == size - 1) {
            u[u.size() - 1] = 0;
            u_new[u_new.size() - 1] = 0;
        }

        for (int i = 1; i <= Nx; i++) {
            u_new[i] = u[i] + (k * tau) / (h * h) * (u[i+1] - 2*u[i] + u[i-1]);
        }
        swap(u, u_new);

        // сравнение Ug и Utrue (соответствует ли многопроцессная версия формуле)
        {
            double maxAbsErrUandUtrue = 0;
            for (int i = 1; i <= Nx; i++) {
                // double t = T;
                double t = j * tau;
                // double x = 0 + (i-1)*h;
                double x = h * (myrank*chuck_size + (i-1));
                double u_true_i = u_true(x, t, k, l, u0, 25);
                // of_u_true << u_true_i << " ";
                maxAbsErrUandUtrue = max(maxAbsErrUandUtrue, fabs(u[i] - u_true_i));
            }
            double maxAbsErrUandUtrueReduced;
            MPI_Reduce(&maxAbsErrUandUtrue, &maxAbsErrUandUtrueReduced, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myrank == 0) {
                if (j % (Nt / 10) == 0) {
                    cout << "j=" << j;
                    cout << " max|U-Utrue|=" << maxAbsErrUandUtrueReduced;
                    cout << endl;
                }
            }
        }

        // расчет 1-проц версии
        for (int i = 1; i <= Nx_global; i++) {
            ug_new[i] = ug[i] + (k * tau) / (h * h) * (ug[i+1] - 2*ug[i] + ug[i-1]);
        }
        swap(ug, ug_new);

        // сравнение Ug и U (соответствует ли 1-процессная версия многопроцессной)
        {
            double maxAbsErrUandUg = 0;
            for (int i = 1; i <= Nx; i++) {
                int ig = 1 + myrank*chuck_size + (i-1);
                maxAbsErrUandUg = max(maxAbsErrUandUg, fabs(u[i] - ug[ig]));
            }
            double maxAbsErrUandUgReduced;
            MPI_Reduce(&maxAbsErrUandUg, &maxAbsErrUandUgReduced, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myrank == 0) {
                if (j % (Nt / 10) == 0) {
                    cout << "j=" << j;
                    cout << " max|U-Ug|=" << maxAbsErrUandUgReduced;
                    cout << endl;
                }
            }
        }

        // сравнение Ug и Utrue (соответствует ли 1-процессная версия формуле)
        // должно быть так же как в (U и Utrue)
        {
            double maxAbsErrUgandUtrue = 0;
            for (int i = 1; i <= Nx; i++) {
                double t = j * tau;
                double x = h * (myrank*chuck_size + (i-1));
                double u_true_i = u_true(x, t, k, l, u0, 25);

                int ig = 1 + myrank*chuck_size + (i-1);

                maxAbsErrUgandUtrue = max(maxAbsErrUgandUtrue, fabs(u_true_i - ug[ig]));
            }
            double maxAbsErrUgandUtrueReduced;
            MPI_Reduce(&maxAbsErrUgandUtrue, &maxAbsErrUgandUtrueReduced, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myrank == 0) {
                if (j % (Nt / 10) == 0) {
                    cout << "j=" << j;
                    cout << " max|Ug-Utrue|=" << maxAbsErrUgandUtrueReduced;
                    cout << endl << endl;
                }
            }
        }
    }
    
    MPI_Finalize();
}