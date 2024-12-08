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

    double tik = MPI_Wtime();

    double u0 = 1;

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

    for (int j = 0; j < Nt; j++) {

        // делаем асинхронный запрос на обмен концами
        MPI_Request req[4];
        // >
        MPI_Isend(
            &u[u.size() - 2], 1, MPI_DOUBLE, (myrank+1 + size)%size, 0,
            MPI_COMM_WORLD, &req[0]
        );
        MPI_Irecv(
            &u[0], 1, MPI_DOUBLE, (myrank-1 + size)%size, 0,
            MPI_COMM_WORLD, &req[1]
        );
        // <
        MPI_Isend(
            &u[1], 1, MPI_DOUBLE, (myrank-1 + size)%size, 1,
            MPI_COMM_WORLD, &req[2]
        );
        MPI_Irecv(
            &u[u.size() - 1], 1, MPI_DOUBLE, (myrank+1 + size)%size, 1,
            MPI_COMM_WORLD, &req[3]
        );

        // пока асинхронно идет обмен концами, считаем то, что от них не зависит
        for (int i = 1; i <= Nx; i++) {
            if ((i == 1) || (i == Nx)) {
                continue;
            }
            u_new[i] = u[i] + (k * tau) / (h * h) * (u[i+1] - 2*u[i] + u[i-1]);
        }

        // завершаем обмен концами
        MPI_Status statuses[4];
        MPI_Waitall(4, req, statuses);
        if (myrank == 0) {
            u[0] = 0;
            u_new[0] = 0;
        }
        if (myrank == size - 1) {
            u[u.size() - 1] = 0;
            u_new[u_new.size() - 1] = 0;
        }

        // обновляем то, что зависит от концов
        for (int i : {1, Nx}) {
            u_new[i] = u[i] + (k * tau) / (h * h) * (u[i+1] - 2*u[i] + u[i-1]);
        }

        swap(u, u_new);

        // сравнение Ug и Utrue (соответствует ли многопроцессная версия формуле)
        {
            double maxAbsErr = 0;
            for (int i = 1; i <= Nx; i++) {
                double t = j * tau;
                double x = h * (myrank*chuck_size + (i-1));
                double u_true_i = u_true(x, t, k, l, u0, 25);
                maxAbsErr = max(maxAbsErr, fabs(u[i] - u_true_i));
            }
            double maxAbsErrReduced;
            MPI_Reduce(&maxAbsErr, &maxAbsErrReduced, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myrank == 0) {
                if (j % (Nt / 10) == 0) {
                    cout << "j=" << j;
                    cout << " max|U-Utrue|=" << maxAbsErrReduced;
                    cout << endl;
                }
            }
        }

    }

    double tok = MPI_Wtime();
    double totalTimeSeconds = tok - tik;
    if (myrank == 0) {
        cout << "totalTimeSeconds=" << totalTimeSeconds << endl;
    }
    
    MPI_Finalize();
}