#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <stdexcept>

#include <mpi.h>

using std::endl;
using std::cout;
using std::swap;
using std::vector;
using std::max;
using std::stoi;
using std::string;
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


void mainloop_sendrecv(int Nt, vector<double> &u, int myrank, int size, int Nx, vector<double> & u_new, double k, double tau, double h)
{
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
    }
}


void mainloop_isend_irecv(int Nt, vector<double> &u, int myrank, int size, int Nx, vector<double> & u_new, double k, double tau, double h)
{
    for (int j = 0; j < Nt; j++) {

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

        for (int i = 1; i <= Nx; i++) {
            u_new[i] = u[i] + (k * tau) / (h * h) * (u[i+1] - 2*u[i] + u[i-1]);
        }
        swap(u, u_new);
    }
}


void mainloop_isend_irecv_1(int Nt, vector<double> &u, int myrank, int size, int Nx, vector<double> & u_new, double k, double tau, double h)
{
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
    }
}


int main(int argc, char **argv)
{
    // инициализация MPI
    int size, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // считывание параметров
    double T = stod(argv[1]);
    int Nt = stoi(argv[2]);
    int Nx_global = stoi(argv[3]); 
    int algo = stoi(argv[4]);

    // пересылка параметров на все процессы
    MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Nx_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&algo, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // определение остальных параметров
    double u0 = 1;
    double l = 1;
    double h = l / (Nx_global - 1);
    double tau = T / (Nt - 1);
    double k = 1;

    // определение размеров блоков
    vector<int> block_sizes(size);
    int full_block_size = (Nx_global + size - 1) / size;
    for (int i = 0; i < size - 1; i++) {
        block_sizes[i] = full_block_size;
    }
    block_sizes[size-1] = (Nx_global - (size-1) * full_block_size);
    int Nx = block_sizes[myrank];
    // if (myrank == 0) {
    //     cout << "block_sizes=";
    //     for (int i = 0; i < size; i++)
    //     {
    //         cout << block_sizes[i] << " ";
    //     }
    //     cout << endl;
    // }

    // массивы значений на итерации t и t+1
    vector<double> u(1 + Nx + 1, -228);
    vector<double> u_new(1 + Nx + 1, -228);
    if (myrank == 0) {
        u[0] = 0;
        u_new[0] = 0;
    }
    if (myrank == size - 1) {
        u[u.size() - 1] = 0;
        u_new[u_new.size() - 1] = 0;
    }

    // пересылка начальных условий
    vector<int> displs(size);
    displs[0] = 0;
    for (int i = 1; i < size; i++) {
        displs[i] = displs[i-1] + block_sizes[i-1];
    }
    // if (myrank == 0 ) {
    //     cout << "displs=";
    //     for (int i = 0; i < size; i++)
    //     {
    //         cout << displs[i] << " ";
    //     }
    //     cout << endl;
    // }

    vector<double> u_init;
    if (myrank == 0) {
        u_init.resize(Nx_global);
        for (int i = 0; i < Nx_global; i++) {
            u_init[i] = 1; // тут можно задать начальное условие
        }
    }

    MPI_Scatterv(
        &u_init[0], &block_sizes[0], &displs[0], MPI_DOUBLE,
        &u[1], Nx, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // for (int i = 0; i < size; i++)
    // {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (myrank == i) {
    //         cout << "r=" << i << " u=";
    //         for (auto ui : u) {
    //             cout << ui << " ";
    //         }
    //         cout << endl;
    //     }
    // }

    // основной цикл
    double tik = MPI_Wtime();
    if (algo == 0)
        mainloop_isend_irecv_1(Nt, u, myrank, size, Nx, u_new, k, tau, h);
    else if (algo == 1)
        mainloop_isend_irecv(Nt, u, myrank, size, Nx, u_new, k, tau, h);
    else if (algo == 2)
        mainloop_sendrecv(Nt, u, myrank, size, Nx, u_new, k, tau, h);
    else {
        if (myrank == 0) {
            cout << "unknown algorithm number " << algo << endl;
        }
        MPI_Finalize();
        return -1;
    }
    double tok = MPI_Wtime();
    double totalTimeSeconds = tok - tik;

    // сравнение численных результатов с формулой
    double maxAbsErr = 0;
    for (int i = 1; i <= Nx; i++) {
        double t = T;
        double x = h * (myrank*full_block_size + (i-1));
        double u_true_i = u_true(x, t, k, l, u0, 25);
        maxAbsErr = max(maxAbsErr, fabs(u[i] - u_true_i));
    }
    double maxAbsErrReduced;
    MPI_Reduce(&maxAbsErr, &maxAbsErrReduced, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // собираем все посчитанное на главном процессе
    MPI_Gatherv(
        &u[1], Nx, MPI_DOUBLE,
        &u_init[0], &block_sizes[0], &displs[0], MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // печать результатов
    if (myrank == 0) {
        cout << "maxAbsErr=" << maxAbsErrReduced << endl;
        cout << "size=" << size << endl;
        cout << "Nx_global=" << Nx_global << endl;
        cout << "Nt=" << Nt << endl;
        cout << "T=" << T << endl;
        cout << "tau=" << tau << endl;
        cout << "h=" << h << endl;
        cout << "totalTime=" << totalTimeSeconds << endl;
        cout << "algo=" << algo << endl;

        cout << "u=";
        for (double ui : u_init) {
            cout << ui << " ";
        }
        cout << endl;

        cout << "u_true=";
        for (int i = 1; i <= Nx_global; i++) {
            double t = T;
            double x = 0 + (i-1)*h;
            double u_true_i = u_true(x, t, k, l, u0, 25);
            cout << u_true_i << " ";
        }
        cout << endl;

        cout << "x=";
         for (int i = 1; i <= Nx_global; i++) {
            double x = 0 + (i-1)*h;
            cout << x << " ";
        }
        cout << endl;
    }
    
    MPI_Finalize();
}