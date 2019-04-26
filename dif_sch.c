#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>

typedef struct netnode_st
{
    double t;
    double x;
    double u;
    double f;
} netnode;

MPI_Datatype netnode_to_mpi();

double func(double t, double x);
double xini(double x);
double tini(double t);

#define TSIZE 10000
#define XSIZE 1000
#define TSTEP 1.0
#define XSTEP 0.001
#define EXCH_COEF (TCOEFF * TSTEP / (XSTEP * XSTEP)) 

netnode net[TSIZE][XSIZE];
#define TCOEFF 0.0000003

void net_init();
void net_printout(FILE* file);

void lcorn(int i, int j);
void cross(int i, int j);
void rcorn(int i, int j);

MPI_Datatype NETNODE;
int np;
int rank = 0;

int main(int argc, char** argv)
{
    NETNODE = netnode_to_mpi();
    int i;
    int j;
    net_init();
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    for(i = 1; i < TSIZE; i++)
    {
        for(j = 1; j < XSIZE - 1; j++)
            cross(i, j);
        //lcorn(i, j);
    }
    net_printout(fopen("output.txt", "w"));
    return 0; 
}

double func(double t, double x)
{
    return t < 5 && x == 0.5 ? 100 : 0;
}

double xini(double x)
{
    return 0;
}

double tini(double t)
{
    return 100.0;
}

void net_init()
{
    int i, j;
    for(i = 0; i < XSIZE; i++)
        net[0][i].u = xini(i * XSTEP);
    for(i = 1; i < TSIZE; i++)
        net[i][0].u = tini(i * TSTEP);
    for(i = 0; i < TSIZE; i++)
    for(j = 0; j < XSIZE; j++)
        {
            net[i][j].x = j * XSTEP;
            net[i][j].t = i * TSTEP;
            net[i][j].f = func(net[i][j].t, net[i][j].x); 
        }
}

void net_printout(FILE* file)
{
    int i, j;
    for(i = 0; i < XSIZE; i++)
    {
        for(j = 0; j < TSIZE; j++)
            fprintf(file, "%-8lg ", net[j][i].u);
        fprintf(file, "\n");
    }
}

void lcorn(int i, int j)
{
    net[i][j].u = net[i - 1][j].u - TCOEFF * TSTEP / XSTEP * (net[i - 1][j].u - net[i - 1][j - 1].u) + net[i][j].f;
}

void rcorn(int i, int j)
{
    net[i][j].u = net[i - 1][j].u - TCOEFF * TSTEP / XSTEP * (net[i - 1][j].u - net[i - 1][j + 1].u) + net[i][j].f;
}

void cross(int i, int j)
{
   net[i][j].u = net[i - 1][j].u + EXCH_COEF * (net[i - 1][j + 1].u + net[i - 1][j - 1].u - 2 * net[i - 1][j].u) + net[i][j].f; 
}
