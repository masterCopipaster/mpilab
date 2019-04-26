#include <stdio.h>
#include <math.h>
//#include "/usr/include/mpich/mpi.h"
#include <mpi.h>
#include <stddef.h>

//#define PARALLEL 
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
#define XSIZE 102
#define TSTEP 0.01
#define XSTEP 0.0001
#define EXCH_COEF (TCOEFF * TSTEP / (XSTEP * XSTEP)) 

netnode net[TSIZE][XSIZE];
#define TCOEFF 0.0000003

void net_init();
void net_printout(FILE* file);

void lcorn(int i, int j);
void cross(int i, int j);
void rcorn(int i, int j);

void sequential(int i, int from, int to);
void parallel(int i);

MPI_Datatype NETNODE;
int np;
int rank = 0;

int main(int argc, char** argv)
{
    int i;
    int j;
    net_init();
    MPI_Init(&argc, &argv);
    NETNODE = netnode_to_mpi();
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    for(i = 1; i < TSIZE; i++)
    {
        #ifdef PARALLEL
        parallel(i);
        #else
        sequential(i, 1, XSIZE - 1);
        #endif
    }
    if(!rank) net_printout(fopen("output.txt", "w"));
    MPI_Finalize();
    return 0; 
}

double func(double t, double x)
{
    return t < 5 && x == 0.5 ? 100 : 0;
}

double xini(double x)
{
    return 100;
}

double tini(double t)
{
    return 0;
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

MPI_Datatype netnode_to_mpi()
{
    // Set-up the arguments for the type constructor
    MPI_Datatype new_type;

    int count = 4;
    int blocklens[] = { 1, 1, 1, 1};

    MPI_Aint indices[4];
    indices[0] = (MPI_Aint) offsetof(netnode, t);
    indices[1] = (MPI_Aint) offsetof(netnode, x);
    indices[2] = (MPI_Aint) offsetof(netnode, u);
    indices[3] = (MPI_Aint) offsetof(netnode, f);

    MPI_Datatype old_types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_struct(count,blocklens,indices,old_types,&new_type);
    MPI_Type_commit(&new_type);

    return new_type;
}

void sequential(int i, int from, int to)
{
    int j;
    for(j = from; j < to; j++)
            cross(i, j);
}

#define min(a, b) (a < b ? a : b)

void parallel(int i)
{
    int abs_start = 1 + ((XSIZE - 2) / np) * rank;
    int len = rank != np - 1 ? ((XSIZE - 2) / np) : XSIZE - 1 - abs_start;
    int abs_end = abs_start + len;
    int stutas;
    MPI_Status st;
    MPI_Request rq;

    netnode* base = &net[i][abs_start];
    int j;
    sequential(i, abs_start, abs_end);
    MPI_Ibcast(base, len, NETNODE, rank, MPI_COMM_WORLD, &rq);
    printf("barrier %d\n", i);
    for(j = 0; j < np - 1; j++)
    {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        printf("lol");
        MPI_Recv(&net[i][1 + ((XSIZE - 2) / np) * st.MPI_SOURCE], st.count_lo, NETNODE, st.MPI_SOURCE, st.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("%d received from %d/n", rank, st.MPI_SOURCE);
    }
}