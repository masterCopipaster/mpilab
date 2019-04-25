#include<mpi.h>
#include<stdio.h>
int main(int argc,char **argv)
{
	int rank,size;
	char A[3][50]={"RVCE","COLLEGE","SIDDA"};
	char B[50]={},C[50]={},D[50]={};
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	if(rank==0)
	{
	printf("Rank %d started \n",rank);
	
	MPI_Recv(B,50,MPI_CHAR,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	printf("Rank %d recieve %s message \n",1,B);
	MPI_Recv(C,50,MPI_CHAR,2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	printf("Rank %d recieve %s message \n",2,C);
	MPI_Recv(D,50,MPI_CHAR,3,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	printf("Rank %d recieve %s message \n",3,D);
	

	}
	else
	{
	printf("Rank %d sends %s message \n",rank,A[rank-1]);
	MPI_Send(A[rank-1],20,MPI_CHAR,0,0,MPI_COMM_WORLD);
	
	}
MPI_Finalize();

	
}