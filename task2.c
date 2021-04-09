
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
 #include <mpi-ext.h>

#define BACKUP_ROUND  100
#define KILL_ROUND 12
int error = 0;
int kill_target = 3;
MPI_Comm main_comm, icomm;
int rank;
char file[12];
char **my_argv;

static void err_handler(MPI_Comm *pcomm, int *perr, ...) {
    error = 1;
	MPI_Comm intercomm;
	printf("Process %d handles an error with code %d\n",rank, *perr);
    MPIX_Comm_shrink(main_comm, &main_comm);
    MPI_Comm_spawn( my_argv[0], &my_argv[1], 1, MPI_INFO_NULL, 0, main_comm, &intercomm, MPI_ERRCODES_IGNORE );
    MPI_Intercomm_merge(intercomm, 0, &main_comm);

    MPI_Comm_rank(main_comm, &rank);
    sprintf(file, "%d.txt", rank);

}

float f(float x)
{
    return(1 / (1 + (x * x)));
}

FILE* exists(const char *fname)
{
    FILE *file;
    if ((file = fopen(fname, "r")))
    {
        return file;
    }
    return NULL;
}
int main(int argc, char **argv)
{
    int i,n;
    float x0 = atof(argv[1]),xn = atof(argv[2]),h = atof(argv[3]);
    int k , id, beg;
    n=(xn - x0) / h;
    my_argv = argv;
	if (n % 2 == 1)
    {
        n = n + 1;
    }
    h=(xn-x0)/n;
    double sum = 0, res = 0;
    MPI_Errhandler errh;
    int err;
 
    if ((err = MPI_Init(&argc, &argv)))  {
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    MPI_Comm_get_parent(&icomm);
    FILE *my_file;

    if (icomm != MPI_COMM_NULL) {
		MPI_Intercomm_merge(icomm, 1, &main_comm);

	} else {	
		main_comm = MPI_COMM_WORLD;
	}
	MPI_Comm_create_errhandler(err_handler, &errh);
    MPI_Comm_set_errhandler(main_comm, errh);
	MPI_Comm_size(main_comm,&k);
	MPI_Comm_rank(main_comm,&rank);
    if (icomm == MPI_COMM_NULL) {
		MPI_Barrier(main_comm);
	}
	sprintf(file, "%d.txt", rank);
    double start;
    if (rank == 0)
    { 
		start = MPI_Wtime();
	}
	int if_checked = 0;
	do {
		error = 0;
		if (my_file = exists(file)) {	
			fscanf(my_file, "%d%lf",  &beg, &sum);
			fclose (my_file);
		} else {
			beg = rank;
			sum = 0;
		}
		for (i = beg; i < n; i += k) {
			if (i % BACKUP_ROUND == 0) {
				my_file = fopen(file, "w+");
				fprintf(my_file, "%d %lf\n", i, sum);
				fclose (my_file);
			}
			if (rank == kill_target &&  icomm == MPI_COMM_NULL && i == KILL_ROUND && if_checked == 0) {
				raise(SIGKILL);
			} 
			sum += f(x0 + i * h);
	   }
	
       MPI_Allreduce(&sum, &res, 1, MPI_DOUBLE, MPI_SUM,  main_comm);
       if_checked = 1;

   } while(error);
   if (!rank) {
	    res += f(x0) / 2 + f(xn) / 2;
		printf("res:%lf\n",res * h);
		start = MPI_Wtime() - start;
		printf("Time:%lf\n",start);
	}
	MPI_Finalize();
	return 0;
}
