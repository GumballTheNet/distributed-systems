
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 4


int best_val;
int min_indeces[2];
int rank, tasks;
MPI_Comm comm;
MPI_Status status;

void spread(int *coords, int num) {
	int offset = (1 + coords[num]) <= SIZE / 2 ? 1 : -1, is_border = (1 + coords[num]) == SIZE / 2 || (1 + coords[num]) == SIZE / 2 + 1;
	int prev_other[2] = {coords[0] - offset * !num, coords[1] - offset * num}, other_num = 0,post_other[2] = {coords[0] + offset* !num, coords[1] + offset * num} ;
	int his_res = 0;
	if (coords[num] && coords[num] != SIZE - 1) {
		MPI_Cart_rank(comm, prev_other, &other_num);
		MPI_Recv(&his_res, 1, MPI_INT, other_num, 0, comm, &status);
		int best_other[2];
		MPI_Recv(best_other, 2, MPI_INT, other_num, 0, comm, &status);
		if (his_res < best_val || ((his_res == best_val) && (best_other[0] < min_indeces[0] ||
		(best_other[0] == min_indeces[0] && best_other[1] < min_indeces[1])))) {
			best_val = his_res;
			min_indeces[0] = best_other[0];
			min_indeces[1] = best_other[1];
		}
	}
	if (!is_border) {
		MPI_Cart_rank(comm, post_other, &other_num);

		MPI_Send(&best_val, 1, MPI_INT, other_num, 0, comm);
		MPI_Send(min_indeces, 2, MPI_INT, other_num, 0, comm);
	}
	 
	MPI_Barrier(comm);
	if (is_border) {
		if ((1 + coords[num]) == SIZE / 2) {
			int other_num, other[2] =  {coords[0] + !num, coords[1] + num} ;
			MPI_Cart_rank(comm, other, &other_num);
			MPI_Send(&best_val, 1, MPI_INT, other_num, 0, comm);
			MPI_Send(min_indeces, 2, MPI_INT, other_num, 0, comm);
			int his_res = 0;
			MPI_Recv(&his_res, 1, MPI_INT, other_num, 0, comm, &status);
			int best_other[2];
			MPI_Recv(best_other, 2, MPI_INT, other_num, 0, comm, &status);

			if (his_res < best_val || ((his_res == best_val) && (best_other[0] < min_indeces[0] ||
			(best_other[0] == min_indeces[0] && best_other[1] < min_indeces[1])))) {
				best_val = his_res;
				min_indeces[0] = best_other[0];
				min_indeces[1] = best_other[1];
			}
		} else {
			int other_num, other[2] =  {coords[0] - !num, coords[1] -  num} ;
			MPI_Cart_rank(comm, other, &other_num);
			
			int his_res = 0;
			MPI_Recv(&his_res, 1, MPI_INT, other_num, 0, comm, &status);
			int best_other[2];
			MPI_Recv(best_other, 2, MPI_INT, other_num, 0, comm, &status);

			if (his_res < best_val || ((his_res == best_val) && (best_other[0] < min_indeces[0] ||
			(best_other[0] == min_indeces[0] && best_other[1] < min_indeces[1])))) {
				best_val = his_res;
				min_indeces[0] = best_other[0];
				min_indeces[1] = best_other[1];
			}
			MPI_Send(&best_val, 1, MPI_INT, other_num, 0, comm);
			MPI_Send(min_indeces, 2, MPI_INT, other_num, 0, comm);
		}
	}
	MPI_Barrier(comm);
	prev_other[0] = coords[0] + offset * !num;
	prev_other[1] = coords[1] + offset * num;
	other_num = 0;
	post_other[0] = coords[0] - offset * !num;
	post_other[1] = coords[1] - offset * num ;
	if (!is_border) {
		his_res = 0;
		MPI_Cart_rank(comm, prev_other, &other_num);
		MPI_Recv(&his_res, 1, MPI_INT, other_num, 0, comm, &status);
		int best_other[2];
		MPI_Recv(best_other, 2, MPI_INT, other_num, 0, comm, &status);

		if (his_res < best_val || ((his_res == best_val) && (best_other[0] < min_indeces[0] ||
		(best_other[0] == min_indeces[0] && best_other[1] < min_indeces[1])))) {
			best_val = his_res;
			min_indeces[0] = best_other[0];
			min_indeces[1] = best_other[1];
		}
	}
	if (coords[num] && coords[num] != SIZE - 1) {
		MPI_Cart_rank(comm, post_other, &other_num);
		MPI_Send(&best_val, 1, MPI_INT, other_num, 0, comm);
		MPI_Send(min_indeces, 2, MPI_INT, other_num, 0, comm);
	}
	
	MPI_Barrier(comm);
}

int main(int argc, char **argv)
{
	int i, best_num;
	MPI_Init(&argc, &argv);
  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);
    int size[2] = {SIZE, SIZE};
    srand(0);

    best_val = rand() % 1000;
    int is_periodic[2] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, is_periodic, 0, &comm);
    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    printf("Coordinates for process %d: (%d, %d)\n", rank, coords[0], coords[1]);
    min_indeces[0] = coords[0];
    min_indeces[1] = coords[1];
    printf("a[%d][%d] = %d\n", coords[0], coords[1],best_val);
    MPI_Barrier(comm);

	spread(coords, 0);
	spread(coords, 1);
	MPI_Cart_rank(comm, min_indeces, &best_num);

	printf("my rank: %d, my res: %d, first best proc: %d\n",rank, best_val, best_num);
	MPI_Finalize();

	return 0;
}
