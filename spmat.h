#ifndef _SPMAT_H
#define _SPMAT_H
#include "group.h"

/*This module is used to handle sparse matrices and store them in memory.
 * the functions in this module are used to calculate and derive data
 * from sparse matrices in order to maximize modularity*/

typedef struct _spmat {
	/* Matrix size (n*n) */
	int		n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void	(*add_row)(struct _spmat *A, const int *row, int i, int n);

	/* Frees all resources used by A */
	void	(*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void	(*mult)(struct _spmat *A, double *v, double *result, group *g, double norm1);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void	*private;
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate_list(int n);

#endif
