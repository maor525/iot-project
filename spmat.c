/*
 * spmat.c
 *
 *  Created on: 9 במאי 2020
 *      Author: itama
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "spmat.h"
#include "group.h"

int* k;
int M = 0;
int is_in(int colind, int *vertices, int n, int* g_idx);


typedef struct _spmat_array{
		int val_ind;
		double* values;
		int* colind;
		int* rowptr;
	} spmat_array;

typedef struct _node {							/*list node struct implementation*/
   int         data;
   int 			colind;
   struct _node*  next;
} node;

typedef struct linked_list{						/*list struct implementation*/
	node *head;
	node *tail;
	int len;
} list;

void add(list *row_list, node *cell){			/*add a new node to a list*/
	if(row_list -> len ==0){
		row_list -> tail = cell;
		row_list -> head = cell;
	}
	else{
		row_list -> tail -> next = cell;
		row_list -> tail = cell;
	}
	row_list -> len++;
}

int* read_i_numbers(FILE* f, int i){	/*read i values from the input file*/
	int *current_read;
	int check;
	current_read = malloc(sizeof(int)*i);
	check = fread(current_read, sizeof(int), i, f);
	if(!(check == i)){
		perror("Couldn't read from input file");
		exit(1);
	}
	return current_read;
}

spmat* build_a(FILE *input){		/*creates the sparse matrix A, allocates memory and updates values*/
	int *read_data, *neighbors;
	int degree, i, n;
	spmat *A;
	read_data = read_i_numbers(input, 1);
	n = *read_data;
	k = malloc(sizeof(int) * n);
	free(read_data);
	A = spmat_allocate_list(n);
	for (i=0; i<n; i++){
		read_data = read_i_numbers(input, 1);
		degree = *read_data;
		k[i] = degree;
		free(read_data);
		neighbors = read_i_numbers(input, degree);
		A -> add_row(A, neighbors, i, degree);
		free(neighbors);
	}
	for(i=0; i< A->n; i++) {M += *(k+i);}
	if(M==0){
		perror("Division by zero(M=0)");
		exit(1);
	}
	return A;
}

void add_row_list(struct _spmat *A, const int *row, int i, int n){			/*add row function for list implementation*/
	int j;
	node* cell;
	list* row_list;
	list** list_arr;
	list_arr = (list**) A -> private;
	row_list = *(list_arr+i);
	for(j=0; j<n ; ++j){
			cell = malloc(sizeof(node));
			cell -> colind = *(row+j);
			cell -> data = 1;
			cell -> next = NULL;
			add(row_list, cell);
	}
}

void add_row_A_g(struct _spmat *A, const double *row, int i){			/*add row function for matrix A[g]*/
	int j,n;
	node* cell;
	list* row_list;
	list** list_arr;
	n = A -> n;
	list_arr = (list**) A -> private;
	row_list = list_arr[i];
	for(j=0; j<n ; ++j){
		if(row[j] != (double) 0){
			cell = malloc(sizeof(node));
			cell -> colind = j;
			cell -> data = row[j];
			cell -> next = NULL;
			add(row_list, cell);
		}
	}
}

void free_list(struct _spmat *A){						/*frees resources used by list implementation*/
	list** list_arr;
	int i,n,len;
	node* curr;
	node* tmp;
	list_arr = A -> private;
	n = A -> n;
	for(i = 0; i<n; i++){
		len = list_arr[i] -> len;
		curr = list_arr[i] -> head;
		while(len > 0){
			tmp = curr -> next;
			free(curr);
			curr = tmp;
			len--;
		}
		free(list_arr[i]);
	}
	free(list_arr);
	free(A);
}

int get_entry(spmat *A, int x, int y){	/*returns the (x,y) entry of the matrix A*/
	list **mat;
	list *A_row;
	node *head;
	int result=0;
	mat = A -> private;
	A_row = *(mat+x);
	head = A_row -> head;
	while (head != NULL){
		if (y==head->colind){
			result = 1;
			break;
		}
		head = head -> next;
	}
	return result;
}

void get_b_g_row(spmat* A, group *g, int i, double* a_row){	/*creates the i'th row of the matrix B[g]*/
	list** list_arr;
	int n,j,len, m;
	node* cell;
	n = g -> len;
	list_arr = (list**) A -> private;
	len = list_arr[g->vertices[i]] -> len;
	cell = list_arr[g->vertices[i]] -> head;
	for(j=0; j<n; j++){a_row[j] = - (((double)k[g->vertices[i]] * (double)k[g->vertices[j]]) / (double)M);}
	for(j=0; j<len; j++){
		if(is_in(cell->colind,g->vertices, n, &m))
			a_row[m] += 1;
		cell = cell ->next;
	}
}

void get_A_row(spmat* A, group *g, int i, double* a_row){	/*creates the i'th row of the matrix A[g]*/
	list** list_arr;
	int n,j,len, k, colind, *vertices;
	node* cell;
	n = g -> len;
	list_arr = (list**) A -> private;
	vertices = g-> vertices;
	len = list_arr[g->vertices[i]] -> len;
	cell = list_arr[g->vertices[i]] -> head;
	for(j=0; j<n; j++){a_row[j] =0;}
	k=0;
	for(j=0; j<len && k<n; j++){
		colind = cell ->colind;
		if(colind > vertices[k]){
			while(colind > vertices[k] && k+1<n){k++;}
		}
		if(colind == vertices[k])
			a_row[k] = 1;
		cell = cell -> next;
	}
}


void get_D_row(group *g, int i, double* d_row){	/*creates the i'th row of the matrix D[g]*/
	int j;
	for(j=0;j<g->len;j++){
		d_row[j] = (((double)k[g->vertices[i]] * (double)k[g->vertices[j]]) / (double)M) ;
	}
}

void get_b_row(spmat* A, group *g, int i, double* b_row, int update){	/*creates the i'th row of the matrix B_hat[g]*/
	int n,j, *vertices;
	double f;
	vertices = g -> vertices;
	n = g -> len;
	f = 0;
	for(j=0; j<n; j++){
		b_row[j] = (((double)-k[vertices[i]] * (double)k[vertices[j]]) / (double)M) ;
	}
	for(j=0; j<n; j++){
		b_row[j] += (double) get_entry(A, *(vertices+i), *(vertices+j));
	}
	for(j=0; j<n; j++){f += *(b_row+j);}
	b_row[i] -= f;
	if(update){
		g->f_g[i] = f;
	}
}

double mult_vec(double* v1, double* v2, int n){		/*performs multiplication of 2 vectors*/
	double sum;
	int i;
	sum = 0;
	for(i=0; i<n; i++){
		sum += *(v1+i) * *(v2+i);
	}
	return sum;
}

int is_in(int colind, int *vertices, int n, int* g_idx){ 	/*checks if a vertex is in a group*/
	int i;
	for(i=0; i<n; i++){
		if(*(vertices+i) == colind){
			*g_idx = i;
			return 1;
		}
	}
	return 0;
}

void get_k_g(double* k_g, group *g){		/*derive k[g] out of k*/
	int i;
	for(i=0; i<g->len; i++){
		k_g[i] = k[g->vertices[i]];
	}
}

void mult_A_g(double *x, double *result, group *g, double norm1, double * k_g){			/*multiply sparse matrix by vec using list implementation*/
	list** list_arr;
	int i,n,j,len;
	double inner_prod;
	node* cell;
	list_arr = (list**) g->A_g -> private;
	n = g -> len;
	inner_prod = (double) mult_vec(k_g, x, n) / (double) M;
	for(i=0; i<n; i++){
		result[i] = (double) 0;
		len = list_arr[i] -> len;
		cell = list_arr[i] -> head;
		for(j=0; j<len; j++){
			result[i] += (cell -> data) * x[cell -> colind];
			cell = cell ->next;
		}
		result[i] +=  ( norm1 - g->f_g[i] )* *(x+i) - inner_prod * *(k_g+i);
	}
}

void temp(struct _spmat *A, double *x, double *result, group *g, double norm1, double * k_g, double *a_row){			/*multiply B_hat[g] + ||B_hat[g]||I  by vec using list implementation*/
	int i,n;
	double inner_prod;
	n = g -> len;
	inner_prod = (double) mult_vec(k_g, x, n) / (double) M;
	for(i=0; i<n; i++){
		result[i] = (double) 0;
		get_A_row(A, g, i, a_row);
		result[i] += mult_vec(a_row,x,n) + ( norm1 - g->f_g[i] )* *(x+i) - inner_prod * *(k_g+i);
	}
}


double get_norm(spmat* A, group *g){	/* calculate the L1-norm of B[g]*/
	double *b_row, max, *sum;
	int n, i, j;
	n = g ->len;
	b_row = malloc(sizeof(double) * n);
	sum = malloc(sizeof(double) * n);
	for(j=0; j<n; j++){
				sum[j] = 0;
		}
	for(i=0; i<n; i++){
		get_b_row(A, g, i, b_row, 1);
		for(j=0; j<n; j++){
			sum[j] += fabs(*(b_row+j));
		}
	}
	max = *sum;
	for(i=1; i<n; i++){
		if(*(sum+i) > max)
			max = *(sum+i);
	}
	free(b_row);
	free(sum);
	return max;
}

void mult_list(struct _spmat *A, double *v, double *result, group *g, double norm1){			/*multiply sparse matrix by vec using list implementation*/
	int i,n;
	double *b_row;
	n = g -> len;;
	b_row = malloc(sizeof(double) * n);
	for(i=0; i<n; i++){
		get_b_row(A, g, i, b_row, 0);
		b_row[i] += norm1;
		result[i] = mult_vec(b_row, v, n);
	}
	free(b_row);
}

void mult_b_hat(struct _spmat *A, double *v, double *result, group *g){			/*multiply b_hat by a vector*/
	int i,n;
	double *b_row;
	n = g -> len;
	b_row = malloc(sizeof(double) * n);
	for(i=0; i<n; i++){
		get_b_row(A, g, i, b_row, 0);
		result[i] = mult_vec(b_row, v, n);
	}
	free(b_row);
}

void norm_vec(double* v, int n){
	int i;
	double sum;
	sum = 0;
	for(i=0; i<n; i++){
		sum += *(v+i)**(v+i);
	}
	sum = sqrt(sum);
	for(i=0; i<n; i++){
		v[i] = *(v+i) / sum;
	}
}

spmat* spmat_allocate_list(int n){				/*allocate resources for a sparse matrix using list implementation*/
	int i;
	spmat* spmat_org;
	list** spmat_list;
	spmat_list = malloc(sizeof(list*) * n);
	for(i=0; i<n; ++i){
		spmat_list[i] = malloc(sizeof(list));
		spmat_list[i] -> head = NULL;
		spmat_list[i] -> tail = NULL;
		spmat_list[i] -> len = 0;
	}
	spmat_org = malloc(sizeof(spmat));;
	spmat_org -> n = n;
	spmat_org -> private = spmat_list;
	spmat_org -> add_row = add_row_list;
	spmat_org -> free = free_list;
	spmat_org -> mult = mult_list;

	return spmat_org;
}

void get_A_g(spmat* A, group* g){
	double *a_row;
	int i, n;
	n =g->len;
	a_row = malloc(sizeof(double) * n);
	g->A_g = spmat_allocate_list(n);
	for(i=0;i<n;i++){
		get_A_row(A,g,i,a_row);
		add_row_A_g(g->A_g,a_row,i);
	}
	free(a_row);
}

void power_iter(spmat* A, group *g, double norm1, double *ret_result){		/*power iteration function*/
	int i, n,cnt, limit;
	double *v,*tmp, *k_g, *result;
	bool advance = true;
	double epsilon = 0.00001;
	n = g-> len;
	k_g = malloc(sizeof(double) * n);
	v = malloc(sizeof(double) * n);
	result = malloc(sizeof(double)* n);
	get_k_g(k_g, g);
	get_A_g(A,g);
	srand(time(NULL));						/*initial vector is randomized and loaded into v*/
	for (i = 0; i < n; ++i){
		v[i] =  rand()% 100;
	}
	cnt = 0;
	limit = ((0.5)*(n)*(n) + 10000*(n) + 300000);
	while(advance){							/*power iteration loop using mult*/
		if(cnt > limit){
			perror("Power Iteration exceeded iteration limit");
			exit(1);
		}
		advance = false;
		mult_A_g(v,result,g,norm1, k_g);
		norm_vec(result, n);
		for(i = 0; i < n; i++){
			if(fabs(*(v+i) - *(result+i)) >= epsilon){
				tmp = result;
				result = v;
				v = tmp;
				advance = true;
				break;
			}
		}
		cnt++;
	}
	for (i=0; i<n; i++) *(ret_result+i) = *(result+i);
	free(result);
	free(k_g);
	free(v);	/*free resources used*/
}

double get_eigen_value(spmat* A, group *g, double* e_vec, double norm1){
	double e_val,*result;
	result = malloc(sizeof(double) * g->len);
	mult_list(A,e_vec,result, g, norm1);
	e_val = mult_vec(e_vec, result, g -> len);
	e_val = e_val / mult_vec(e_vec,e_vec, g -> len);
	free(result);
	return e_val;
}

group_node** split_by_s(group_node** split, group* g, double *s){	/*splits a group g into 2 groups according to the vector s*/
	group_node *g1, *g2;
	int i, n, g1_len, g2_len, g1_i, g2_i;
	g1_len = 0; g2_len = 0; g1_i = 0; g2_i = 0;
	n = g -> len;
	g1 = malloc(sizeof(group_node));
	g2 = malloc(sizeof(group_node));
	g1 -> data = malloc(sizeof(group));
	g2 -> data = malloc(sizeof(group));
	for(i = 0; i < n; i++){
		if(*(s+i) > 0)
			g1_len++;
		else
			g2_len++;
	}
	g1 -> data -> vertices = malloc(sizeof(int) * g1_len);
	g2 -> data -> vertices = malloc(sizeof(int) * g2_len);
	for(i = 0; i < n; i++){
		if(s[i] > 0){
			g1 -> data -> vertices[g1_i] = g -> vertices[i];
			g1_i++;
		}
		else{
			g2 -> data -> vertices[g2_i] = g -> vertices[i];
			g2_i++;
		}
	}
	g1->data->len = g1_len;
	g2->data->len = g2_len;
	g1->data->f_g = malloc(sizeof(double) * g1_len);
	g2->data->f_g = malloc(sizeof(double) * g2_len);
	*split = g1; *(split+1) = g2;
	return split;
}

void build_s(double* e_vector, int n, double *s){	/*derive vector s out of the eigen vector*/
	int i;
	for (i=0; i<n; i++){
		if (*(e_vector+i) >=0){
			*(s+i) = 1;
		}
		else{
			*(s+i) = -1;
		}
	}
}

double calc_score(double *s,double *x, int i, group *g){ /*calculates the initial i'th entry of the score vec*/
	return (-2 * (s[i] * x[i] + (((double)k[g->vertices[i]] * (double)k[g->vertices[i]]) / (double)M)) );
}

double update_score(double *s, int i, int l, group *g, spmat *A){ /*calculates the value to update the i'th entry of the score vec*/
	double entry;
	entry = (double) get_entry(A, *(g->vertices+i), *(g->vertices+l)) - (((double)k[g->vertices[i]] * (double)k[g->vertices[l]]) / (double)M);
	return (entry * s[i] *s[l] *4);
}

double* mult(spmat* A, group *g, double *s){ /*multiply A by vec s*/
	double *row, *res;
	int i,n;
	n = g -> len;
	row = malloc(sizeof(double) * n);
	res = malloc(sizeof(double) * n);
	for(i = 0;i<n; i++){
		get_b_g_row(A,g,i,row);
		res[i] = mult_vec(row,s,n);
	}
	free(row);
	return res;
}

void free_k(){
	free(k);
}
