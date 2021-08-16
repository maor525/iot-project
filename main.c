/*
 * main.c
 *
 *  Created on: Aug 22, 2020
 *      Author: Maor
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include "spmat.h"
#include "group.h"
#define epsilon 0.00001


void Algorithm4(double *s, group *g, spmat *A);
void power_iter(spmat* A, group *g, double norm1, double *result);
double get_eigen_value(spmat* A, group *g, double* e_vec, double norm1);
void build_s(double* e_vector, int n, double *s);
void mult_b_hat(struct _spmat *A, double *v, double *result, group *g);
double mult_vec(double* v1, double* v2, int n);
group_node** split_by_s(group_node** split, group* g, double *s);
double get_norm(spmat* A, group *g);
spmat* build_a(FILE *input);
void free_group_list(group_list *g);
double calc_score(double *s,double *x, int i, group *g);
double* mult(spmat* A, group *g, double *s);
double update_score(double *s, int i, int l, group *g, spmat *A);
void free_k();


group_node** Algorithm2(spmat* A, group *g){
	double *e_vec, e_val, *tmp, Q, *s, norm1;
	group_node **split;


	tmp = malloc(sizeof(double) * g-> len);
	split = malloc(sizeof(group_node*) * 2);
	norm1 = get_norm(A, g);
	e_vec = malloc(sizeof(double) * g -> len);
	power_iter(A,g, norm1, e_vec);
	e_val = get_eigen_value(A, g, e_vec, norm1) - norm1;

	if(e_val <= epsilon){
		free(split);
		free(e_vec);
		free(tmp);
		return NULL;
	}
	s = malloc(sizeof(double)* g-> len);
	build_s(e_vec, g -> len, s);
	mult_b_hat(A, s, tmp, g);
	Q = mult_vec(s, tmp, g -> len);
	if (Q <= epsilon){
		free(split);
		free(e_vec);
		free(s);
		free(tmp);
		return NULL;
	}
	else{
		Algorithm4(s,g,A);
		split_by_s(split, g, s);

	}
	free(e_vec);
	free(s);
	free(tmp);
	return split;
}

group_list* Algorithm3(spmat* A){
	group_list *O, *P;
	int i, cnt;
	group *g; group_node *g_n; group_node *curr;
	group_node** split;


	P = malloc(sizeof(group_list));
	O = malloc(sizeof(group_list));
	g = (group*) malloc(sizeof(group));
	g_n = malloc(sizeof(group_node));
	g->vertices = malloc(sizeof(int) * A->n);
	g->f_g = malloc(sizeof(double) * A->n);
	g->len = A->n;
	for(i=0; i< A->n; i++){
		g->vertices[i] = i;
	}
	P-> len = 0; O->len = 0;
	g_n ->data = g;
	add_group(P,g_n);
	cnt = 0;
	while(P->len != 0){
		cnt++;
		curr = pop_group(P);
		curr->next = NULL;

		split = Algorithm2(A, curr->data);
		curr->data->A_g->free(curr->data->A_g);
		if( (split == NULL)){
			add_group(O, curr);
		}
		else{
			free_group_node(curr);
			if(split[0] ->data-> len == 1){
				add_group(O, *split);
			}
			else{
				add_group(P, *split);
			}
			if(split[1] ->data-> len == 1){
				add_group(O, *(split+1));
			}
			else{
				add_group(P, *(split+1));
			}
		}
		free(split);

	}
	free_group_list(P);
	return O;
}

void Algorithm4(double *s, group *g, spmat *A){
	double delta_q, *score, k, *improve, *x;
	int n, i, i_tag, j, j_tag, *unmoved, *indices, unm_len, cnt=0, tmp, l;
	n = g->len;
	improve = malloc(sizeof(double) * n);
	score = malloc(sizeof(double) * n);
	indices = malloc(sizeof(int) * n);
	unmoved = malloc(sizeof(int) * n);
	do{
		unm_len = n;
		for(i=0; i<n ; i++){
			*(unmoved+i) = i;
			improve[i] = 0;
		}
		for(i=0; i<n ; i++){
			if(i == 0){
			x = mult(A,g,s);
			for(l=0;l<n;l++){
				score[l] = calc_score(s,x, l, g);
			}
			free(x);
			}
			else{
				for(l=0;l<n;l++){
					if(l == indices[i-1])
						score[l] = -score[l];
					else
						score[l] = score[l] - update_score(s, l, indices[i-1], g, A);
				}
			}

			j_tag = unmoved[0];
			k = score[unmoved[0]];
			for(j=1; j< unm_len; j++){
				tmp = *(unmoved+j);
				if(k < score[tmp]){
					k = score[tmp];
					j_tag = tmp;
				}
			}
			s[j_tag] = -*(s+j_tag);
			indices[i] = j_tag;
			if(i == 0)
				improve[i] = *(score+j_tag);
			else
				improve[i] = *(improve+i-1) + *(score+j_tag);
			del_val(unmoved, j_tag, unm_len);
			unm_len--;
		}
		i_tag = arg_max(improve, n);
		for(i = n-1; i> i_tag; i--){
			j = *(indices+i);
			s[j] = -*(s+j);
		}
		if(i_tag == n-1)
			delta_q = 0;
		else
			delta_q = *(improve+i_tag);
		cnt++;
	}while(delta_q > epsilon);
	free(unmoved);
	free(improve);
	free(score);
	free(indices);
}

void write_to_out(FILE *out, group_list *O){
	int n, check;
	group_node *g;
	n = O->len;
	check = fwrite(&n, sizeof(int), 1, out);
	if(!(check == 1)){
		perror("Couldn't write to output file");
		exit(1);
	}
	g = O->head;
	while (g!=NULL){
		n = g ->data->len;
		check = fwrite(&n, sizeof(int), 1, out);
		if(!(check == 1)){
				perror("Couldn't write to output file");
				exit(1);
			}
		check = fwrite(g ->data-> vertices, sizeof(int), n, out);
		if(!(check == n)){
				perror("Couldn't write to output file");
				exit(1);
			}
		g = g->next;
	}
	fclose(out);
}

int main(int argc, char **argv) {
	group_list* O;

	FILE *input, *output;
	spmat *A;
	if(!(argc==3)){
		perror("Invalid Input");
		exit(1);
	}

	input = fopen(argv[1], "rb");
	output = fopen(argv[2], "wb");
	A = build_a(input);
	fclose(input);
	O=Algorithm3(A);
	write_to_out(output, O);
	free_group_list(O);
	A->free(A);
	free_k();
	return 0;
}
