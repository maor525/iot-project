/*
 * group.c
 *
 *  Created on: Sep 9, 2020
 *      Author: Maor
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


void free_group_list(group_list *g){	/*frees the resources used by a group node*/
	group_node *tmp_g_node;
	while (g->head != NULL){
		tmp_g_node = g->head->next;
		free_group_node(g->head);
		g->head = tmp_g_node;
	}
	free(g);
}


void free_group_node(group_node *g){ /*frees the resources used by a group node*/
	free(g->data->vertices);
	free(g->data->f_g);
	free(g->data);
	free(g);

}

void add_group(group_list* g_list, group_node *g){ /*adds a group node to a group list*/
	if(g_list -> len ==0){
		g_list -> tail = g;
		g_list -> head = g;
	}
	else{
		g_list -> tail -> next = g;
		g_list -> tail = g;
	}
	g_list -> len++;
}

group_node* pop_group(group_list* g_list){ /* pops the first group node in a group list*/
	group_node *popped;
	popped = g_list -> head;
	if(g_list->len == 1){
		g_list -> head = NULL;
		g_list -> tail = NULL;
	}
	else{
		g_list -> head =g_list -> head -> next;
	}
	g_list -> len--;
	return popped;
}

void del_val(int *array, int val, int n){ /*deletes a value of an array and group the rest of the values*/
	bool swap = false;
	int i;
	for (i=0; i<n; i++){
		 if (array[i] == val){
			 swap = true;
		 }
		 if (swap) {
			 if (i==n-1){
				 continue;
			 }
			 array[i] = array[i+1];
		 }
	}
}

int arg_max(double* array, int n){ /*calculate the index of the entry with the largest value*/
	double max=array[0];
	int max_ind=0, i;
	for (i=1;i<n;i++){
		if (*(array+i) > max){
			max = *(array+i);
			max_ind = i;
		}
	}
	return max_ind;
}

