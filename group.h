/*
 * group.h
 *
 *  Created on: Sep 9, 2020
 *      Author: Maor
 */

#ifndef GROUP_H_
#define GROUP_H_
#include "spmat.h"

/*this module is used to handle groups and store them in memory,
 * the functions in this module are used to manage the splitting
 * and moving of groups*/

/*group struct, containing the number of vertices in the group,
 * array of vertices, A[g] subMatrix, and f[g]*/
typedef struct _group {
	struct _spmat			*A_g;
	double 			*f_g;
	int 			*vertices;
	int  			len;
} group;

/*a node in a list which contains groups, each node contains
 * a group and a pointer to the next group node */
typedef struct _group_node {
   group          *data;
   struct _group_node  *next;
} group_node;

/*a list containing groups, includes pointers to the start and
 * the end of the list aswell as its length*/
typedef struct group_linked_list{
	group_node *head;
	group_node *tail;
	int len;
} group_list;

void free_group_list(group_list *g);
void free_group_node(group_node *g);
void add_group(group_list* g_list, group_node *g);
group_node* pop_group(group_list* g_list);
void del_val(int *array, int val, int n);
int arg_max(double* array, int n);
#endif /* GROUP_H_ */
