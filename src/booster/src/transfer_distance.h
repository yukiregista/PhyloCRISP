
#ifndef _TRANSFER_DISTANCE
#define _TRANSFER_DISTANCE

# include "tree.h"

typedef struct __transfer_distance{
    int dist;
    unsigned int h1;
    unsigned int h2;
    int node_id;
} transferDistance;


typedef struct __transfer_distance_all{
    transferDistance** td;
    int n_elements;
    unsigned int node_h1;
    unsigned int node_h2;
    int topo_depth;
} transferDistanceAll;


typedef struct __transfer_distance_all_list{
    int num_entries;
    int num_memory;
    transferDistanceAll** tda_list; 
} transferDistanceAllList;

/* Global var to store transfer distance all list*/
extern transferDistanceAllList* TDA_LIST;


transferDistanceAll* transfer_distance_all_fast(Edge* edge, Tree* tree, int n_taxa, unsigned int node_h1, unsigned int node_h2);

void set_leaf_ones(Node* current, Node* origin, int* countOnesSubtree);

int update_transfer_distance(Node* current, Node* origin, int* countOnesSubtree, int* td_array, int subtreesize, int n_taxa);

void simple_post_order_traversal_with_intarray_recur(Node* current, Node* origin, int* intarray, void (*func)(Node*, Node*, int*));

void simple_post_order_traversal2(Node* current, Node* origin, int* intarray1, int* intarray2, int subtreesize, int n_taxa,  void (*func)(Node*, Node*, int*, int*, int, int));


void free_transfer_distance_all(transferDistanceAll* tda);

void init_TDA_LIST();
void TDA_LIST_append(transferDistanceAll* tda);


#endif