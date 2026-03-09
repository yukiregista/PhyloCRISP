

# include "tree.h"
# include "transfer_distance.h"

# include <assert.h>



/* Global Var*/
transferDistanceAllList* TDA_LIST = NULL;

/* Compute all transfer distance between a node and a tree*/
transferDistanceAll* transfer_distance_all_fast(Edge* edge, Tree* tree, int n_taxa, unsigned int node_h1, unsigned int node_h2){

    Node* child_node = edge->right;
    /* For simplicity, we always assign zero to the subtree of child_node.*/

    int m = tree->nb_edges;
    int* countOnesSubtree = calloc(n_taxa, sizeof(int));

    // Post-order traversal of the child_node to set 1 to the countOnesSubtree of leaf ids in the subtree;
    simple_post_order_traversal_with_intarray_recur(child_node, edge->left, countOnesSubtree, set_leaf_ones);

    // Post-order traversal of tree and compute transfer distance;
    int* transfer_distance = calloc(tree->nb_nodes, sizeof(int));
    update_transfer_distance(tree->node0, NULL, countOnesSubtree, transfer_distance, child_node->subtreesize, n_taxa);

    // Construct transfer_distance_all instance
    transferDistanceAll* tda = malloc(sizeof(transferDistanceAll));
    tda->td = calloc(tree->nb_nodes - tree->nb_taxa, sizeof(transferDistance*));
    tda->topo_depth = min(child_node->subtreesize, n_taxa - child_node->subtreesize);
    tda->node_h1 = node_h1;
    tda->node_h2 = node_h2;
    int count = 0;
    for (int i=0; i<tree->nb_nodes; i++){
        if(tree->a_nodes[i]!= NULL && tree->a_nodes[i]->nneigh!=1){
            // Internal node present in the tree
            tda->td[count] = malloc(sizeof(transferDistance));
            tda->td[count]->dist = transfer_distance[tree->a_nodes[i]->id]; // tree->a_nodes[i]->id should be equal to i.
            tda->td[count]->h1 = tree->a_nodes[tree->a_nodes[i]->id]->h1;
            tda->td[count]->h2 = tree->a_nodes[tree->a_nodes[i]->id]->h2;
            tda->td[count]->node_id = i;
            count ++;
        }
    }
    tda->n_elements = count;
    tda->td = realloc(tda->td, count * sizeof(transferDistance*));

    free(transfer_distance);
    free(countOnesSubtree);

    return tda;
}


void set_leaf_ones(Node* current, Node* origin, int* countOnesSubtree){
    if (current->nneigh==1){
        countOnesSubtree[current->leaf_id] = 1;
    }

}
int update_transfer_distance(Node* current, Node* origin, int* countOnesSubtree, int* td_array, int subtreesize, int n_taxa){
    if (current->nneigh == 1){
        int rtd = subtreesize - countOnesSubtree[current->leaf_id] + (1-countOnesSubtree[current->leaf_id]); // rooted transfer distance
        td_array[current->id] = min(rtd, n_taxa - rtd);
        return countOnesSubtree[current->leaf_id];
    }else{
        int start_index = 1;
        if (origin==NULL) start_index=0; // root case.
        int countOnes = 0;
        for (int i=start_index; i<current->nneigh; i++){
            /* Compute all transfer distance to the subtree, increment the countOnes of the child.*/
            countOnes += update_transfer_distance(current->neigh[i], current, countOnesSubtree, td_array, subtreesize, n_taxa);
        }
        int rtd = subtreesize - countOnes + (current->subtreesize - countOnes);
        td_array[current->id] = min(rtd, n_taxa-rtd);
        return countOnes;
    }
}


void simple_post_order_traversal_with_intarray_recur(Node* current, Node* origin, int* intarray, void (*func)(Node*, Node*, int*)) {
	/* does the post order traversal on current Node and its "descendants" (i.e. not including origin, who is a neighbour of current */
	int i, n = current->nneigh;

	for(i=0; i < n; i++){
		if(current->neigh[i] != origin){
			simple_post_order_traversal_with_intarray_recur(current->neigh[i], current, intarray, func);
		}
	}

	/* and then in any case, call the function on the current node */
	func(current, origin, intarray);
}


void simple_post_order_traversal2(Node* current, Node* origin, int* intarray1, int* intarray2, int subtreesize, int n_taxa,  void (*func)(Node*, Node*, int*, int*, int, int)){
    /* does the post order traversal on current Node and its "descendants" (i.e. not including origin, who is a neighbour of current */
	int i, n = current->nneigh;

	for(i=0; i < n; i++){
		if(current->neigh[i] != origin){
			simple_post_order_traversal2(current->neigh[i], current, intarray1, intarray2, subtreesize, n_taxa, func);
		}
	}

	/* and then in any case, call the function on the current node */
	func(current, origin, intarray1, intarray2, subtreesize, n_taxa);
}

void free_transfer_distance_all(transferDistanceAll* tda){
    for (int i=0; i<tda->n_elements; i++){
        free(tda->td[i]);
    }
    free(tda->td);
    free(tda);
}

void init_TDA_LIST(){
    TDA_LIST = malloc(sizeof(transferDistanceAllList));
    TDA_LIST->num_entries = 0;
    int num_memory = 100; // memory allocated.
    TDA_LIST->num_memory = num_memory;
    TDA_LIST->tda_list = malloc(num_memory * sizeof(transferDistanceAll*));
    for (int i=0; i<num_memory; i++) TDA_LIST->tda_list[i] = NULL; //initialize with NULL pointer
}

void TDA_LIST_append(transferDistanceAll* tda){
    if (TDA_LIST->num_entries == TDA_LIST->num_memory){
        int new_memory_blocks = 2 * TDA_LIST->num_memory;
        TDA_LIST->tda_list = realloc(TDA_LIST->tda_list, new_memory_blocks * sizeof(transferDistanceAll*)); // Double memory.
        for (int i=TDA_LIST->num_memory; i<new_memory_blocks; i++) TDA_LIST->tda_list[i] = NULL; //initialize with NULL pointer
        TDA_LIST->num_memory = new_memory_blocks;
    }
    TDA_LIST->tda_list[TDA_LIST->num_entries] = tda;
    TDA_LIST->num_entries ++;
}