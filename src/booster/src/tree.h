/*

BOOSTER: BOOtstrap Support by TransfER: 
BOOSTER is an alternative method to compute bootstrap branch supports 
in large trees. It uses transfer distance between bipartitions, instead
of perfect match.

Copyright (C) 2017 Frederic Lemoine, Jean-Baka Domelevo Entfellner, Olivier Gascuel

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef _TREE_H_
#define _TREE_H_

#include "hashtables_bfields.h"	/* for the hashtables to store taxa names on the branches */
#include "hashmap.h"
#include "heavy_paths.h"
#include "io.h"
#include "M1M2_hashmap.h"
#include "split.h"
#include <ctype.h>
#include <stdbool.h>
#include <string.h>

#define	TRUE	1
#define	FALSE	0

#define	MIN_BRLEN	1e-8
#define MAX_TREELENGTH	100000000 /* more or less 10MB for a tree file in NH format */
#define MAX_MHEIGHT	100000 /* max mheight for nodes in the tree */

#define MAX_NAMELENGTH		255	/* max length of a taxon name */
#define MAX_COMMENTLENGTH	255	/* max length of a comment string in NHX format */

#define INCLUDE_EXCLUDE_SIZE 256 //The default size of the include/exlude arrrays.



/* TYPES */
typedef struct __NodeArray NodeArray;
typedef struct __Path Path;

/* Every node in our binary trees have several neighbours with indices  0, 1, 2.... We allow polytomies of any degree.
   An internal node with no multifurcation has 3 outgoing directions/neighbours.

   In rooted trees, the interpretation is the following:
   - for internal nodes, direction 0 is to the father, other directions are to the sons
   - for tips (leaves), direction 0 is to the father, no other neighbours
   - the root has only two neighbours.
   So it's easy to check whether a node is a leaf (one neighbour) or the root
   of a rooted tree (two neighbours).

   For unrooted trees, the interpretation is the same, except that no node has two neighbours.
   The pseudo-root (from the NH import) is three-furcated, so behaves exactly like any other internal node.

   It is not advisable to have several nodes of degree two in the same tree. 

   The root or pseudo-root is ALWAYS assigned id 0 at beginning. May change later, upon rerooting.
   In any case, it is always pointed to by tree->node0. */

typedef struct __Node Node;
typedef struct __Node {
	char* name;       /* Only set if this is a leaf node. */
	char* comment;		/* for further use: store any comment (e.g. from NHX format) */
	int id;			         /* unique id attributed to the node (index of node into a_nodes array*/
                           // Debugging Hint: ids are assigned according to to a pre-order traversal
                           // where the leaves are arranged in the order they come in the input file.
	short int nneigh;	      /* number of neighbours */
	short int nneigh_space; /* Currently allocated neighbors */
	struct __Node** neigh;	/* neighbour nodes */
	struct __Edge** br;	   /* corresponding branches going from this node */
	double mheight;	/* the height of a node is its min distance to a leaf */
  
        // Variables used for rapid transfer index calculation on alt_tree:
        // (used only in absence of heavypath decomposition in the alt_tree)
  int subtreesize; // Number of leaves in subtree rooted at this node (assume rooted)
  int depth;       // The depth of the node (from the root)
  int d_lazy;      // The lazily updated transfer distance
  int diff;        // For a node v, td(u,v) = d_lazy + Sum_{n \in Pv} diff_n
                   // (Pv is the path from v to the root)
  int d_min;       // Minimum TI found in this subtree
  int d_max;       // Maximum TI found in this subtree (used for unrooted TI)
  NodeArray* include; // Include these leaves in the transfer set for the subtree
  NodeArray* exclude; // Exclude these leaves from the transfer set for this node
  bool exclude_this;  // Used for leaf nodes when calculating the transfer set


  // Variables used for computation of transfer median
  split subtree_split; //Bitvector corresponding to the subtree.
  int leaf_id; //Ordered id of the leaf node.

        // Variables used for rapid transfer index calculation on the heavypath
        // decomposition for alt_tree:
  Path* path;      // Corresponding Path from the heavypath tree for this leaf.
                   // (allocated and freed in compute_transfer_indices_fast())
  

        // Variables used for rapid transfer index calculation on ref_tree:
  int k; // The number of best match
  int* ti_min;       // The (rooted) transfer index for this node.
  int* ti_max;       // The (rooted) maximum transfer distance for this node.
  Node** min_node;   // alt_tree node with minimum transfer distance to this one.
  Node** max_node;   // alt_tree node with maximum transfer distance to this one.
  NodeArray* lightleaves;   // The leaves in the light children.
  Node* heavychild; // The heaviest child
  Node* other;      // Corresponding leaf in another tree (see set_leaf_bijection())

  // Hash values
  unsigned int h1;
  unsigned int h2;
} Node;


/* Every edge connects two nodes.
   By convention, all terminal branches will have the tip on their RIGHT end */


typedef struct __Edge {
	int id;
	Node *left, *right; /* in rooted trees the right end will always be the descendant.
			       In any case, a leaf is always on the right side of its branch. */
	double brlen;
	double branch_support;
	int* subtype_counts[2];			/* first index is 0 for the left of the branch, 1 for its right side */
	id_hash_table_t *hashtbl;		/* bit vector containing the ids of the taxa in each subtree */
   
						/* index 0 corresponds to the left of the branch, index 1 to its right.
						   following our implementation, we only keep hashtbl[1] populated. */
	short int had_zero_length; 		/* set at the moment when we read the tree, even though
				      		   we then immediately set the branch length to MIN_BRLEN */
	short int has_branch_support; 	
	int topo_depth;				/* the topological depth is the number of taxa on the lightest side of the bipar */

         // Variables used for rapid transfer index calculation on ref_tree:
   int transfer_index;  // The (unrooted) transfer index for this edge.
} Edge;


typedef struct __Tree {
   int id; /* Unique ID of the tree. */
	Node** a_nodes;	/* array of node pointers */
	Edge** a_edges;	/* array of edge pointers */
	Node* node0;	/* the root or pseudo-root node */
   Node* outgroup_node; /* the node that is above the root.*/
   Edge* outgroup_edge; 
	int nb_nodes;
	int nb_edges;
	int nb_taxa;
  	int nb_taxa_space;  // Space currently allocated for taxa
  	int nb_edges_space; // Space currently allocated for edges
  	int nb_nodes_space; // Space currently allocated for nodes
  
	char** taxa_names; /* store only once the taxa names */
	char** taxname_lookup_table;   // maps taxa id [0,nb_taxa-1] to name

   //struct IntSet** M1M2hash;

         // Variables used for rapid transfer index calculation:
	NodeArray* leaves;       // Array of Node pointers sorted by name
   //bool* exclude_vector;    // Array used to mark leaf ids not in transfer set
   NodeArray* transfer_set; // The transfer set for the node with the min transfer index

   M1M2Hashtable* hashtable; //Pointer to the hashtable

   short int bipar_count; // If bipartition count will be incremented in the hashtable for this tree.
} Tree;


/* GLOBAL VARIABLE*/
extern int TREE_COUNT;
extern Tree** treelist;	

/* FUNCTIONS */

/* UTILS/DEBUG: counting specific branches or nodes in the tree */
int count_zero_length_branches(Tree* tree);
int count_leaves(Tree* tree);
int count_roots(Tree* tree);
int count_multifurcations(Tree* tree);
int dir_a_to_b(Node* a, Node* b);

/* various statistics on branch support */
double mean_bootstrap_support(Tree* tree);
double median_bootstrap_support(Tree* tree);
int summary_bootstrap_support(Tree* tree, double* result);

/* parsing utils: discovering and dealing with tokens */
int index_next_toplevel_comma(char* in_str, int begin, int end);
int count_outer_commas(char* in_str, int begin, int end);
void strip_toplevel_parentheses(char* in_str, int begin, int end, int* pair);
int index_toplevel_colon(char* in_str, int begin, int end);
void parse_double(char* in_str, int begin, int end, double* location);

/* creating a node, a branch, a tree: to create a tree from scratch, not from parsing */
Node* newNode(Tree* t, int k);
Node* new_node(const char* name, Tree* t, int degree);
Edge* new_edge(Tree* t);
Tree* new_tree(const char* name);
void addTip(Tree *t, char* name);
Edge* connect_to_father(Node* father, Node* son, Tree* current_tree);
Node* graft_new_node_on_branch(Edge* target_edge, Tree* tree, double ratio_from_left, double new_edge_length, char* node_name);

void addTip(Tree *t, char* name);
bool isNewickChar(char ch);
char parse_iter(Tree* t, char* in_str, int* position, int in_length, int* level, int k);


/* Replicate only the parts of the given tree important to the computation of
the rapid Transfer Index.
*/
Tree* copy_tree_rapidTI(Tree* t);
/* Copy the children of the old Node to the new Node.
@warning  assumes new is already innitialized with copy_node_rapidTI()
*/
void copy_tree_rapidTI_rec(Tree* t, Node* oldn, Node* newn);
/* Copy the Edge data essential to the rapid Transfer Index calculations.
*/
Edge* copy_edge_rapidTI(Edge *old, Node *parent, Node *child);
/* Copy the Node data essential to the rapid Transfer Index calculations.
*/
Node* copy_node_rapidTI(Node* old);


/* collapsing a branch */
void collapse_branch(Edge* branch, Tree* tree);

/**
   This function removes a taxon from the tree (identified by its taxon_id)
   And recomputed the branch length of the branch it was branched on.

   Also, the bootstrap support (if any) is recomputed, taking the maximum of the
   supports of the joined branches

   Be careful: The taxnames_lookup_table is modified after this function!
   Do not use this function if you share the same taxnames_lookup_table in
   several trees.
*/
void remove_taxon(int taxon_id, Tree* tree);

void addback_root_leaf(Tree* tree);


/**
   If there remains 2 neighbors to connect_node
   We connect them directly and delete connect_node
   We keep l_edge and delete r_edge
   -> If nneigh de connect node != 2 : Do nothing
   This function deletes a node like that:
              connect_node
             l_edge   r_edge
     l_node *-------*--------* r_node
   => Careful: After this function, you may want to call 
   => recompute_identifiers()

*/
void remove_single_node(Tree *tree, Node *connect_node);
/**
   This method recomputes all the identifiers 
   of the nodes and of the edges
   for which the tree->a_nodes is not null
   or tree->a_edges is not null
   It also recomputes the total number of edges 
   and nodes in the tree
 */
void recompute_identifiers(Tree *tree);

/**
   This function shuffles the taxa of an input tree 
*/
void shuffle_taxa(Tree *tree);

/* (re)rooting a tree */
void reroot_acceptable(Tree* t);
void unrooted_to_rooted(Tree* t);

/* To be called after a reroot*/
void reorient_edges(Tree *t);
void reorient_edges_recur(Node *n, Node *prev, Edge *e);

/* utility functions to deal with NH files */
unsigned int tell_size_of_one_tree(char* filename);
int copy_nh_stream_into_str(FILE* nh_stream, char* big_string);

/* actually parsing a tree */
void process_name_and_brlen(Node* son_node, Edge* edge, Tree* current_tree, char* in_str, int begin, int end);
Node* create_son_and_connect_to_father(Node* current_node, Tree* current_tree, int direction, char* in_str, int begin, int end);
void parse_substring_into_node(char* in_str, int begin, int end, Node* current_node, int has_father, Tree* current_tree);
Tree* parse_nh_string(char* in_str, int k);
char parse_recur(Tree* t, char* in_str, int* position, int in_length, Node* node, Edge* edge, int* level);
/* complete parse tree: parse NH string, update hashtables and subtype counts */
Tree *complete_parse_nh(char* big_string, char*** taxname_lookup_table,
                        bool skip_hashtables, int k, M1M2Hashtable* M1M2hashtable, int add_hash_count);


/* Functions for finding edge corresponding to hash*/
Edge* find_edge_from_hash_recur(Tree* tree, Node* current, unsigned int h1, unsigned int h2);
Edge* retrieve_edge_from_hash(Tree* tree, unsigned int h1, unsigned int h2);
Edge* retrieve_edge(unsigned int h1, unsigned int h2, M1M2Hashtable* hashtable, Tree** treelist);

/* taxname lookup table functions */
char** build_taxname_lookup_table(Tree* tree);
map_t build_taxid_hashmap(char** taxname_lookup_table, int nb_taxa);
void free_taxid_hashmap(map_t taxmap);
int free_hashmap_data(any_t arg,any_t key, any_t elemt);

char** get_taxname_lookup_table(Tree* tree);
Taxon_id get_tax_id_from_tax_name(char* str, char** lookup_table, int length);

/* (unnecessary/deprecated) multifurcation treatment */
void regraft_branch_on_node(Edge* branch, Node* target_node, int dir);

/***************************************************************
  ******************* neatly implementing tree traversals ******
***************************************************************/

void post_order_traversal_recur(Node* current, Node* origin, Edge *e, Tree* tree, void (*func)(Node*, Node*, Edge*, Tree*));
void post_order_traversal(Tree* t, void (*func)(Node*, Node*, Edge *e, Tree*));

/* post order traversal with any data passed to the function call */
void post_order_traversal_data_recur(Node* current, Node* origin, Edge* e, Tree* tree, void*, void (*func)(Node*, Node*, Edge*, Tree*, void*));
void post_order_traversal_data(Tree* t, void*, void (*func)(Node*, Node*, Edge*, Tree*, void*));


void pre_order_traversal_recur(Node* current, Node* origin, Edge *e, Tree* tree, void (*func)(Node*, Node*, Edge*, Tree*));
void pre_order_traversal(Tree* t, void (*func)(Node*, Node*, Edge*, Tree*));
/* pre order traversal starting at the given node */
void pre_order_traversal_subtree(Tree* t, Node* startnode, void (*func)(Node*, Node*, Edge*, Tree*));

/* pre order traversal with any data passed to the function call */
void pre_order_traversal_data_recur(Node* current, Node* origin, Edge *e, Tree* tree, void* data, void (*func)(Node*, Node*, Edge*, Tree*, void*));
void pre_order_traversal_data(Tree* t, void* data, void (*func)(Node*, Node*, Edge*, Tree*, void*));

/* bootstrap values */
void update_bootstrap_supports_from_node_names(Tree* tree);
void update_bootstrap_supports_doer(Node* current, Node* origin, Edge * e, Tree* tree);


/* node heights */

void update_node_heights_post_doer(Node* target, Node* orig, Edge *e, Tree* t);
void update_node_heights_pre_doer(Node* target, Node* orig, Edge *e, Tree* t);
void update_node_heights_post_alltree(Tree* tree);
void update_node_heights_pre_alltree(Tree* tree);

/* tree depth */
/* Set the depth of all the nodes of the tree. */
void update_node_depths_pre_doer(Node* target, Node* orig, Edge *e, Tree* t);
void prepare_rapid_TI_pre(Tree* tree);


/* topological depths */
void update_all_topo_depths_from_hashtables(Tree* tree);
int greatest_topo_depth(Tree* tree);

/* WORKING WITH HASHTABLES */

void update_hashtables_post_doer(Node* current, Node* orig, Edge *e, Tree* t);

void update_hashtables_post_alltree(Tree* tree);
void update_hashtables_pre_alltree(Tree* tree);


/* UNION AND INTERSECT CALCULATIONS FOR THE TRANSFER METHOD (from Bréhélin/Gascuel/Martin 2008) */
void update_i_c_post_order_ref_tree(Tree* ref_tree, Node* orig, Node* target, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix);
void update_all_i_c_post_order_ref_tree(Tree* ref_tree, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix);

void update_i_c_post_order_boot_tree(Tree* ref_tree, Tree* boot_tree, Node* orig, Node* target, short unsigned** i_matrix, short unsigned** c_matrix, short unsigned** hamming, short unsigned* min_dist, short unsigned* min_dist_edge);
void update_all_i_c_post_order_boot_tree(Tree* ref_tree, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix, short unsigned** hamming, short unsigned* min_dist, short unsigned* min_dist_edge);


/*Generate Random Tree*/
/**
   - nbr_taxa: Number of leaves in the output tree
   - taxa_names: array of leaf names: must be NULL if you do not want to use it 
   (names will be numbered in this case)

   - The output tree has branch lengths attributed using a normal distribution  N(0.1,0.05), and any br len < 0 is set to 0
*/
Tree * gen_rand_tree(int nbr_taxa, char **taxa_names);

/* writing a tree */

void write_nh_tree(Tree* tree, FILE* stream, bool newline);
void write_subtree_to_stream(Node* node, Node* node_from, Edge *e, FILE* stream);

/* freeing stuff */

void free_edge(Edge* edge);
void free_node(Node* node);
void free_tree(Tree* tree);


/* ____________________________________________________________ */
/* Functions added for rapid computation of the Transfer Index. */

/* NodeArray - - - - - - */
/*
An array of Node* along with its length.
Double the size of the array dynamically.
*/
typedef struct __NodeArray {
  Node** a;      //The array of node pointers
  int n;         //The length of the leaf array
  int i;         //The left-most unused index
} NodeArray;

/* Allocate a NodeArray of this size.
*/
NodeArray* allocateNA(int n);

/* Return a copy of the given NodeArray.
*/
NodeArray* copyNA(NodeArray *la);

/* Add a leaf to the leaf array.
*/
void addNodeNA(NodeArray *la, Node *u);

/* Remove the last leaf from the array.
*/
void removeNodeNA(NodeArray *la);

/* Clear the array.
*/
void clearNA(NodeArray *la);

/* Free the array in the NodeArray.
*/
void freeNA(NodeArray *la);

/* Print the nodes in the NodeArray.
*/
void printNA(NodeArray *la);

/* Sort by the taxa names.
*/
void sortNA(NodeArray *la);

/* Concatinate the given NodeArray. Free the memory of la1 and la2 if freemem
is true.

@note  user responsible for memory.
*/
NodeArray* concatinateNA(NodeArray *la1, NodeArray *la2, bool freemem);

/* Append the elemnts of la2 to la1.
*/
void appendNA(NodeArray *la1, NodeArray *la2);

/* - - - - - - - - - - - - - */


/* Do everything necessary to prepare for rapid Transfer Index (TI) computation.
*/
void prepare_rapid_TI(Tree* mytree);

/* Set the .other members for the leaves of the trees.

@warning  depends on leaves being in the same order for the two trees
*/
void set_leaf_bijection(Tree* tree1, Tree* tree2);


/*
Return all leaves coming from the light subtree of this node.

@warning  user responsible for memory
*/
NodeArray* get_leaves_in_light_subtree(Node *u);

/*
Find the heaviest child of this node (set u->heavychild), set u->lightleaves
to point to a NodeArray with all leaves not in the heavychild.

@warning  user responsible for memory of u->lightleaves (use freeLA())
*/
void setup_heavy_light_subtrees(Node *u);


/* Return a list of Node pointers to the leaves of this subtree.

@warning  user responsible for memory
*/
NodeArray* get_leaves_in_subtree(Node *u);


/* Add a leave to the leafarray, otherwise recurse.
*/
void add_leaves_in_subtree(Node *u, NodeArray *leafarray);


/* Return an array of indices to leaves in the node list.

@warning  user responsible for the memory of the returned array.
*/
int* get_leaf_indices(const Tree* tree);

/* Return a Node* array of all the leaves in the Tree.

@warning  user responsible for the memory of the returned array.
*/
Node** get_leaves(const Tree* tree);

/* Update the subtree size of the target node.

@warning  assumes binary rooted tree.
*/
void prepare_rapid_TI_doer(Node* target, Node* orig, Edge *e, Tree* t);

/* Set subtree size for all nodes.

@warning  assumes binary rooted tree.
*/
void prepare_rapid_TI_post(Tree* tree);


/* Return true if the given Node is the right child of its parent.

@warning  assume u is not the root.
*/
bool is_right_child(const Node* u);


/* Return true if the given pair of nodes represent the same taxon (possibly in
different trees).
*/
bool same_taxon(const Node *l1, const Node *l2);

/* Compare Nodes, using the string name of the node, so that leaves can be sorted.
Return <0 if n1 should go before, 0 if equal, and 1 if n1 should go after n2.

@warning  assume we are given leaves
*/
int compare_nodes(const void *l1, const void *l2);

/* Compare Nodes, using the bitarray on the edge, so that leaves can be sorted.
Return <0 if n1 should go before, 0 if equal, and 1 if n1 should go after n2.

@warning  assume we are given leaves
*/
int compare_nodes_bitarray(const void *l1, const void *l2);


/* Return the sibling to this Node.

@warning  assume the node is not the root and the tree is binary.
*/
Node* get_sibling(Node* u);

/* Return the siblings to this Node. Return NULL if n is the root.
*/
NodeArray* get_siblings(Node* n);

/* Return the sibling to this Node that is not the given Node sib.
Returns NULL if there is not another sibling (the root is not a pseudo-root).

@warning  assume the node is not the root
*/
Node* get_other_sibling(Node *n, Node *sib);

/* Return the node of alt_tree that has the current minimum transfer distance.
*/
Node* get_min_node(Tree* t);

/* Return the node of alt_tree that has the current maximum transfer distance.
*/
Node* get_max_node(Tree* t);

/* Return the Node that has the minimum transfer distance, if min is true.
Otherwise return the Node with the max. Work our way from the root to a leaf
and stop if we are at a root or all subtrees are not as good.
*/
Node* get_x_node(Tree* t, bool min);

/* Return the transfer set for the node. 

@note  user responsible for the memory
*/
NodeArray* get_transfer_set(Tree* t);

/* Return the transfer index on the tree.
*/
int transfer_index(Tree* t);

/* Return true if the min value represents the transfer index for the given
tree, rather than the max value.
*/
bool transfer_index_is_min(Tree* t);

/* Return the transfer set for the given node. If usemax is true, then
get the leaf set corresponding to the max value.
node).

@note  you must free any memory associated to t->transfer_set before calling
       this function
*/
NodeArray* get_transfer_set_for_node(Tree* t, Node* n, bool usemax);

/* Return the complement the transfer_set.

@note  user responsible for the memory
*/
NodeArray* get_complement_tset(Tree* t);

/* Add to n->transfer_set all of those leaves in the subtree that are not
in the n->exclude array.

@note  allocateLA() must already have been called on n->transfer_set 
*/
void add_transferset_from_subtree(Tree* t, Node* n, NodeArray* la);

/* Descend until the node with the best transfer index, adding the included
leaves to the given NodeArray.
*/
Node* collect_included(Node* n, NodeArray* includearray);

/* Get the leaves from all of the subtrees ABOVE this node, while subtracting
the leaves from the included array for that subtree (the leaves that should
be include in this subtree transfer set (TS) should be excluded from the set
above this node).
*/
Node* collect_included_above(Tree* t, Node* n, NodeArray* includearray);

/* Add all leaves from the subtree to the transfer_set, except those with the
id marked true in the exclude_vector.
*/
void include_subtree(Node* current, Node* previous, Edge *e, Tree* tree);

/* Return the transfer distance of the node. If min is true, then get the
value d_min + Sum_{n \in Pv} diff_n. Otherwise use d_max.
*/
int transfer_distance(Tree*t, Node* n);

/* Return the minimum of two integers.
*/
int min(int i1, int i2);
/* Return the maximum of two integers.
*/
int max(int i1, int i2);
/* Return the minimum of three integers.
*/
int min3(int i1, int i2, int i3);
/* Return the minimum of four integers.
*/
int min4(int i1, int i2, int i3, int i4);
/* Return the maximum of three integers.
*/
int max3(int i1, int i2, int i3);
/* Return the maximum of four integers.
*/
int max4(int i1, int i2, int i3, int i4);


/* Return the argmin of two integers
*/
int argmin(int i1, int i2);
/* Return the argmax of two integers
*/
int argmax(int i1, int i2);

/* Return the argmin of three integers
*/
int argmin3(int i1, int i2, int i3);
/* Return the argmax of three integers
*/
int argmax3(int i1, int i2, int i3);


/* - - - - - - - - - - - - - Printing Trees - - - - - - - - - - - - - - - */


/* Print all nodes of the tree in a post-order traversal.
*/
void print_nodes_post_order(Tree* t);
void print_node_callback(Node* n, Node* m, Edge *e, Tree* t);
void print_node(const Node* n);
void print_node_TI(const Node* n);
/* Print TI variables for a node in alt_tree.
*/
void print_node_TIvars(const Node* n);

/* Print the nodes from the given Node* array.
*/
void print_nodes(Node** nodes, const int n);
/* Print the nodes from the given Node* array (with the transfer index).
*/
void print_nodes_TI(Node** nodes, const int n);
/* Print the TI variables for the given nodes from alt_tree.
*/
void print_nodes_TIvars(Node** nodes, const int n);


/* Print the tree in dot format to the given filename.
*/
void print_tree_dot(Tree* t, char* filename, bool is_reftree);
void rec_print_ref_tree_dot(Node* n, FILE* f);
void rec_print_alt_tree_dot(Node* n, FILE* f);

/*
Print the alt_tree with the given index appended to the given filename.
*/
void print_alt_tree_dot(Tree* t, char* fileprefix, int index);

/* Print a node that describes the values in the positions of the alt_tree.
*/
void print_alttree_keynode_dot(FILE* f);


/* - - - - - - - - - - - - - Using Heavy Paths - - - - - - - - - - - - - - - */

/* Verify that all the leaves were reached in the heavypath decomposition.
*/
void verify_all_leaves_touched(Tree* t);


/* newick writer*/
void write_nh_tree_to_buffer(Tree* tree, char** buffer, size_t* capacity, size_t* offset, bool newline);

/* Pruning */
Tree* prune_branches(Tree* t, int* node_ids, int n_prune);


/* hash computation*/

unsigned int compute_hash_recur(Node* current);

/* reallocation of k mins and maxs*/
void realloc_nodes_k(Tree *t, int k);


#endif
