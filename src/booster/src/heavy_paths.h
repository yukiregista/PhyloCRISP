#ifndef __HEAVY_PATHS_H__
#define __HEAVY_PATHS_H__
/*
  This is the code that uses heavy paths to help when the alternate tree
  (T_o in the paper) is not balanced. We use bookkeeping, instead of a
  balanced datastructure, to maintain the minimum transfer distance over the
  whole tree. This shaves a logarithmic factor off of the running time
  when traveling up the alternate tree from leaf to root.
 
  In the alternate tree, a heavy path is split into a tree ("Path Tree" or
  PT) where each node of the tree represents a subpath of the heavypath. Each
  subpath knows its min/max value on the subpath, and its min/max value for
  all the subtrees hanging off the subpath. This way, when a heavy path is
  updated the min value on the path will decrease by 1, but the min value off
  the path increase by 1.

  The PTs are joined together into a "HeavyPath Tree" (HPT), where a root of
  a PT is connected to a leaf of another PT using the parent_heavypath and
  child_heavypath pointers.
*/

#include "tree.h"
#include "limits.h"
#include "debug.h"
#include "stdbool.h"
#include "split.h"


typedef struct __Node Node;

/*
  We use the following conventions:
  - the original alt_tree is decomposed into heavypaths
  - a heavypath is represented by a binary tree of Path objects (called the
    PathTree PT) where:
    * the root of the tree contains information for the entire heavypath
    * the leaves of the tree contain node information for the corresponding
      node in alt_tree, along with pointers (child_heavypaths) to the
      pendant heavypaths (in the case of a binary tree there is only one
      pendant heavypath).

  In other words, a tree of Paths (PT) represents a heavypath, where the Path
  object can be an internal subpath node, the root (which contains the
  summary bookkeeping for the entire heavypath), or a leaf.
  In the case of a PT leaf, the Path represents the node in alt_tree on the
  current heavypath, pointed to by node, and the pendant heavypath is pointed
  to by child_heavypath.

  The entire group of PTs, that are glued together by child_heavypath pointers,
  is called the HeavyPathTree (HPT).

  If the PT leaf is also a leaf of alt_tree, then n->node and
  n->child_heavypaths are NULL and the path_to_root is an array of Path objects
  that represent the path in the HPT from the leaf to the root of the HPT.
*/
typedef struct __Path Path;
typedef struct __Path {
  int id;          //A unique identifier.

    //The structure of the heavypath tree:
  Path* left;      //Left child.
  Path* right;     //Right child.
  Path* parent;
  Path* sibling;

  Node* node;               //The node of alt_tree corresponding to this Path.
                            // (this applies only to leaves of the PT)

  Path** child_heavypaths;  //The array of Path trees pendant to this Path.
  int num_child_paths;      //The size of the child_heavypaths array.
  Path* parent_heavypath;   //The Path that this PT hangs on.

  int total_depth;          //# of Path structs to root through all PTs
                            // (i.e. # nodes in path to HPT root).

  Path*** path_to_root_p;   //Pointer to array of Path objects to root.
                            // (this points to memory shared between all HPT
                            //  leaves -- it is set once when used)
                            // When calling add_leaf_HPT(leaf), this will be renewed to represent
                            // the path from that leaf to the root.

    //The transfer index (TI) values:
  int diff_path;     //Diff to add to subtree rooted on path.
  int diff_subtree;  //Diff to add to pendant subtrees.


  int k;             // The number of minimums and maximums that it holds. Specified at runtime.
  int* order_array; // Auxiliary array to be used for storing order of picking.

  int* d_min_path;    //First k Min value for the subpath. For a leaf node v in PT, it also contains D[v] value.
  int* d_min_subtree;  //First k Min value over all pendant subtrees for this (sub)Path. 


  int* d_max_path;    //First k Max value for the subpath.
  int* d_max_subtree; //First k Max value over all pendant subtrees for this (sub)Path. 


  Node** d_min_path_node; // Array of k Nodes corresponding to first k min_path
  Node** d_min_subtree_node; // Array of k Nodes corresponding to first k min_subtree

  Node** d_max_path_node; // Array of k Nodes corresponding to first k max_path
  Node** d_max_subtree_node; // Array of k Nodes corresponding to first k max_subtree
} Path;


/*
Safe addition.
*/
int safe_add(int a, int b);

/*
Compare two sorted ascending arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_ascending(int* a, int*b, int* result, int a_length, int b_length, int result_length);

/*
Compare two sorted ascending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_ascending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* result, int a_length, int b_length, int result_length);


/*
Compare two sorted descending arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_descending(int* a, int*b, int* result, int a_length, int b_length, int result_length);

/*
Compare two sorted descending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_descending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* result, int a_length, int b_length, int result_length);

/*
Set the array from two candidates, using the order created by 'compare_two_sorted_...' functions.
It will not do additional checking of validity of order.
*/
void set_int_array_from_order_2(int* a, int* b, int* result, int* order_array, int result_length);
void set_int_array_from_order_2_with_diffs(int* a, int a_diff, int* b, int b_diff, int* result, int* order_array, int result_length);
void set_node_array_from_order_2(Node** a, Node** b, Node** result, int* order_array, int result_length);
void set_split_array_from_order_2(Node**a, Node**b, split* result, int* order_array, int result_length, int split_len);

void set_int_array_from_order_3(int* a, int* b, int* c, int* result, int* order_array, int result_length);
void set_int_array_from_order_3_with_diffs(int* a, int a_diff, int* b, int b_diff, int* c, int c_diff, int* result, int* order_array, int result_length);
void set_node_array_from_order_3(Node** a, Node** b, Node** c, Node** result, int* order_array, int result_length);


/*
Compare three sorted ascending arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_ascending(int* a, int*b, int* c, int* result, int a_length, int b_length, int c_length, int result_length);

/*
Compare three sorted ascending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_ascending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* c, int c_diff, int* result, int a_length, int b_length, int c_length, int result_length);


/*
Compare three sorted descending arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_descending(int* a, int*b, int* c, int* result, int a_length, int b_length, int c_length, int result_length);

/*
Compare three sorted descending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_descending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* c, int c_diff, int* result, int a_length, int b_length, int c_length, int result_length);

/*
Set the array from three candidates, using the order created by 'compare_three_sorted_...' functions.
It will not do additional checking of validity of order.
*/
void set_int_array_from_order_3(int* a, int* b, int* c,  int* result, int* order_array, int result_length);
void set_node_array_from_order_3(Node** a, Node** b, Node** c, Node** result, int* order_array, int result_length);

void update_int_array(int* a, int* result, int d, int result_length);
void copy_node_array(Node** a, Node** result, int result_length);


/*
Allocate a new Path, setting all the default values.
*/
Path* new_Path(int k);

/* Recursiveley decompose the alternate tree into heavy paths according to
the scheme described in the definition of the Path struct. Return the root
Path of the Path tree.

@note   each heavypath corresponds to a tree of Paths we call the PathTree (PT)
        leaves of the PTs are glued the roots of other PTs (using the
        child_heavypath pointer).
        We call the entire tree the HeavyPathTree (HPT)
*/
Path* do_heavy_decomposition(Node* root, int k);
Path* heavy_decomposition(Node* root, int depth, int* maxdepth, int k,
                          Path*** path_to_root_pointer, int* order_array);

/* Free the memory for the HeavyPathTree (allocated in heavy_decomposition).
*/
void free_HPT(Path* root);
void free_HPT_rec(Path* node);

/* For the given heavypath, create a Path structure that represents the path.
Split the path in half and create a Path for each half.  If a half is a single
node, then hang the next heavy path off of it. If it's a leaf of alt_tree, then
link the Path to the corresponding leaf in alt_tree.
*/
Path* partition_heavypath(Node** n, int length, int depth, int* maxdepth, int k,
                          Path*** path_to_root_pointer, int* order_array);

/*
Return a Path for the given node of alt_tree.  The Path will be a leaf
node of the path tree.
Either 1) the leaf will point to a leaf node of alt_tree,
or     2) child_heavypath will point to a heavypath representing the
          descendant of the alt_tree node.
*/
Path* heavypath_leaf(Node* node, int depth, int* maxdepth, int k,
                     Path*** path_to_root_pointer, int* order_array);

/* Return the heavypath rooted at the node.
   ** user responsible for memory of returned heavypath.
*/
Node** get_heavypath(Node* root, int* length);

/* Return the length of the heavy path rooted at this node.
*/
int get_heavypath_length(Node* n);

/* Return True if the (sub)Path corresponds to a leaf of alt_tree (i.e. it is a
leaf of the HPT).
*/
bool is_HPT_leaf(Path* n);

/* Print the given heavypath.
*/
void print_heavypath(Node** heavypath, int length);


/*
Descend to the leaves of the HPT. Once there, create the path_to_root array
for that Path object.
*/
//void set_paths_to_root(Path* node);

/* Build a path from this Path leaf up to the root of the HPT, following each
PT to it's root in turn.
*/
Path** get_path_to_root_HPT(Path* leaf);

/* Build a path from this Path leaf up to the root of the HPT, following each
PT to it's root in turn.
*/
void set_path_to_root_HPT(Path* leaf, Path** path_to_root);


/* Return the root of the HPT with the given alt_tree leaf.
*/
Path* get_HPT_root(Node* leaf);

/* Reset the path and subtree min and max, along with the diff values for the
path from the given leaf to the root of the HPT.
*/
void reset_leaf_HPT(Node* leaf);

/*--------------------- OUTPUT FUNCTIONS -------------------------*/


/* Print the Heavy Path Tree (HPT) in dot format.
*/
void print_HPT_dot(Path* hproot, Node* altroot, int repid);

/* Print the given Path node.
*/
void print_HPT_node(const Path* n);

/* Recursively print the subPath.
*/
void print_HPT_subpath_dot(Path* n, FILE* f);

/* Recursively print the subtree.
*/
void print_HPT_subtree_dot(Node* n, FILE* f);



/* Print a string representing this Path node formatted for dot output.
*/
void print_HPT_node_dot(Path* n, FILE* f);

/* Print a string that formats a heavypath node in a PT.
*/
void print_HPT_hpnode_dot(Path* n, FILE* f);

/* Print a string that formats a PT node.
*/
void print_HPT_ptnode_dot(Path* n, FILE* f);

/* Print a descriptive node.
*/
void print_HPT_keynode_dot(FILE* f);

#endif
