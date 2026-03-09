/*
Here we implement the heavy path decomposition of the alt_tree.
*/
#include "heavy_paths.h"

/*
Safe addition.
*/
int safe_add(int a, int b){
  int res;
  if (b>0 && a > INT_MAX-b) res=INT_MAX;
  else if (b<0 && a < INT_MIN - b) res=INT_MIN;
  else res=a+b;
  return res;
}

/*
Compare two sorted (ascending) arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_ascending(int* a, int*b, int* result, int a_length, int b_length, int result_length)
{
  if (result_length > a_length + b_length || a_length==0 || b_length==0) {
        fprintf(stderr, "Error: result_length must be less than or equal to the sum of a_length and b_length, and a_length, b_length should be more than zero.\n");
        exit(EXIT_FAILURE);
    }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
   while(i<result_length){
    if(b[k]<a[j]){
      result[i++]=1;
      ++k;
    }
    else {
      result[i++]=0;
      ++j;
    }
    if(j==a_length) {while(i<result_length) result[i++]=1;}//Now i=result_length, so it will exit the while loop
    else if(k==b_length) {while (i<result_length) result[i++]=0;}//Now i=result_length, so it will exit the while loop
   }
}

/*
Compare two sorted ascending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_ascending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* result, int a_length, int b_length, int result_length)
{
  if (result_length > a_length + b_length || a_length==0 || b_length==0) {
        fprintf(stderr, "Error: result_length must be less than or equal to the sum of a_length and b_length, and a_length, b_length should be more than zero.\n");
        exit(EXIT_FAILURE);
    }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
   while(i<result_length){
    if(safe_add(b[k],b_diff)<safe_add(a[j],a_diff)){
      result[i++]=1;
      ++k;
    }
    else {
      result[i++]=0;
      ++j;
    }
    if(j==a_length) {while(i<result_length) result[i++]=1;}//Now i=result_length, so it will exit the while loop
    else if(k==b_length) {while (i<result_length) result[i++]=0;}//Now i=result_length, so it will exit the while loop
   }
}

/*
Compare two sorted (ascending) arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_descending(int* a, int*b, int* result, int a_length, int b_length, int result_length)
{
  if (result_length > a_length + b_length || a_length==0 || b_length==0) {
        fprintf(stderr, "Error: result_length must be less than or equal to the sum of a_length and b_length, and a_length, b_length should be more than zero.\n");
        exit(EXIT_FAILURE);
    }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
   while(i<result_length){
    if(b[k]>a[j]){
      result[i++]=1;
      ++k;
    }
    else {
      result[i++]=0;
      ++j;
    }
    if(j==a_length) {while(i<result_length) result[i++]=1;}//Now i=result_length, so it will exit the while loop
    else if(k==b_length) {while (i<result_length) result[i++]=0;}//Now i=result_length, so it will exit the while loop
   }
}


/*
Compare two sorted descending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_two_sorted_descending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* result, int a_length, int b_length, int result_length)
{
  if (result_length > a_length + b_length || a_length==0 || b_length==0) {
        fprintf(stderr, "Error: result_length must be more than or equal to the sum of a_length and b_length, and a_length, b_length should be more than zero.\n");
        exit(EXIT_FAILURE);
    }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
   while(i<result_length){
    if(safe_add(b[k],b_diff)>safe_add(a[j],a_diff)){
      result[i++]=1;
      ++k;
    }
    else {
      result[i++]=0;
      ++j;
    }
    if(j==a_length) {while(i<result_length) result[i++]=1;}//Now i=result_length, so it will exit the while loop
    else if(k==b_length) {while (i<result_length) result[i++]=0;}//Now i=result_length, so it will exit the while loop
   }
}


/*
Update d_min or d_max arrays, caring for the case when they include INT_MIN or INT_MAX
*/

void update_int_array(int* a, int* result, int d, int result_length)
{
  if (d>0){
    for(int i=0; i<result_length;++i){
    if(!(a[i]>INT_MAX-d))result[i] = a[i] + d;
    else result[i]=a[i];
    }
  }
  else if (d<0){
    for(int i=0; i<result_length;++i){
    if(!(a[i]<INT_MIN-d))result[i] = a[i] + d;
    else result[i]=a[i];
    }
  }
  else{
    assert(d==0);
    for(int i=0; i<result_length;++i) result[i]=a[i];
  }
}

void copy_node_array(Node** a, Node** result, int result_length){
  for(int i=0; i<result_length; ++i) result[i]=a[i];
}

/*
Set the array from two candidates, using the order created by 'compare_two_sorted_...' functions.
It will not do additional checking of validity of order.
*/
void set_int_array_from_order_2(int* a, int* b, int* result, int* order_array, int result_length)
{
  int j=0; //index of a
  int k=0; //index of b 

  for (int i=0; i<result_length; ++i){
    if(order_array[i] == 0) result[i] = a[j++];
    else result[i]=b[k++];
  }
}

void set_int_array_from_order_2_with_diffs(int* a, int a_diff, int* b, int b_diff, int* result, int* order_array, int result_length)
{
  int j=0; //index of a
  int k=0; //index of b

  for (int i=0; i<result_length; ++i){
    if(order_array[i] == 0){
      result[i] = safe_add(a[j++], a_diff);
    }
    else{
      result[i] = safe_add(b[k++], b_diff);
    }
  }
}

/*
Set the array from two candidates, using the order created by 'compare_two_sorted_...' functions.
It will not do additional checking of validity of order.
*/
void set_node_array_from_order_2(Node** a, Node** b, Node** result, int* order_array, int result_length)
{
  int j=0; //index of a
  int k=0; //index of b

  for (int i=0; i<result_length; ++i){
    if(order_array[i] == 0) result[i] = a[j++];
    else result[i]=b[k++];
  }
}


/*
Set the array from two candidates, using the order created by 'compare_two_sorted_...' functions.
It will not do additional checking of validity of order.
*/
void set_split_array_from_order_2(Node** a, Node** b, split* result, int* order_array, int result_length, int split_len)
{
  int j=0; //index of a
  int k=0; //index of b

  for (int i=0; i<result_length; ++i){
    if(order_array[i] == 0) memcpy(result[i], a[j++]->subtree_split, split_len * sizeof(unsigned int));
    else memcpy(result[i], b[k++]->subtree_split, split_len * sizeof(unsigned int));
  }
}


/*
Compare three sorted ascending arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_ascending(int* a, int*b, int* c, int* result, int a_length, int b_length, int c_length, int result_length)
{
  if (result_length > a_length + b_length + c_length || a_length==0 || b_length==0 || c_length==0) {
        fprintf(stderr, "Error: result_length must be more than or equal to the sum of lengths, and each length should be more than zero.\n");
        exit(EXIT_FAILURE);
    }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
  int l=0; // index of l
   while(i<result_length){
    if(b[k]<a[j]){
      if (c[l]<b[k]){
        result[i++] = 2;
        ++l;
      }else{
        result[i++] = 1;
        ++k;
      }
    }
    else {
      if(c[l]<a[j]){
        result[i++] = 2;
        ++l;
      }else{
        result[i++]=0;
        ++j;
      }
    }
    if(j==a_length) {compare_two_sorted_ascending(&b[k], &c[l], &result[i], b_length-k, c_length-l, result_length-i); break; }
    else if(k==b_length) {compare_two_sorted_ascending(&a[j], &c[l], &result[i], a_length-j, c_length-l, result_length-i); break; }
    else if(l==c_length) {compare_two_sorted_ascending(&a[j], &b[k], &result[i], a_length-j, b_length-k, result_length-i); break; }
   }
}

/*
Compare three sorted ascending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_ascending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* c, int c_diff, int* result, int a_length, int b_length, int c_length, int result_length)
{
  if (result_length > a_length + b_length + c_length || a_length==0 || b_length==0 || c_length==0) {
      fprintf(stderr, "Error: result_length must be more than or equal to the sum of lengths, and each length should be more than zero.\n");
      exit(EXIT_FAILURE);
  }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
  int l=0; // index of l
   while(i<result_length){
    if(safe_add(b[k],b_diff)<safe_add(a[j],a_diff)){
      if (c[l]+c_diff<b[k]+b_diff){
        result[i++] = 2;
        ++l;
      }else{
        result[i++] = 1;
        ++k;
      }
    }
    else {
      if(safe_add(c[l],c_diff)<safe_add(a[j],a_diff)){
        result[i++] = 2;
        ++l;
      }else{
        result[i++]=0;
        ++j;
      }
    }
    if(j==a_length) {compare_two_sorted_ascending_with_diffs(&b[k], b_diff, &c[l], c_diff, &result[i], b_length-k, c_length-l, result_length-i); break; }
    else if(k==b_length) {compare_two_sorted_ascending_with_diffs(&a[j], a_diff, &c[l], c_diff, &result[i], a_length-j, c_length-l, result_length-i); break; }
    else if(l==c_length) {compare_two_sorted_ascending_with_diffs(&a[j], a_diff, &b[k], b_diff, &result[i], a_length-j, b_length-k, result_length-i); break; }
   }
}

/*
Compare three sorted descending arrays and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_descending(int* a, int*b, int* c, int* result, int a_length, int b_length, int c_length, int result_length)
{
  if (result_length > a_length + b_length + c_length || a_length==0 || b_length==0 || c_length==0) {
        fprintf(stderr, "Error: result_length must be more than or equal to the sum of lengths, and each length should be more than zero.\n");
        exit(EXIT_FAILURE);
    }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
  int l=0; // index of l
   while(i<result_length){
    if(b[k]>a[j]){
      if (c[l]>b[k]){
        result[i++] = 2;
        ++l;
      }else{
        result[i++] = 1;
        ++k;
      }
    }
    else {
      if(c[l]>a[j]){
        result[i++] = 2;
        ++l;
      }else{
        result[i++]=0;
        ++j;
      }
    }
    if(j==a_length) {compare_two_sorted_descending(&b[k], &c[l], &result[i], b_length-k, c_length-l, result_length-i); break; }
    else if(k==b_length) {compare_two_sorted_descending(&a[j], &c[l], &result[i], a_length-j, c_length-l, result_length-i); break; }
    else if(l==c_length) {compare_two_sorted_descending(&a[j], &b[k], &result[i], a_length-j, b_length-k, result_length-i); break; }
   }
}


/*
Compare three sorted descending arrays with diffs and return the order of picking from the two arrays that produce sorted array.
*/
void compare_three_sorted_descending_with_diffs(int* a, int a_diff, int*b, int b_diff, int* c, int c_diff, int* result, int a_length, int b_length, int c_length, int result_length)
{
    if (result_length > a_length + b_length + c_length || a_length==0 || b_length==0 || c_length==0) {
        fprintf(stderr, "Error: result_length must be more than or equal to the sum of lengths, and each length should be more than zero.\n");
        exit(EXIT_FAILURE);
    }
  int i=0; // index of result
  int j=0; // index of a
  int k=0; // index of b
  int l=0; // index of l
   while(i<result_length){
    if(safe_add(b[k],b_diff)>safe_add(a[j],a_diff)){
      if (c[l]+c_diff>b[k]+b_diff){
        result[i++] = 2;
        ++l;
      }else{
        result[i++] = 1;
        ++k;
      }
    }
    else {
      if(safe_add(c[l],c_diff)>safe_add(a[j],a_diff)){
        result[i++] = 2;
        ++l;
      }else{
        result[i++]=0;
        ++j;
      }
    }
    if(j==a_length) {compare_two_sorted_descending_with_diffs(&b[k], b_diff, &c[l], c_diff, &result[i], b_length-k, c_length-l, result_length-i); break; }
    else if(k==b_length) {compare_two_sorted_descending_with_diffs(&a[j], a_diff,  &c[l], c_diff, &result[i], a_length-j, c_length-l, result_length-i); break; }
    else if(l==c_length) {compare_two_sorted_descending_with_diffs(&a[j], a_diff, &b[k], b_diff, &result[i], a_length-j, b_length-k, result_length-i); break; }
   }
}

/*
Set the array from three candidates, using the order created by 'compare_three_sorted_...' functions.
It will not do additional checking of validity of order.
*/
void set_int_array_from_order_3(int* a, int* b, int* c,  int* result, int* order_array, int result_length)
{
  int j=0; //index of a
  int k=0; //index of b
  int l=0; //index of c

  for (int i=0; i<result_length; ++i){
    if(order_array[i]==0) result[i] = a[j++];
    else if(order_array[i]==1) result[i]=b[k++];
    else result[i]=c[l++];
  }
}
void set_node_array_from_order_3(Node** a, Node** b, Node** c, Node** result, int* order_array, int result_length)
{
  int j=0; //index of a
  int k=0; //index of b
  int l=0; //index of c

  for (int i=0; i<result_length; ++i){
    if(order_array[i]==0) result[i]=a[j++];
    else if(order_array[i]==1) result[i]=b[k++];
    else result[i]=c[l++];
  }
}

void set_int_array_from_order_3_with_diffs(int* a, int a_diff, int* b, int b_diff, int* c, int c_diff, int* result, int* order_array, int result_length)
{
  int j=0; //index of a
  int k=0; //index of b
  int l=0; //index of c

  for (int i=0; i<result_length; ++i){
    if(order_array[i]==0){
      result[i] = safe_add(a[j++],a_diff);
    }
    else if(order_array[i]==1){
      result[i]= safe_add(b[k++] , b_diff);
    }else{
      result[i] = safe_add(c[l++] , c_diff);
    }
  }
}


/*
Allocate a new Path, setting all the default values.
*/
int id_counter = 0;          //Global id counter for Path structs.
Path* new_Path(int k)
{
  Path *newpath = malloc(sizeof(Path));

  newpath -> k = k; // Number of minimums/maximums to store.

  newpath->id = id_counter++;

  newpath->left = NULL;
  newpath->right = NULL;
  newpath->parent = NULL;
  newpath->sibling = NULL;

  newpath->node = NULL;

  newpath->child_heavypaths = NULL;
  newpath->num_child_paths = 0;
  newpath->parent_heavypath = NULL;
  newpath->path_to_root_p = NULL;

  newpath->total_depth = 0;

  newpath->diff_path = newpath->diff_subtree = 0;


  newpath->d_min_path = (int*)malloc(k*sizeof(int));
  newpath->d_min_subtree = (int*)malloc(k*sizeof(int));
  newpath->d_max_path = (int*)malloc(k*sizeof(int));
  newpath->d_max_subtree = (int*)malloc(k*sizeof(int));

  // Initialize all values.
  for(int i=0; i < k; i++){
    newpath->d_min_path[i] = INT_MAX;
    newpath->d_min_subtree[i] = INT_MAX;
    newpath->d_max_path[i] = INT_MIN;
    newpath->d_max_subtree[i] = INT_MIN;
  }

  // newpath->d_min_path = 1; // minimum of d in the subtree of this pathtree
  // newpath->d_min_subtree = 1; /* minimum of d in the subtree of the whole HPT; 
  // After the initialization, this is always equal to 1 (since d[x,x] for any node x is 1.); 
  // This is why we don't update d_min_subtree when constructing heavy path tree.
  
  // Note that d_min_subtree for the HPT_leaf is meaningless, and will never be used.
  // (Because there is no pending subtree for HPT leaf).
  // */

  // newpath->d_max_path = 0;
  // newpath->d_max_subtree = 1;

  newpath->d_min_path_node = (Node**)calloc(k, sizeof(Node*)); // Initialized with NULL pointers
  newpath->d_min_subtree_node = (Node**)calloc(k, sizeof(Node*));

  newpath->d_max_path_node = (Node**)calloc(k, sizeof(Node*));
  newpath->d_max_subtree_node = (Node**)calloc(k, sizeof(Node*));

  return newpath;
}


/*
Recursiveley decompose the alternative tree into heavy paths according to
the scheme described in the definition of the Path struct. Return the root
Path of the Path tree.

@note   Each heavypath corresponds to a tree of Paths we call the PathTree (PT).
        Leaves of the PTs are glued to the roots of other PTs (using the
        child_heavypath pointer).
        We call the entire tree the HeavyPathTree (HPT)

@note   Allocate a pointer to an array that will hold a path to the root for
        a leaf in the HPT, for use whe calculating add_leaf_HPT(). Each leaf
        of the HPT will have path_to_root_p set to this value, as well as
        the root of the HPT.
*/
Path* do_heavy_decomposition(Node *root, int k)
{
  int maxdepth = 0;
  Path*** path_to_root_pointer = malloc(sizeof(Path**)); // Pointer to the pointer vector (pointer to the pointer of the first Node).
  int* order_array = calloc(k, sizeof(int)); // The array to use for storing order of picking from several arrays.
  Path* heavy_path_tree_root = heavy_decomposition(root, 0, &maxdepth, k, 
                                                   path_to_root_pointer, order_array); // Returns pointer to the root of HPT.

  *path_to_root_pointer = calloc(maxdepth+1, sizeof(Path*)); //path to root: Pointer array of pointers.
  heavy_path_tree_root->path_to_root_p = path_to_root_pointer;
  //set_paths_to_root(heavy_path_tree_root);

  return heavy_path_tree_root;
}

Path* heavy_decomposition(Node *root, int depth, int *maxdepth, int k,
                          Path*** path_to_root_pointer, int* order_array)
{
  int length;
  Node** heavypath = get_heavypath(root, &length);

  Path *path_root;
  if(length == 1)
    path_root = heavypath_leaf(heavypath[0], depth, maxdepth, k,
                               path_to_root_pointer, order_array);
  else
    path_root = partition_heavypath(heavypath, length, depth, maxdepth, k,
                                    path_to_root_pointer, order_array);

  free(heavypath);

  return path_root;
}

/*
Free the memory for the HeavyPathTree (allocated in heavy_decomposition).
*/
void free_HPT(Path* root)
{
  free(*root->path_to_root_p); // Pointer to the first Node pointer.
  free(root->path_to_root_p); // Pointer to the pointer array of Node pointers
  free(root->order_array);
  free_HPT_rec(root);
}
void free_HPT_rec(Path* node)
{
  if(node->child_heavypaths)                    //PT leaf with descendents
  {
    for(int i=0; i < node->num_child_paths; i++)
      free_HPT_rec(node->child_heavypaths[i]);  //decend to next PT

    free(node->child_heavypaths);
  }

  else if(!node->node)                      //not leaf of HPT (internal PT node)
  {
    free_HPT_rec(node->left);
    free_HPT_rec(node->right);
  }

  //else                                    //a leaf of the HPT
  //{
  //  free(node->path_to_root);
  //}

  //Free dynamically added arrays for minimums and maximums.

  free(node->d_min_path);
  free(node->d_min_subtree);
  free(node->d_max_path);
  free(node->d_max_subtree);

  free(node->d_min_path_node);
  free(node->d_min_subtree_node);
  free(node->d_max_path_node);
  free(node->d_max_subtree_node);


  free(node);
}

/*
For the given heavypath, create a Path structure that represents the path.
Split the path in half and create a Path for each half.  If a half is a single
node, then hang the next heavy path off of it. If it's a leaf of alt_tree, then
link the Path to the corresponding leaf in alt_tree.
*/
Path* partition_heavypath(Node **heavypath, int length, int depth, int *maxdepth, 
                      int k, Path ***path_to_root_pointer, int* order_array)
{
  Path* newpath = new_Path(k);
  newpath->order_array = order_array;
  newpath->total_depth = depth;

    //Split the heavy path into two equal-length subpaths:
  int l1 = ceil(length/2);
  if(l1 == 1)
    newpath->left = heavypath_leaf(heavypath[0], depth+1, maxdepth, k,
                                   path_to_root_pointer, order_array);
  else
    newpath->left = partition_heavypath(heavypath, l1, depth+1, maxdepth, k,
                                        path_to_root_pointer, order_array);
  newpath->left->parent = newpath;

  int l2 = length - l1;
  if(l2 == 1)
    newpath->right = heavypath_leaf(heavypath[l1], depth+1, maxdepth, k,
                                    path_to_root_pointer, order_array);
  else
    newpath->right = partition_heavypath(&heavypath[l1], l2, depth+1, maxdepth, k,
                                         path_to_root_pointer, order_array);
  newpath->right->parent = newpath;

  newpath->right->sibling = newpath->left;
  newpath->left->sibling = newpath->right;

  /* update d_min_path and d_min_path_node */
  compare_two_sorted_ascending(newpath->left->d_min_path, newpath->right->d_min_path, 
        newpath->order_array, newpath->k, newpath->k, newpath-> k); //order_array will be updated.
  set_int_array_from_order_2(newpath->left->d_min_path, newpath->right->d_min_path, newpath->d_min_path,
        newpath->order_array, newpath->k); //update newpath->d_min_path using left and right.
  set_node_array_from_order_2(newpath->left->d_min_path_node, newpath->right->d_min_path_node, newpath->d_min_path_node,
        newpath->order_array, newpath->k); //update newpath->d_min_path_node using left and right.

  /* update d_max_path and d_max_path_node */
  compare_two_sorted_descending(newpath->left->d_max_path, newpath->right->d_max_path, 
        newpath->order_array, newpath->k, newpath->k, newpath-> k); //order_array will be updated.
  set_int_array_from_order_2(newpath->left->d_max_path, newpath->right->d_max_path, newpath->d_max_path,
        newpath->order_array, newpath->k); //update newpath->d_max_path using left and right.
  set_node_array_from_order_2(newpath->left->d_max_path_node, newpath->right->d_max_path_node, newpath->d_max_path_node,
        newpath->order_array, newpath->k); //update newpath->d_max_path_node using left and right.

  /* update d_min_subtree and d_min_path_subtree */
  compare_two_sorted_ascending(newpath->left->d_min_subtree, newpath->right->d_min_subtree, 
        newpath->order_array, newpath->k, newpath->k, newpath-> k); //order_array will be updated.
  set_int_array_from_order_2(newpath->left->d_min_subtree, newpath->right->d_min_subtree, newpath->d_min_subtree,
        newpath->order_array, newpath->k); //update newpath->d_min_subtree using left and right.
  set_node_array_from_order_2(newpath->left->d_min_subtree_node, newpath->right->d_min_subtree_node, newpath->d_min_subtree_node,
        newpath->order_array, newpath->k); //update newpath->d_min_subtree_node using left and right.


  /* update d_max_subtree and d_max_path_subtree */
  compare_two_sorted_descending(newpath->left->d_max_subtree, newpath->right->d_max_subtree, 
        newpath->order_array, newpath->k, newpath->k, newpath-> k); //order_array will be updated.
  set_int_array_from_order_2(newpath->left->d_max_subtree, newpath->right->d_max_subtree, newpath->d_max_subtree,
        newpath->order_array, newpath->k); //update newpath->d_max_subtree using left and right.
  set_node_array_from_order_2(newpath->left->d_max_subtree_node, newpath->right->d_max_subtree_node, newpath->d_max_subtree_node,
        newpath->order_array, newpath->k); //update newpath->d_max_subtree_node using left and right.


  // int min_path_arg = argmin(newpath->left->d_min_path, newpath->right->d_min_path);
  // if (min_path_arg == 0){
  //   newpath->d_min_path = newpath->left->d_min_path;
  //   newpath->d_min_path_node = newpath->left->d_min_path_node;
  // }else{
  //   newpath->d_min_path = newpath->right->d_min_path;
  //   newpath->d_min_path_node = newpath->right->d_min_path_node;
  // }

  // int max_path_arg = argmax(newpath->left->d_max_path, newpath->right->d_max_path);
  // if (max_path_arg == 0){
  //   newpath->d_max_path = newpath->left->d_max_path;
  //   newpath->d_max_path_node = newpath->left->d_max_path_node;
  // }else{
  //   newpath->d_max_path = newpath->right->d_max_path;
  //   newpath->d_max_path_node = newpath->right->d_max_path_node;
  // }

  // int max_subtree_arg = argmax(newpath->left->d_max_subtree,
  //                              newpath->right->d_max_subtree);
  // if (max_subtree_arg==0){
  //   newpath->d_max_subtree = newpath->left->d_max_subtree;
  //   newpath -> d_max_subtree_node = newpath -> left -> d_max_subtree_node;
  // }else{
  //   newpath->d_max_subtree = newpath->right->d_max_subtree;
  //   newpath -> d_max_subtree_node = newpath -> right -> d_max_subtree_node;
  // }

  // assert((newpath->left->d_min_subtree[0] == newpath->right->d_min_subtree[0]) && (newpath->left->d_min_subtree == 1));
  // newpath->d_min_subtree_node = newpath->left->d_min_subtree_node; // Choose the left one.       
                      
  return newpath;
}

/*
Return a Path for the given node of alt_tree.  The Path will be a leaf
node of the path tree.
Either 1) the leaf will point to a leaf node of alt_tree,
or     2) child_heavypath will point to a heavypath representing the
          descendant of the alt_tree node.
*/
Path* heavypath_leaf(Node *node, int depth, int *maxdepth, int k,
                     Path ***path_to_root_pointer, int* order_array)
{
  Path* newpath = new_Path(k);
  newpath->order_array = order_array;

  newpath->total_depth = depth;
  newpath->node = node;           //attach the path to the node
  node->path = newpath;           //attach the node to the path

  newpath->d_min_path[0] = newpath->d_max_path[0] = node->subtreesize;
  newpath->d_min_path_node[0] = newpath->d_max_path_node[0] = node;

    //Handle an internal alt_tree node with pendant heavypath:
  if(node->nneigh > 1)
  {
    newpath->num_child_paths = node->nneigh - 2;
    int i_neigh = 1;              //don't look at the parent
    if(node->depth == 0)          //If we are the root of alt_tree
    {
      i_neigh = 0;                //then we have no parent.
      newpath->num_child_paths += 1;
    }

    newpath->child_heavypaths = calloc(newpath->num_child_paths, sizeof(Path*));
    // newpath->d_min_subtree = INT_MAX;
    // newpath->d_max_subtree = INT_MIN;
    int* temporary_array = calloc(k, sizeof(int));
    Node** temporary_node_array = calloc(k, sizeof(Node*));

    int j = 0;                    //index the child heavypaths
    while(i_neigh < node->nneigh)
    {
      if(node->neigh[i_neigh] != node->heavychild)
      {
        newpath->child_heavypaths[j] = heavy_decomposition(node->neigh[i_neigh],
                                                           depth+1, maxdepth, k,
                                                           path_to_root_pointer, order_array);
        newpath->child_heavypaths[j]->parent_heavypath = newpath;
        

        /* set d_min_subtree and d_min_subtree_node*/
        compare_three_sorted_ascending(newpath->d_min_subtree, newpath->child_heavypaths[j]->d_min_path, newpath->child_heavypaths[j]->d_min_subtree, 
                newpath->order_array, k, k, k, k);
        set_int_array_from_order_3(newpath->d_min_subtree, newpath->child_heavypaths[j]->d_min_path, newpath->child_heavypaths[j]->d_min_subtree,
                temporary_array, newpath->order_array, k);
        update_int_array(temporary_array, newpath->d_min_subtree, 0, k);
        set_node_array_from_order_3(newpath->d_min_subtree_node, newpath->child_heavypaths[j]->d_min_path_node, newpath->child_heavypaths[j]->d_min_subtree_node,
                temporary_node_array, newpath->order_array, k);
        copy_node_array(temporary_node_array, newpath->d_min_subtree_node, k);
        
        
        /* set d_max_subtree and d_max_subtree_node*/
        compare_three_sorted_descending(newpath->d_max_subtree, newpath->child_heavypaths[j]->d_max_path, newpath->child_heavypaths[j]->d_max_subtree, 
                newpath->order_array, k, k, k, k);
        set_int_array_from_order_3(newpath->d_max_subtree, newpath->child_heavypaths[j]->d_max_path, newpath->child_heavypaths[j]->d_max_subtree,
                temporary_array, newpath->order_array, k);
        update_int_array(temporary_array, newpath->d_max_subtree, 0, k);
        set_node_array_from_order_3(newpath->d_max_subtree_node, newpath->child_heavypaths[j]->d_max_path_node, newpath->child_heavypaths[j]->d_max_subtree_node,
                temporary_node_array, newpath->order_array, k);
        copy_node_array(temporary_node_array, newpath->d_max_subtree_node, k);

        j++;
      }

      i_neigh++;
    }

    free(temporary_array);
    free(temporary_node_array);

    // assert(newpath->d_min_path = newpath->d_max_path);
    // assert(newpath->d_max_path = node->subtreesize);
  }
  else                            //HPT leaf (corresponds to alt_tree leaf)
  {
    *maxdepth = max(depth, *maxdepth);
    newpath->path_to_root_p = path_to_root_pointer;
    // assert(newpath->d_min_path == newpath->node->subtreesize);
    // assert(newpath->d_max_path == newpath->node->subtreesize);
  }

  return newpath;
}


/*
Descend to the leaves of the HPT. Once there, create the path_to_root vector
for that Path object.
*/
//void set_paths_to_root(Path* node)
//{
//  if(node->child_heavypaths)
//    for(int i=0; i < node->num_child_paths; i++)
//      set_paths_to_root(node->child_heavypaths[i]);
//
//  else if(node->left)
//  {
//    set_paths_to_root(node->left);
//    set_paths_to_root(node->right);
//  }
//
//  else                                //a leaf of the HPT tree
//    node->path_to_root = get_path_to_root_HPT(node);
//}

/*
Build a path (vector of Path*) from this Path leaf up to the root of the HPT,
following each PT to it's root in turn.
*/
void set_path_to_root_HPT(Path* leaf, Path** path_to_root)
{
  int i_path = 0;
  Path* w = leaf;
  while(w != NULL)                    //traverse up between PTs
  {
    while(1)                          //traverse up each PT
    {
      path_to_root[i_path++] = w;
      if(w->parent == NULL)
        break;

      w = w->parent;
    }

    w = w->parent_heavypath;
  }
}

/*
Build a path (vector of Path*) from this Path leaf up to the root of the HPT,
following each PT to it's root in turn.

@warning    user reponsible for memory
*/
Path** get_path_to_root_HPT(Path* leaf)
{
  int pathlen = leaf->total_depth+1;
  Path** path_to_root = calloc(pathlen, sizeof(Path*));
  int i_path = 0;

  Path* w = leaf;
  while(w != NULL)                    //traverse up between PTs
  {
    while(1)                          //traverse up each PT
    {
      path_to_root[i_path++] = w;
      if(w->parent == NULL)
        break;

      w = w->parent;
    }

    w = w->parent_heavypath;
  }
  assert(i_path == pathlen);

  return path_to_root;
}



/*
Return the heavypath rooted at the node.
   - user responsible for memory of returned heavypath.
*/
Node** get_heavypath(Node* root, int* length)
{
  *length = get_heavypath_length(root);

  Node** heavypath = calloc(*length, sizeof(Node*));
  Node* current = root;
  for(int i=0; i < *length; i++)
  {
    heavypath[i] = current;
    current = current->heavychild;
  }

  return heavypath;
}


/*
Return the number of nodes in the heavy path rooted at this node.
*/
int get_heavypath_length(Node *n)
{
    //Get the heavypath length:
  int length = 1;
  while(n->nneigh != 1)      //not a leaf
  {
    length++;
    n = n->heavychild;
  }
  return length;
}

/*
Return True if the (sub)Path corresponds to a leaf of alt_tree (i.e. it is a
leaf of the HPT).
*/
bool is_HPT_leaf(Path *n)
{
  return n->node && !n->child_heavypaths; //or (n->node && n->node->nneigh == 1)
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - Update Heavypath Tree - - - - - - - - - - - - - -


/*
Add the given leaf (from alt_tree) to the set L(v) for all v on a path from
leaf to the root.
*/
void add_leaf_HPT(Node* leaf)
{
  // fprintf(stderr, "leaf added %d\n", leaf->id);
  DB_TRACE(0, "alt_tree leaf - "); DB_CALL(0, print_node(leaf)); 
  // fprintf(stderr, "Leaf: %d\n", leaf->other->leaf_id);
    //Go down the path from the root to leaf, pushing down and modifying diff
    //values as we go. Subtract 1 for nodes on the path, and nodes above (to
    //the left) in the heavypath, and add 1 for nodes to the right in the
    //heavypath:
  set_path_to_root_HPT(leaf->path, *leaf->path->path_to_root_p);
  Path** path = *leaf->path->path_to_root_p;
  //Path** path = leaf->path->path_to_root;
  int pathlen = path[0]->total_depth + 1; //length (in number of nodes)
  for(int i = pathlen-1; i > 0; i--)      //go from root to parent of leaf.
  {
    // if (path[i]->id == 4) fprintf(stderr, "here %d\n", path[i]->num_child_paths);
    if(path[i]->node)                 //leaf of PT (points to node in alt_tree)
    {                                 //(child is root of heavypath)
      for(int j=0; j < path[i]->num_child_paths; j++)
      {
        path[i]->child_heavypaths[j]->diff_path += path[i]->diff_subtree;
        path[i]->child_heavypaths[j]->diff_subtree += path[i]->diff_subtree;

        if(path[i]->child_heavypaths[j] != path[i-1])
        {
          path[i]->child_heavypaths[j]->diff_path += 1;
          path[i]->child_heavypaths[j]->diff_subtree += 1;
        }
      }

      path[i]->d_min_path[0] += path[i]->diff_path - 1;
      path[i]->d_max_path[0] = path[i]->d_min_path[0];
    }
    else                              //internal node of PT
    {
      path[i-1]->diff_path += path[i]->diff_path;
      path[i-1]->diff_subtree += path[i]->diff_subtree;

      if(path[i-1] == path[i]->right) //right child of i is in the path
      {
        path[i]->left->diff_path += path[i]->diff_path - 1;
        path[i]->left->diff_subtree += path[i]->diff_subtree + 1;
      }
      else                            //left child of i is in the path
      {
        assert(path[i-1] == path[i]->left);
        path[i]->right->diff_path += path[i]->diff_path + 1;
        path[i]->right->diff_subtree += path[i]->diff_subtree + 1;
      }
    }

    path[i]->diff_path = path[i]->diff_subtree = 0;
  }
  assert(path[0]->node && path[0]->child_heavypaths == NULL &&
         path[0]->left == NULL && path[0]->right == NULL);       //HPT leaf
  path[0]->d_min_path[0] += path[0]->diff_path - 1;
  path[0]->d_max_path[0] = path[0]->d_min_path[0];
  path[0]->diff_path = path[0]->diff_subtree = 0;


    //Go up the path, updating the min and max values along the way:

  // Allocate memory for temporary array used for update
  int* temporary_array = calloc(path[0]->k, sizeof(int));
  Node** temporary_node_array = calloc(path[0]->k, sizeof(Node*));
  int* zero_array = calloc(path[0]->k, sizeof(int));
  Node** null_array = calloc(path[0]->k, sizeof(Node*));
  for(int i=1; i < pathlen; i++)
  {
    if(path[i]->child_heavypaths)               //leaf of a PT. (But not HPT leaf since we start from the index 1.)
    { 
      // fprintf(stderr, "%d\n", path[i]->num_child_paths);
      assert(path[i]->child_heavypaths[0]->diff_path == path[i]->child_heavypaths[0]->diff_subtree);
      //This is to confirm that for the root node of PT, diff_path and diff_subtree always matches.
                                                //d_min_path = d_max_path already set

      // Initialize d_min_subtree and d_min_subtree_node
      update_int_array(zero_array, path[i]->d_min_subtree, INT_MAX, path[i]->k);
      copy_node_array(null_array, path[i]->d_min_subtree_node, path[i]->k);

      // Initialize d_max_subtree and d_max_subtree_node
      update_int_array(zero_array, path[i]->d_max_subtree, INT_MIN, path[i]->k);
      copy_node_array(null_array, path[i]->d_max_subtree_node, path[i]->k);


      for(int j=0; j < path[i]->num_child_paths; j++)
      {
        if(is_HPT_leaf(path[i]->child_heavypaths[j])) 
        /*leaf of alt_tree as well. 
      -> d_min_subtree of path[i]->child_heavypath is meaningless so we don't use it for updating d_min_subtree of this path[i]..*/
        {
          // Update d_min_subtree
          compare_two_sorted_ascending_with_diffs(path[i]->d_min_subtree, 0, path[i]->child_heavypaths[j]->d_min_path, 
                    path[i]->child_heavypaths[j]->diff_path, path[i]->order_array, path[i]->k, path[i]->k, path[i]->k);
          set_int_array_from_order_2_with_diffs(path[i]->d_min_subtree, 0, path[i]->child_heavypaths[j]->d_min_path,
                    path[i]->child_heavypaths[j]->diff_path, temporary_array, path[i]->order_array, path[i]->k);
          update_int_array(temporary_array, path[i]->d_min_subtree, 0, path[i]->k);
          set_node_array_from_order_2(path[i]->d_min_subtree_node, path[i]->child_heavypaths[j]->d_min_path_node, temporary_node_array,
                    path[i]->order_array, path[i]->k);
          copy_node_array(temporary_node_array, path[i]->d_min_subtree_node, path[i]->k);

          // Update d_max_subtree
          compare_two_sorted_descending_with_diffs(path[i]->d_max_subtree, 0, path[i]->child_heavypaths[j]->d_max_path, 
                    path[i]->child_heavypaths[j]->diff_path, path[i]->order_array, path[i]->k, path[i]->k, path[i]->k);
          set_int_array_from_order_2_with_diffs(path[i]->d_max_subtree, 0, path[i]->child_heavypaths[j]->d_max_path,
                    path[i]->child_heavypaths[j]->diff_path, temporary_array, path[i]->order_array, path[i]->k);
          update_int_array(temporary_array, path[i]->d_max_subtree, 0, path[i]->k);
          // if (path[i]->d_max_subtree[0] > 10000) fprintf(stderr, "d_max_subtree_exceeded: %d\n", path[i]->d_max_subtree[0]); -> not here
          set_node_array_from_order_2(path[i]->d_max_subtree_node, path[i]->child_heavypaths[j]->d_max_path_node, temporary_node_array,
                    path[i]->order_array, path[i]->k);
          copy_node_array(temporary_node_array, path[i]->d_max_subtree_node, path[i]->k);
          // if (path[i]->id==4) fprintf(stderr, "d_max_subtree_updated %d\n", path[i]->d_max_subtree[0]);
          // if(path[i]->id==4) fprintf(stderr, "id4 at j=%d: %d, %d, %d\n", j, path[i]->d_max_subtree[0], path[i]->child_heavypaths[j]->d_max_path[0], path[i]->child_heavypaths[j]->diff_path);

          
        }
        else
        {
          // update d_min_subtree

          compare_three_sorted_ascending_with_diffs(path[i]->d_min_subtree, 0, 
                                                    path[i]->child_heavypaths[j]->d_min_path, path[i]->child_heavypaths[j]->diff_path,
                                                    path[i]->child_heavypaths[j]->d_min_subtree, path[i]->child_heavypaths[j]->diff_subtree,
                                                    path[i]->order_array, path[i]->k, path[i]->k,path[i]->k,path[i]->k);
          set_int_array_from_order_3_with_diffs(path[i]->d_min_subtree, 0,
                                                path[i]->child_heavypaths[j]->d_min_path, path[i]->child_heavypaths[j]->diff_path,
                                                path[i]->child_heavypaths[j]->d_min_subtree, path[i]->child_heavypaths[j]->diff_subtree,
                                                temporary_array, path[i]->order_array, path[i]->k);
          update_int_array(temporary_array, path[i]->d_min_subtree, 0, path[i]->k);
          set_node_array_from_order_3(path[i]->d_min_subtree_node, path[i]->child_heavypaths[j]->d_min_path_node, 
                                      path[i]->child_heavypaths[j]->d_min_subtree_node, temporary_node_array, path[i]->order_array, path[i]->k);
          copy_node_array(temporary_node_array, path[i]->d_min_subtree_node, path[i]->k);

          // update d_max_subtree
          compare_three_sorted_descending_with_diffs(path[i]->d_max_subtree, 0, 
                                                    path[i]->child_heavypaths[j]->d_max_path, path[i]->child_heavypaths[j]->diff_path,
                                                    path[i]->child_heavypaths[j]->d_max_subtree, path[i]->child_heavypaths[j]->diff_subtree,
                                                    path[i]->order_array, path[i]->k, path[i]->k,path[i]->k,path[i]->k);
          set_int_array_from_order_3_with_diffs(path[i]->d_max_subtree, 0,
                                                path[i]->child_heavypaths[j]->d_max_path, path[i]->child_heavypaths[j]->diff_path,
                                                path[i]->child_heavypaths[j]->d_max_subtree, path[i]->child_heavypaths[j]->diff_subtree,
                                                temporary_array, path[i]->order_array, path[i]->k);
          update_int_array(temporary_array, path[i]->d_max_subtree, 0, path[i]->k);
          // if (path[i]->d_max_subtree[0] > 10000 && path[i]->child_heavypaths[j]->d_max_subtree[0] < 10000) fprintf(stderr, "d_max_subtree_exceeded: %d \n", path[i]->child_heavypaths[j]->diff_subtree); //-> here.
          set_node_array_from_order_3(path[i]->d_max_subtree_node, path[i]->child_heavypaths[j]->d_max_path_node, 
                                      path[i]->child_heavypaths[j]->d_max_subtree_node, temporary_node_array, path[i]->order_array, path[i]->k);
          copy_node_array(temporary_node_array, path[i]->d_max_subtree_node, path[i]->k);

          // if (path[i]->id==4) fprintf(stderr, "d_max_subtree_updated %d\n", path[i]->d_max_subtree[0]);
        }
      }
    }
    else // internal node of a PT
    {
      assert(path[i]->left && path[i]->right);  //internal PT node

      // update d_min_path
      compare_two_sorted_ascending_with_diffs(path[i]->left->d_min_path, path[i]->left->diff_path,
                                              path[i]->right->d_min_path, path[i]->right->diff_path,
                                              path[i]->order_array, path[i]->k, path[i]->k, path[i]->k);
      set_int_array_from_order_2_with_diffs(path[i]->left->d_min_path, path[i]->left->diff_path,
                                            path[i]->right->d_min_path, path[i]->right->diff_path,  
                                            path[i]->d_min_path, path[i]->order_array, path[i]->k);
      set_node_array_from_order_2(path[i]->left->d_min_path_node, path[i]->right->d_min_path_node, path[i]->d_min_path_node, 
                                  path[i]->order_array, path[i]->k);
      

      // update d_max_path
      compare_two_sorted_descending_with_diffs(path[i]->left->d_max_path, path[i]->left->diff_path,
                                              path[i]->right->d_max_path, path[i]->right->diff_path,
                                              path[i]->order_array, path[i]->k, path[i]->k, path[i]->k);
      set_int_array_from_order_2_with_diffs(path[i]->left->d_max_path, path[i]->left->diff_path,
                                            path[i]->right->d_max_path, path[i]->right->diff_path,  
                                            path[i]->d_max_path, path[i]->order_array, path[i]->k);
      set_node_array_from_order_2(path[i]->left->d_max_path_node, path[i]->right->d_max_path_node, path[i]->d_max_path_node, 
                                  path[i]->order_array, path[i]->k);


      if(is_HPT_leaf(path[i]->left))        //left is leaf of alt_tree
      {
        assert(0); //THIS CASE SHOULD NOT HAPPEN.
        fprintf(stderr, "Internal error: Unexpected condition encountered. Exiting.");
        exit(0);
        // path[i]->d_min_subtree = path[i]->right->d_min_subtree +
        //                          path[i]->right->diff_subtree;
        // path[i]->d_max_subtree = path[i]->right->d_max_subtree +
        //                          path[i]->right->diff_subtree;
      }
      else if(is_HPT_leaf(path[i]->right))  //right is leaf of alt_tree
      {
        //update d_min_subtree
        update_int_array(path[i]->left->d_min_subtree, path[i]->d_min_subtree, path[i]->left->diff_subtree, path[i]->k);
        copy_node_array(path[i]->left->d_min_subtree_node, path[i]->d_min_subtree_node, path[i]->k);

        //update d_max_subtree
        update_int_array(path[i]->left->d_max_subtree, path[i]->d_max_subtree, path[i]->left->diff_subtree, path[i]->k);
        copy_node_array(path[i]->left->d_max_subtree_node, path[i]->d_max_subtree_node, path[i]->k);
        // if (path[i]->id==4) fprintf(stderr, "d_max_subtree_updated %d\n", path[i]->d_max_subtree[0]);
      }
      else
      {
        
        //update d_min_subtree
        compare_two_sorted_ascending_with_diffs(path[i]->left->d_min_subtree,path[i]->left->diff_subtree, 
                                                path[i]->right->d_min_subtree, path[i]->right->diff_subtree,
                                                path[i]->order_array, path[i]->k, path[i]->k, path[i]->k);
        set_int_array_from_order_2_with_diffs(path[i]->left->d_min_subtree,path[i]->left->diff_subtree, 
                                              path[i]->right->d_min_subtree, path[i]->right->diff_subtree, 
                                              path[i]->d_min_subtree, path[i]->order_array, path[i]->k);
        set_node_array_from_order_2(path[i]->left->d_min_subtree_node, path[i]->right->d_min_subtree_node, 
                                    path[i]->d_min_subtree_node, path[i]->order_array, path[i]->k);

        //update d_max_subtree
        compare_two_sorted_descending_with_diffs(path[i]->left->d_max_subtree, path[i]->left->diff_subtree, 
                                                path[i]->right->d_max_subtree, path[i]->right->diff_subtree,
                                                path[i]->order_array, path[i]->k, path[i]->k, path[i]->k);
        set_int_array_from_order_2_with_diffs(path[i]->left->d_max_subtree, path[i]->left->diff_subtree, 
                                              path[i]->right->d_max_subtree, path[i]->right->diff_subtree, 
                                              path[i]->d_max_subtree, path[i]->order_array, path[i]->k);
        set_node_array_from_order_2(path[i]->left->d_max_subtree_node, path[i]->right->d_max_subtree_node, 
                                    path[i]->d_max_subtree_node, path[i]->order_array, path[i]->k);        
        // if (path[i]->id==4) fprintf(stderr, "d_max_subtree_updated %d\n", path[i]->d_max_subtree[0]);             
      }
    }
    // if (path[i]->d_max_subtree[0]!=debug_arr[global_count]){
    //   fprintf(stderr, "%d\n",i);
    //   fprintf(stderr, "%d and %d \n", path[i]->d_max_subtree[0], debug_arr[global_count]);
    //   // if(is_HPT_leaf(path[i]->right)) fprintf(stderr, "Right is HPT leaf\n"); -> never!
    //   fprintf(stderr, "left_dmax_subtree + diff: %d + %d\n", path[i]->left->d_max_subtree[0] , path[i]->left->diff_subtree);
    //   fprintf(stderr, "right_dmax_subtree + diff: %d + %d\n", path[i]->right->d_max_subtree[0] , path[i]->right->diff_subtree);
    //   fprintf(stderr, "order array 0: %d\n", path[i]->order_array[0]);
    //   if(path[i]->left->child_heavypaths) fprintf(stderr, "Left is leaf %d\n",path[i]->left->id);
    // }
    // global_count += 1;
  }
  free(temporary_array);
  free(temporary_node_array);
  free(zero_array);
  free(null_array);
}



/* Return the root of the HPT with the given alt_tree leaf.
*/
Path* get_HPT_root(Node* leaf)
{
  assert(leaf->nneigh == 1);
  Path* node = leaf->path;
  Path* w = node;
  while(node != NULL)                 //traverse up between PTs
  {
    while(node->parent != NULL)       //traverse up each PT
      node = node->parent;

    w = node;
    node = node->parent_heavypath;
  }

  return w;
}



bool check_unique_ptr_array(Node** ptr_array, int k){
    Node **non_null_pointers = malloc(k * sizeof(void*));
    int non_null_count = 0;

    // Filter out NULL pointers and store the non-NULL pointers in non_null_pointers array
    for (int i = 0; i < k; i++) {
        if (ptr_array[i] != NULL) {
            non_null_pointers[non_null_count++] = ptr_array[i];
        }
    }

    // Check for duplicates among the non-NULL pointers
    for (int i = 0; i < non_null_count; i++) {
        for (int j = i + 1; j < non_null_count; j++) {
            if (non_null_pointers[i] == non_null_pointers[j]) {
                fprintf(stderr, "\n\n%d %d \n\n", non_null_pointers[i]->id, non_null_pointers[i]->id);
                free(non_null_pointers);
                return false;  // Duplicate found
            }
        }
    }

    free(non_null_pointers);
    return true;  // No duplicates found
}


/*
Reset the path and subtree min and max, along with the diff values for the
path from the given leaf to the root of the HPT.
*/
void reset_leaf_HPT(Node *leaf)
{
    //Follow the path from the leaft to the root, resetting the values along
    //the way:
    // Instead of doing the reverse transformation of add_leaf_HPT, 
    // It does simply try to reset all the values to the initial value
    // using the fact that diff's are only changed for those nodes in the HPT paths or sisters of the HPT paths,
    // and that min/max's are only changed for those nodes in the HPT paths.
  Path* w = leaf->path;
  // Path* lastw = w;

  // Allocate memory for temporary arrays
  int* temporary_array = calloc(w->k, sizeof(int));
  Node** temporary_node_array = calloc(w->k, sizeof(Node*));
  int* zero_array = calloc(w->k, sizeof(int));
  Node** null_array = calloc(w->k, sizeof(Node*));

  while(w != NULL)                    //traverse up between PTs
  {
    // Check if k minimums are all pointed to different node
    assert(check_unique_ptr_array(w->d_min_path_node, w->k));
    assert(check_unique_ptr_array(w->d_max_path_node, w->k));
    assert(check_unique_ptr_array(w->d_min_subtree_node, w->k));
    assert(check_unique_ptr_array(w->d_max_subtree_node, w->k));

    // w here is a leaf node of a PT tree.
    w->diff_path = w->diff_subtree = 0;
    w->d_min_path[0] = w->d_max_path[0] = w->node->subtreesize; // d_min(max)_path[1:] should be still kept as INT_MAX(MIN).
    assert(w->d_min_path_node[0] == w->node);
    // printf("%i\n", w->d_min_path[1]);
    assert(w->k==1 ||  w->d_min_path[1]==INT_MAX);
    assert(w->k==1 || w->d_max_path[1]==INT_MIN); 
    if(!is_HPT_leaf(w))               //this PT leaf is not a HPT leaf; corresponding to an internal node in the alt_tree.
    { 
      // Initialize d_min_subtree  
      update_int_array(zero_array, w->d_min_subtree, INT_MAX, w->k);
      copy_node_array(null_array, w->d_min_subtree_node, w->k);

      // Initialize d_max_subtree
      update_int_array(zero_array, w->d_max_subtree, INT_MIN, w->k);
      copy_node_array(null_array, w->d_max_subtree_node, w->k);


      // compare_two_sorted_ascending(w->child_heavypaths[0]->d_min_path, w->child_heavypaths[0]->d_min_subtree, w->order_array, w->k, w->k, w->k);
      // set_int_array_from_order_2(w->child_heavypaths[0]->d_min_path, w->child_heavypaths[0]->d_min_subtree, w->d_min_subtree, w->order_array, w->k);
      // set_node_array_from_order_2(w->child_heavypaths[0]->d_min_path_node, w->child_heavypaths[0]->d_min_subtree_node, w->d_min_subtree_node, w->order_array, w->k);
    
      // Initialize d_max_subtree
      // compare_two_sorted_descending(w->child_heavypaths[0]->d_max_path, w->child_heavypaths[0]->d_max_subtree, w->order_array, w->k, w->k, w->k);
      // set_int_array_from_order_2(w->child_heavypaths[0]->d_max_path, w->child_heavypaths[0]->d_max_subtree, w->d_max_subtree, w->order_array, w->k);
      // set_node_array_from_order_2(w->child_heavypaths[0]->d_max_path_node, w->child_heavypaths[0]->d_max_subtree_node, w->d_max_subtree_node, w->order_array, w->k);

      // w->child_heavypaths[0]->diff_path = 0;
      // w->child_heavypaths[0]->diff_subtree = 0;


      for(int j=0; j < w->num_child_paths; ++j){
        if(is_HPT_leaf(w->child_heavypaths[j])) 
        {
          // Update d_min_subtree
          compare_two_sorted_ascending_with_diffs(w->d_min_subtree, 0, w->child_heavypaths[j]->d_min_path, 
                    0, w->order_array, w->k, w->k, w->k);
          set_int_array_from_order_2_with_diffs(w->d_min_subtree, 0, w->child_heavypaths[j]->d_min_path,
                    0, temporary_array, w->order_array, w->k);
          update_int_array(temporary_array, w->d_min_subtree, 0, w->k);
          set_node_array_from_order_2(w->d_min_subtree_node, w->child_heavypaths[j]->d_min_path_node, temporary_node_array,
                    w->order_array, w->k);
          copy_node_array(temporary_node_array, w->d_min_subtree_node, w->k);

          // Update d_max_subtree
          compare_two_sorted_descending_with_diffs(w->d_max_subtree, 0, w->child_heavypaths[j]->d_max_path, 
                    0, w->order_array, w->k, w->k, w->k);
          set_int_array_from_order_2_with_diffs(w->d_max_subtree, 0, w->child_heavypaths[j]->d_max_path,
                    0, temporary_array, w->order_array, w->k);
          update_int_array(temporary_array, w->d_max_subtree, 0, w->k);
          set_node_array_from_order_2(w->d_max_subtree_node, w->child_heavypaths[j]->d_max_path_node, temporary_node_array,
                    w->order_array, w->k);
          copy_node_array(temporary_node_array, w->d_max_subtree_node, w->k);
        }else
        {
          // update d_min_subtree
          compare_three_sorted_ascending_with_diffs(w->d_min_subtree, 0, 
                                                    w->child_heavypaths[j]->d_min_path, 0,
                                                    w->child_heavypaths[j]->d_min_subtree, 0,
                                                    w->order_array, w->k, w->k,w->k,w->k);
          set_int_array_from_order_3_with_diffs(w->d_min_subtree, 0,
                                                w->child_heavypaths[j]->d_min_path, 0,
                                                w->child_heavypaths[j]->d_min_subtree, 0,
                                                temporary_array, w->order_array, w->k);
          update_int_array(temporary_array, w->d_min_subtree, 0, w->k);
          set_node_array_from_order_3(w->d_min_subtree_node, w->child_heavypaths[j]->d_min_path_node, 
                                      w->child_heavypaths[j]->d_min_subtree_node, temporary_node_array, w->order_array, w->k);
          copy_node_array(temporary_node_array, w->d_min_subtree_node, w->k);

          // update d_max_subtree
          compare_three_sorted_descending_with_diffs(w->d_max_subtree, 0, 
                                                    w->child_heavypaths[j]->d_max_path, 0,
                                                    w->child_heavypaths[j]->d_max_subtree, 0,
                                                    w->order_array, w->k, w->k,w->k,w->k);
          set_int_array_from_order_3_with_diffs(w->d_max_subtree, 0,
                                                w->child_heavypaths[j]->d_max_path, 0,
                                                w->child_heavypaths[j]->d_max_subtree, 0,
                                                temporary_array, w->order_array, w->k);
          update_int_array(temporary_array, w->d_max_subtree, 0, w->k);
          set_node_array_from_order_3(w->d_max_subtree_node, w->child_heavypaths[j]->d_max_path_node, 
                                      w->child_heavypaths[j]->d_max_subtree_node, temporary_node_array, w->order_array, w->k);
          copy_node_array(temporary_node_array, w->d_max_subtree_node, w->k);
        }
        w->child_heavypaths[j]->diff_path = 0;
        w->child_heavypaths[j]->diff_subtree = 0;
      }
      // if (w->id == 4) fprintf(stderr, "id 4 reset: d_max_subtree: %d\n", w->d_max_subtree[0]);
    }

    while(w->parent != NULL)          //traverse up each PT
    {
      w = w->parent;

      w->diff_path = w->diff_subtree = 0;

      // update d_min_path
      compare_two_sorted_ascending(w->left->d_min_path, w->right->d_min_path, w->order_array, w->k, w->k, w->k);
      set_int_array_from_order_2(w->left->d_min_path, w->right->d_min_path, w->d_min_path, w->order_array, w->k);
      set_node_array_from_order_2(w->left->d_min_path_node, w->right->d_min_path_node, w->d_min_path_node, w->order_array, w->k);


      // update d_max_path
      compare_two_sorted_descending(w->left->d_max_path, w->right->d_max_path, w->order_array, w->k, w->k, w->k);
      set_int_array_from_order_2(w->left->d_max_path, w->right->d_max_path, w->d_max_path, w->order_array, w->k);
      set_node_array_from_order_2(w->left->d_max_path_node, w->right->d_max_path_node, w->d_max_path_node, w->order_array, w->k);


      // update d_min_subtree
      compare_two_sorted_ascending(w->left->d_min_subtree, w->right->d_min_subtree, w->order_array, w->k, w->k, w->k);
      set_int_array_from_order_2(w->left->d_min_subtree, w->right->d_min_subtree, w->d_min_subtree, w->order_array, w->k);
      set_node_array_from_order_2(w->left->d_min_subtree_node, w->right->d_min_subtree_node, w->d_min_subtree_node, w->order_array, w->k);

      // update d_max_subtree
      compare_two_sorted_descending(w->left->d_max_subtree, w->right->d_max_subtree, w->order_array, w->k, w->k, w->k);
      set_int_array_from_order_2(w->left->d_max_subtree, w->right->d_max_subtree, w->d_max_subtree, w->order_array, w->k);
      set_node_array_from_order_2(w->left->d_max_subtree_node, w->right->d_max_subtree_node, w->d_max_subtree_node, w->order_array, w->k);

      // int min_path_arg = argmin(w->left->d_min_path, w->right->d_min_path);
      // if (min_path_arg==0){
      //   w->d_min_path = w->left->d_min_path;
      //   w->d_min_path_node = w->left->d_min_path_node;
      // }else{
      //   w->d_min_path = w->right->d_min_path;
      //   w->d_min_path_node = w->right->d_min_path_node;
      // }
      // // w->d_min_path = min(w->left->d_min_path, w->right->d_min_path);
      // int max_path_arg = argmax(w->left->d_max_path, w->right->d_max_path);
      // if (max_path_arg==0){
      //   w->d_max_path = w->left->d_max_path;
      //   w->d_max_path_node = w->left->d_max_path_node;
      // }else{
      //   w->d_max_path = w->right->d_max_path;
      //   w->d_max_path_node = w->right->d_max_path_node;
      // }
      // w->d_max_path = max(w->left->d_max_path, w->right->d_max_path);


      // w->d_min_subtree = 1;
      // w->d_min_subtree_node = w->left->d_min_subtree_node;

      // int max_subtree_arg = argmax(w->left->d_max_subtree, w->right->d_max_subtree);
      // if (max_subtree_arg==0){
      //   w->d_max_subtree = w->left->d_max_subtree;
      //   w->d_max_subtree_node = w->left->d_max_subtree_node;
      // }else{
      //   w->d_max_subtree = w->right->d_max_subtree;
      //   w->d_max_subtree_node = w->right->d_max_subtree_node;
      // }
      // w->d_max_subtree = max(w->left->d_max_subtree, w->right->d_max_subtree);

      w->left->diff_path = w->left->diff_subtree = 0;
      w->right->diff_path = w->right->diff_subtree = 0;
    }

    // lastw = w;
    w = w->parent_heavypath;
  }

  free(temporary_array);
  free(temporary_node_array);
  free(zero_array);
  free(null_array);

}




// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - I/O - - - - - - - - - - - - - - - - - - - - -

/*
Print the given heavypath.
*/
void print_heavypath(Node **heavypath, int length)
{
  for(int i=0; i < length; i++)
    printf("%i ", heavypath[i]->id);
  printf("\n");
}

/*
Print the given Path node.
*/
void print_HPT_node(const Path* n)
{
  char *name = "----";
  if(n->node && n->node->nneigh == 1)  //a leaf
    name = n->node->name;
  fprintf(stderr, "node id: %i name: %s\n", n->id, name);
}


/*
Print the Heavy Path Tree (HPT) in dot format. HPT edges will be dashed, while
edges of alt_tree (not in the HPT) will be solid.
*/
void print_HPT_dot(Path* hproot, Node* altroot, int repid)
{
  char filename[15];
  sprintf(filename, "hptree_%i.dot", repid);
  fprintf(stderr, "Creating DOT: %s\n", filename);
  FILE *f = fopen(filename, "w");
  if(f == NULL)
  {
    fprintf(stderr, "Can't open file for writing!");
    exit(0);
  }

  fprintf(f, "digraph HPT\n  {\n  center=true;\n");
  print_HPT_keynode_dot(f);
  print_HPT_subpath_dot(hproot, f);
  print_HPT_subtree_dot(altroot, f);
  fprintf(f, "  }\n");
  fclose(f);
}

/*
Recursively print the subPath (for the tree structure on the Path).
*/
void print_HPT_subpath_dot(Path* n, FILE *f)
{
  if(n->left)
  {
    print_HPT_ptnode_dot(n, f);       //print node formatting information.

    fprintf(f, "  ");
    print_HPT_node_dot(n, f);
    fprintf(f, " -> ");
    print_HPT_node_dot(n->left, f);
    fprintf(f, " [style=dashed];\n");

    print_HPT_subpath_dot(n->left, f);
  }
  if(n->right)
  {
    //print_HPT_ptnode_dot(n, f);     //already printed this for n->left.

    fprintf(f, "  ");
    print_HPT_node_dot(n, f);
    fprintf(f, " -> ");
    print_HPT_node_dot(n->right, f);
    fprintf(f, " [style=dashed];\n");

    print_HPT_subpath_dot(n->right, f);
  }
  if(n->node)                         //node of alt_tree
  {
    print_HPT_hpnode_dot(n, f);       //print node formatting information.

    if(n->child_heavypaths)           //not a leaf of HPT (and alt_tree)
    {
      for(int i=0; i < n->num_child_paths; i++)
      {
        fprintf(f, "  ");
        print_HPT_node_dot(n, f);
        fprintf(f, " -> ");
        print_HPT_node_dot(n->child_heavypaths[i], f);
        fprintf(f, " [style=dashed color=gray];\n");

        print_HPT_subpath_dot(n->child_heavypaths[i], f);
      }
    }
  }
}


/*
Recursively print the subtree (for alt_tree edges).
*/
void print_HPT_subtree_dot(Node* node, FILE *f)
{
  if(node->nneigh > 1)    //not leaf
  {
    int firstchild = 1;
    if(node->depth == 0)
      firstchild = 0;

    for(int i = firstchild; i < node->nneigh; i++)
    {
      fprintf(f, "  ");
      print_HPT_node_dot(node->path, f);
      fprintf(f, " -> ");
      print_HPT_node_dot(node->neigh[i]->path, f);
      fprintf(f, ";\n");

      print_HPT_subtree_dot(node->neigh[i], f);
    }
  }
}


/*
Print a string representing this Path node formatted for dot output.
*/
void print_HPT_node_dot(Path* n, FILE *f)
{
  fprintf(f, "%i", n->id);
}

/*
Print a string that formats a heavypath (alt_tree) node in a PT (PathTree).

A heavypath node will be circular and have the following format:
  1 (2)                  id (alt_id)
p: 0 1 1    diff_path    d_min_path    d_max_path
s: 0 1 1    diff_subtree d_min_subtree d_max_subtree

if it's an internal node of alt_tree, and:

 3 (2):a              id (alt_id): name
P: 0 1 1    diff_path    d_min_path    d_max_path
s: 0 1 1    diff_subtree d_min_subtree d_max_subtree

if it's a leaf of alt_tree.
*/
void print_HPT_hpnode_dot(Path* n, FILE *f)
{
  if(n->node->nneigh == 1)              //a leaf of alt_tree
    fprintf(f, "  %i [label=\"%i (%i): %s",
            n->id, n->id, n->node->id, n->node->name);
  else
    fprintf(f, "  %i [label=\"%i (%i)", n->id, n->id, n->node->id);

  fprintf(f, "\np: %i\nd_min_path:", n->diff_path);
  for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_min_path[i]);
  fprintf(f, "\nd_max_path: ");
  for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_max_path[i]);

  // fprintf(f, "\np: %i %i %i\n", n->diff_path, n->d_min_path, n->d_max_path);

  if(n->node->nneigh == 1)              //a leaf of alt_tree
    fprintf(f, "s: %i x x", n->diff_subtree);
  else{
    fprintf(f, "s: %i x x\nd_min_subtree: ", n->diff_subtree);
    for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_min_subtree[i]);
    fprintf(f, "\nd_max_subtree: ");
    for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_max_subtree[i]);
  }
    // fprintf(f, "s: %i %i %i", n->diff_subtree, n->d_min_subtree, n->d_max_subtree);

  fprintf(f, "\"];\n");
}

/* Print a string that formats a PT (PathTree) node.

A pathtree node will be rectangular and have the following format:
     1                      id
p: 0 1 1    diff_path    d_min_path    d_max_path
s: 0 1 1    diff_subtree d_min_subtree d_max_subtree
*/
void print_HPT_ptnode_dot(Path* n, FILE *f)
{
  fprintf(f, "  %i [shape=rectangle ", n->id);
  fprintf(f, "label=\"%i", n->id);

  fprintf(f, "\np: %i %i\nd_min_path:", n->diff_path, n->diff_subtree);
  for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_min_path[i]);
  fprintf(f, "\nd_max_path");
  for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_max_path[i]);
  fprintf(f, "\nd_min_subtree");
  for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_min_subtree[i]);
  fprintf(f, "\nd_max_subtree");
  for(int i=0;i<n->k;i++) fprintf(f, "%i ", n->d_max_subtree[i]);
  


  // fprintf(f, "\np: %i %i %i\ns: %i %i %i", n->diff_path, n->d_min_path,
  //         n->d_max_path, n->diff_subtree, n->d_min_subtree, n->d_max_subtree);
  fprintf(f, "\"];\n");
}

/* Print a node that describes the values in the positions.
*/
void print_HPT_keynode_dot(FILE *f)
{
  fprintf(f, "  keynode [shape=record ");
  fprintf(f, "label=\"{node id|{");
  fprintf(f, "{diff_path|diff_subtree} | ");
  fprintf(f, "{min_path|min_subtree} | ");
  fprintf(f, "{max_path|max_subtree}}}");
  fprintf(f, "\"];\n");
}
