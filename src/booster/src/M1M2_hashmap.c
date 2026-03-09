#include "M1M2_hashmap.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#define MAX_UINT32 4294967295U




const int M1 = 524287;  /*prime*/
const int M2 = 998244353;   /*prime*/


/*ここから使い方
struct IntSet** M1M2hash;

 if((M1M2hash = (struct IntSet **) malloc(M1*sizeof(*M1M2hash))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 
ここまで*/


unsigned int* HASH_INTEGERS1 = NULL;
unsigned int* HASH_INTEGERS2 = NULL;

// This function allocates memory for the hashtable. Freeing outside the function is necessary.
M1M2Hashtable* new_M1M2hashtable(){
  M1M2Hashtable* hashtable = malloc(sizeof(M1M2Hashtable));
  hashtable->input_unique_internal_bipar_count = 0;
  hashtable->unique_bipar_count = 0;
  if((hashtable->intset = malloc(M1*sizeof(IntSet*))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 
  for(int h = 0 ; h < M1; h ++) hashtable->intset[h] = NULL;   
  return hashtable; 
}

// This function allocates memory for the hashtable. Freeing outside the function is necessary.
void add_M1M2hashtable(M1M2Hashtable* M1M2hash,unsigned int h1, unsigned int h2,int treeid, int isinitial, int is_external, int ref_id){//mod by M1,M2
  IntSet* tmp;
  IntSet* tmp_origin = NULL;

  assert(h1<M1);

  int debug=0;
  tmp = M1M2hash->intset[h1];

  bool newsplit = true;
  while(tmp != NULL){
    debug = 1;
      if(tmp ->h2 == h2){
        newsplit = false;
        break;
      }
      tmp_origin = tmp;
      tmp = tmp->next;    
  }

  
  if(newsplit == true){
    tmp = (IntSet *) malloc(sizeof(IntSet));

    tmp->h2 = h2;
    tmp->input_count = 0;
    tmp -> treeid = -1;//initial treeしかなかったらtreeidは-1.
    tmp->bipartition_id = M1M2hash->unique_bipar_count; 
    tmp->input_internal_bipartition_id = -1;//initial treeしかなかったらinput_bipartition_idは-1.
    tmp -> isinitial = 0;
    tmp->ref_id = -1;
    tmp->next = NULL;
    tmp->is_external = is_external;  

    M1M2hash->unique_bipar_count += 1;

    if(tmp_origin == NULL) M1M2hash->intset[h1] = tmp;//
    else{
      tmp_origin -> next = tmp;
    }
  }

  if(isinitial) {
    tmp -> isinitial = 1;
    tmp -> ref_id = ref_id;
  }
  else{
    if (tmp->input_count == 0 && !is_external) {
      tmp->input_internal_bipartition_id = M1M2hash->input_unique_internal_bipar_count;
      M1M2hash->input_unique_internal_bipar_count += 1; // This happens when ( (!isinitial) && ( (newsplit==true) || (newsplit==false && tmp->input_count == 0) ) )
    }
    tmp -> input_count ++;
    tmp -> treeid = treeid;
  }
}


IntSet* find_element_M1M2hashtable(M1M2Hashtable* M1M2Hash, unsigned int h1, unsigned int h2){
  if (h1 >= M1 || h2 >= M2){
    fprintf(stderr, "Hash value too big. Returning NULL pointer.\n");
    return NULL;
  }

  IntSet* hash_loc = M1M2Hash->intset[h1];
  while (hash_loc!=NULL){
    if (hash_loc->h2 == h2) break;
    hash_loc = hash_loc->next;
  }
  if (hash_loc == NULL) fprintf(stderr, "Hash value could not be found.\n");
  return hash_loc;
}

void free_M1M2hashtable(M1M2Hashtable* M1M2hash){
    for(int i = 0; i < M1; i++){
        IntSet *current = M1M2hash->intset[i];
        while(current != NULL){
            IntSet *next = current->next;
            free(current);
            current = next;
        }
    }
    free(M1M2hash->intset);
    free(M1M2hash);
}


simpleHashTable* new_simple_hashtable(unsigned int n_elements){
  simpleHashTable* simple_hash_table = malloc(sizeof(simpleHashTable));
  simple_hash_table->table_length = n_elements;
  simple_hash_table->simple_hash = malloc(n_elements * sizeof(simpleHash*));
  for (int i=0; i<n_elements; i++) simple_hash_table->simple_hash[i]=NULL;
  return simple_hash_table;
}

void add_new_simple_hash(unsigned int value, simpleHashTable* simple_hashtable) {
    unsigned int hash = value % simple_hashtable->table_length;
    
    // Allocate memory for the new node
    simpleHash* new_hash = malloc(sizeof(simpleHash));
    if (new_hash == NULL) {
        // Handle allocation failure (e.g., log error)
        return; // Or handle appropriately
    }
    
    new_hash->value = value;
    new_hash->next = NULL;

    // Check if the bucket is empty
    if (simple_hashtable->simple_hash[hash] == NULL) {
        // If empty, set the new node as the first entry
        simple_hashtable->simple_hash[hash] = new_hash;
    } else {
        // Otherwise, traverse to the end of the list
        simpleHash* current = simple_hashtable->simple_hash[hash];
        while (current->next != NULL) {
            current = current->next;
        }
        // Link the new node at the end
        current->next = new_hash;
    }
}

int exist_value(unsigned int value, simpleHashTable* simple_hashtable){
  unsigned int hash = value % simple_hashtable->table_length;
  simpleHash* hash_loc = simple_hashtable->simple_hash[hash];
  int exist = 0;
  while (hash_loc !=NULL){
    if (hash_loc->value == value){
      exist = 1;
      break;
    }
    hash_loc = hash_loc->next;
  }
  return exist;
}

void free_simple_hashtable(simpleHashTable* simple_hashtable){
  for (int i=0; i<simple_hashtable->table_length; i++){
    simpleHash* tmp = simple_hashtable->simple_hash[i];
    while (tmp!=NULL) {
      simpleHash* to_free = tmp;
      tmp = tmp->next;
      free(to_free);
    }
  }
  free(simple_hashtable->simple_hash);
  free(simple_hashtable);
}


unsigned int generate_random_32bit() {
    return ((unsigned int)rand() << 16) | rand();
}


/* Generates random array of size `n_elements` with unique elements in range [0,N-1] using Knuth's algorithm*/
void Knuth_random_gen(unsigned int n_elements, unsigned int N, unsigned int* vector){

  unsigned int m=0;

  for (unsigned int n = 0; n < N && m < n_elements; n++) {
    unsigned int rn = N - n;
    unsigned int rm = n_elements - m;
    if (generate_random_32bit() % rn < rm)    
      vector[m++] = n + 1;
  }
}

void Floyd_random_gen(unsigned int n_elements, unsigned int N, unsigned int* vector){

  simpleHashTable* SHT = new_simple_hashtable(10007);

  if (N > MAX_UINT32){
    fprintf(stderr, "The range of random number is too large.\n");
    exit(1);
  }
  unsigned int in, im;

  im = 0;

  for (in = N - n_elements; in < N && im < n_elements; ++in) {
    unsigned int r = generate_random_32bit() % (in + 1); /* generate a random number 'r' */

    if (exist_value(r, SHT))
      /* we already have 'r' */
      r = in; /* use 'in' instead of the generated number */

    assert(!exist_value(r, SHT));
    vector[im++] = r;
    add_new_simple_hash(r, SHT);
  }

  assert(im == n_elements);
  free_simple_hashtable(SHT);
}


void print_hashtable(M1M2Hashtable* M1M2hash){
  int one_count = 0;
  int multiple_count = 0;
  int double_count = 0;
  for (unsigned int i=0; i<M1; i++){
    IntSet* elem = M1M2hash->intset[i];
    if (elem!=NULL){
      int count = 0;
      while(elem!=NULL){
        elem = elem->next;
        count ++;
      }
      // fprintf(stderr, "%d elements found at h1=%u\n", count, i);
      if (count==1) one_count +=1;
      else if (count==2) double_count +=1;
      else if (count>2) multiple_count += 1;
    }
  }
  fprintf(stderr, "Uniquely hashed elements: %d\n", one_count);
  fprintf(stderr, "Slots that have two elements: %d\n", double_count);
  fprintf(stderr, "Slots that have two elements: %d\n", multiple_count);
  fprintf(stderr, "Number of unique elements in the hashtable: %d\n", M1M2hash->unique_bipar_count);
  fprintf(stderr, "Number of unique internal elements of the input trees in the hashtable: %d\n", M1M2hash->input_unique_internal_bipar_count);
}