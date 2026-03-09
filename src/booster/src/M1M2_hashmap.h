#ifndef __M1M2_HASHMAP_H__
#define __M1M2_HASHMAP_H__



// セット構造体
typedef struct __IntSet{
    int is_external; //If it corresponds to an external branch.
    int input_count;//bipartitionを持つtreeの数
    int isinitial;//0(not have) or 1(have)
    int treeid;//input_count分あるtreeidのうちの1つ,initial treeしかなかったらtreeidは-1
    int bipartition_id; // ID of bipartition; Consecutive numbers are set during the adding process.
    int input_internal_bipartition_id; // Id within input trees; Consecutive numbers are set during the adding process; If only present in the initial, it is set to -1.
    int ref_id; // Corresponding node id in the reference tree; Set to -1 if it does not exist in the ref tree.
    unsigned int h2; //mod by M2
    struct __IntSet *next;
} IntSet;

typedef struct __M1M2Hashtable{
    IntSet** intset; //hash elements
    int input_unique_internal_bipar_count; //Number of unique bipartitions in the input trees. Should be equal to (maximum of input_bipartition_id -1)
    int unique_bipar_count; //Number of unique bipartitions in the whole set of trees.
} M1M2Hashtable;


typedef struct __simpleHash{
    unsigned int value;
    struct __simpleHash* next;
} simpleHash;


typedef struct __simpleHashTable{
    simpleHash** simple_hash;
    unsigned int table_length;
} simpleHashTable;

M1M2Hashtable* new_M1M2hashtable();
void add_M1M2hashtable(M1M2Hashtable* M1M2hash,unsigned int h1, unsigned int h2,int treeid,int isinitial, int is_external, int ref_id);
void free_M1M2hashtable(M1M2Hashtable* M1M2hash);
IntSet* find_element_M1M2hashtable(M1M2Hashtable* M1M2Hash, unsigned int h1, unsigned int h2);

void print_hashtable(M1M2Hashtable* M1M2hash);


simpleHashTable* new_simple_hashtable(unsigned int n_elements);

void Knuth_random_gen(unsigned int n_elements, unsigned int N, unsigned int* vector);
void Floyd_random_gen(unsigned int n_elements, unsigned int N, unsigned int* vector);


/* Global Variables*/
extern unsigned int* HASH_INTEGERS1;
extern unsigned int* HASH_INTEGERS2;
extern const int M1;
extern const int M2;

#endif