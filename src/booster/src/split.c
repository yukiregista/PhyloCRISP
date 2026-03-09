
#include "split.h"

unsigned int split_length(unsigned int n_tips)
{
    assert(n_tips>0);
    return (n_tips-1)/(sizeof(unsigned int) * CHAR_BIT) + 1;
}


int split_match(split a, split b, unsigned int spl){
    int res = 1;
    for (int i=0; i<spl; i++){
        if (a[i]!=b[i]){
            res=0;
            break;
        }
    }
    return res;
}

/*
Set the bit at tip_id 1.
*/
void set_bit(split bitvector, unsigned int tip_id)
{
    unsigned int bits_per_int = CHAR_BIT * sizeof(unsigned int);
    unsigned int index = tip_id / bits_per_int;
    unsigned int bit_position = tip_id % bits_per_int;
    bitvector[index] |= (1U << bit_position);
}

/*
Take bitwise or
*/
void bitwise_or(split bitvector1, split bitvector2, split result, unsigned int split_length)
{
    for(int i=0; i<split_length; ++i) result[i] = bitvector1[i] | bitvector2[i];
}

// Function to print the bit vector for debugging purposes
void print_bitvector(split bitvector, unsigned int n_tips) {
    unsigned int bits_per_int = sizeof(unsigned int) * CHAR_BIT;
    for (unsigned int i = 0; i < n_tips; i++) {
        unsigned int index = i / bits_per_int;
        unsigned int bit_position = i % bits_per_int;
        if (bitvector[index] & (1U << bit_position))
            printf("1");
        else
            printf("0");
        if ((i + 1) % bits_per_int == 0)
            printf(" ");
    }
    printf("\n");
}