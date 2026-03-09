#ifndef _SPLIT_H_
#define _SPLIT_H_


#include <limits.h>
#include <assert.h>
#include "debug.h"


typedef unsigned int* split;

unsigned int split_length(unsigned int n_tips);
int split_match(split a, split b, unsigned int spl);


void set_bit(split bitvector, unsigned int tip_id);
void bitwise_or(split bitvector1, split bitvector2, split result, unsigned int split_length);
void print_bitvector(split bitvector, unsigned int n_tips);


#endif