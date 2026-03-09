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

#ifndef _EDGE_STACK_H_
#define _EDGE_STACK_H_

#include "tree.h"

/* The NodeStack elements are constituted of Node/Edge pairs */
typedef struct __NodeStackElt NodeStackElt;
typedef struct __NodeStackElt {
  Node *node;
  Edge *edge;
  struct __NodeStackElt *prev;
} NodeStackElt;


typedef struct __NodeStack NodeStack;
typedef struct __NodeStack {
  struct __NodeStackElt *head;
} NodeStack;


/* Initialize a new Node Stack*/
NodeStack *new_nodestack();
/* Pushes a new node/edge pair to the Stack */
void nodestack_push(NodeStack *ns, Node *n, Edge *e);
/* Pops and returns the head of the Stack. The calling function is
 responsible for freeing the elt: free(elt).*/
NodeStackElt * nodestack_pop(NodeStack *ns);
/* frees the whole stack and all its elements */
void nodestack_free(NodeStack *ns);

#endif
