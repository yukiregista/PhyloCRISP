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

#include "node_stack.h"

/*
This simple implementation of a stack is used for iterative parsing of 
Newick format
 */

NodeStack *new_nodestack(){
  NodeStack * ns = malloc(sizeof(NodeStack));
  ns->head = NULL;
  return(ns);
}

void nodestack_push(NodeStack *ns, Node *n, Edge *e){
  NodeStackElt *new_elt = malloc(sizeof(NodeStackElt));
  new_elt->node = n;
  new_elt->edge = e;
  new_elt->prev = ns->head;
  ns->head = new_elt;
}

NodeStackElt * nodestack_pop(NodeStack *ns){
  NodeStackElt *head = ns->head;
  if(head != NULL){
    ns->head = head->prev;
  }
  return head;
}

void nodestack_free(NodeStack *ns){
  NodeStackElt *elt = nodestack_pop(ns);
  while(elt!=NULL){
    free(elt);
    elt = nodestack_pop(ns);
  }
  free(ns);
}
