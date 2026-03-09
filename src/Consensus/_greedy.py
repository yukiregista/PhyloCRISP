from ._consensus import TreeList_with_support, Tree_with_support, _compatible, SQD_loss
from bitstring import Bits
import numpy as np
from itertools import combinations
import distance
from collections import OrderedDict
import dendropy

def _create_refinfo_splitkey(bipartitions, n_taxa):
    """Helper function for TBE

    Parameters
    ----------
    bipartitions : List of integer representation of bipartitions
        _description_
    n_taxa : int
        _description_

    Returns
    -------
    dict
        bipar_int as key an tuple of (bitstr, p) as value.
    """
    
    # bipartitions: iterable of dendropy.Bipartition
    refinfo = dict()
    for bipartition in bipartitions:
        edge_bitstr = Bits(uint = bipartition, length=n_taxa)
        counts = (edge_bitstr.count(0), edge_bitstr.count(1))
        if counts[0] < counts[1]:
            # zero lready assigned to light side
            bitstr = edge_bitstr
            p = counts[0]
        else:
            # reverse bitstring
            bitstr = (~edge_bitstr)
            p = counts[1]
        refinfo[bipartition] = (bitstr, p)
    # O( len(bipartitions) *  n_taxa )
    return refinfo



def _MinHammingDist(a : Bits, b : Bits):
    return min(distance.hamming(a.bin,b.bin), distance.hamming((~a).bin,b.bin))


class GreedyConsensusBase():
    
    def _Compatible_with_bipar_ints(self, bipar_int, tree_bipar_ints):
        
        if len(tree_bipar_ints) == 0:
            return True
        elif len(tree_bipar_ints) == self.n_taxa - 3:
            return False
        else:
            # check compatibility with each edge.
            tree_bipar_bits = [Bits(uint=key, length=self.n_taxa) for key in tree_bipar_ints]
            compatible = True
            cand_bit = Bits(uint=bipar_int,length=self.n_taxa)
            for bit in tree_bipar_bits:
                if not _compatible(cand_bit, bit):
                    return False
        return compatible
    
    def _CreateBipartitionListAndCounts(self):
        # create BipartitionDict
        BipartitionCountDict = OrderedDict()
        BipartitionDict = OrderedDict()
        for tree in self.input_trees:
            dict_keys = BipartitionCountDict.keys() # no duplicate keys in one tree
            tree.encode_bipartitions()
            for branch in tree.internal_edges(exclude_seed_edge=True):
                key = branch.bipartition.split_as_int()
                BipartitionDict.setdefault(key, branch.bipartition)
                if key in dict_keys:
                    BipartitionCountDict[key] += 1
                else:
                    BipartitionCountDict[key] = 1
        # Create BipartitionList, BipartitionCounts, Index2Bipartition and BipartitionP
        BipartitionList = []; BipartitionCounts = []; BipartitionP = []; BipartitionBitsList = []
        
        for key, value in BipartitionCountDict.items():
            BipartitionList.append(key)
            BipartitionCounts.append(value)
            bitstr = Bits(uint = key, length=self.n_taxa)
            BipartitionBitsList.append(bitstr)
            BipartitionP.append(min(bitstr.count(0), bitstr.count(1))) # size of the light side
        
        Bipartition2Index = {BipartitionList[i]:i for i in range(len(BipartitionList))}
        return np.array(BipartitionList), np.array(BipartitionCounts), BipartitionBitsList, np.array(BipartitionP), Bipartition2Index, BipartitionDict
    
    def __init__(self, input_trees : TreeList_with_support):
        self.input_trees = input_trees
        self.taxon_namespace = input_trees.taxon_namespace
        self.n_taxa = len(input_trees.taxon_namespace)

        self.BipartitionList, self.BipartitionCounts, self.BipartitionBitsList, self.BipartitionP, self.Bipartition2Index, self.BipartitionDict \
            = self._CreateBipartitionListAndCounts()
        
        self.n_bipartitions = len(self.BipartitionList)
        
        # Initialize current tree
        self.current_tree_included = np.zeros(self.n_bipartitions)
    
    
    def is_compatible(self, BiparBit):
        
        included_indices = np.nonzero(self.current_tree_included)[0]
        
        if len(included_indices) == 0:
            return True
        elif len(included_indices) == self.n_taxa - 3:
            return False
        else:
            # check compatibility with each edge.
            compatible = True
            for index in included_indices:
                current = self.BipartitionBitsList[index]
                if not _compatible(BiparBit, current):
                    return False
        return compatible
    
    def specify_initial_tree(self, initial_tree: Tree_with_support, *args, **kwargs):
        raise NotImplementedError()

    def reset_initial_tree(self, *args, **kwargs):
        raise NotImplementedError()

    def return_current_tree(self):
        bipartition_keys = self.BipartitionList[np.nonzero(self.current_tree_included)[0]]
        bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.n_taxa - 1) for item in bipartition_keys]
        dendropy_tree = dendropy.Tree.from_bipartition_encoding(bipartitions, self.taxon_namespace)
        return Tree_with_support(dendropy_tree)  
    
    def current_loss(self):
        raise NotImplementedError()
    
    def add(self, branch_index, *args, **kwargs):
        return NotImplementedError()
    
    def remove(self, branch_index, *args, **kwargs):
        return NotImplementedError()




class STDGreedyPrune(GreedyConsensusBase):
    pass

class SUTDGreedyPrune(GreedyConsensusBase):
    pass

class SQDGreedyPRune(GreedyConsensusBase):
    pass

class STDGreedyConsensus(GreedyConsensusBase):           
    def _ComputeDIST(self, method="simple"):
        DIST = np.zeros((self.n_bipartitions+1, self.n_bipartitions+1))
        if method == "simple":
            for i, j in combinations([i for i in range(self.n_bipartitions)], 2):
                bit_i = self.BipartitionBitsList[i]
                bit_j = self.BipartitionBitsList[j]
                DIST[i,j] = DIST[j,i] = _MinHammingDist(bit_i, bit_j)
        elif method == "efficient":
            refinfos = _create_refinfo_splitkey(self.BipartitionList, self.n_taxa)  
            for ind, bipar_key in enumerate(self.BipartitionList):
                refinfo_b = refinfos[bipar_key]
                p = refinfo_b[1]
                b_bitstr = refinfo_b[0]
                d = p - 1
                for tree in self.input_trees:
                    # WE ASSUME THAT bipartition and tree has the exact same TAXON_NAMESPACE 
                    taxon_namespace = tree.taxon_namespace
                    taxon_labels = [taxon.label for taxon in taxon_namespace]
                    n_taxa = len(taxon_labels)
                    taxon_labels_to_bitstr_digit = dict(zip(taxon_labels, [i for i in range(n_taxa)]))
                    node_biparint_to_postorder_index = dict()

                    if tree.bipartition_encoding is None:
                        tree.encode_bipartitions()
                    
                    m = len(tree.bipartition_encoding)
                    CountOnesSubtree = np.zeros(m)
                    CountZerosSubtree = np.zeros(m)

                    for i, node in enumerate(tree.postorder_node_iter()):
                        node_biparint_to_postorder_index[node.bipartition.split_as_int()] = i
                        if node.is_leaf():
                            digit = taxon_labels_to_bitstr_digit[node.taxon.label]
                            CountOnesSubtree[i] =  int(b_bitstr[- digit - 1])
                            CountZerosSubtree[i] = 1 - CountOnesSubtree[i]
                        else:
                            for child in node.child_node_iter():
                                CountOnesSubtree[i] += CountOnesSubtree[ node_biparint_to_postorder_index[child.bipartition.split_as_int()] ]
                                CountZerosSubtree[i] += CountZerosSubtree[ node_biparint_to_postorder_index[child.bipartition.split_as_int()] ]
                            actDist = p - CountZerosSubtree[i] + CountOnesSubtree[i]
                            if actDist > n_taxa / 2:
                                actDist = n_taxa - actDist
                            
                            if not node == tree.seed_node:
                                node_ind = self.Bipartition2Index[node.bipartition.split_as_int()]
                                DIST[ind, node_ind] = DIST[node_ind, ind] = actDist                            
                            if actDist < d:
                                d = actDist  
        else:
            raise ValueError(f"Invalid method: {method}")
        DIST[self.n_bipartitions, :-1] = self.BipartitionP - 1
        DIST[:-1, self.n_bipartitions] = self.BipartitionP - 1
        return DIST

    
    def _ComputeTransferDissimilarityCost(self):
        TD = np.zeros(self.n_bipartitions)
        for tree in self.input_trees:
            internal_bipartitions = [item.bipartition.split_as_int() for item in tree.internal_edges(exclude_seed_edge=True)]
            indices = [self.Bipartition2Index[item] for item in internal_bipartitions]
            for bipar in [i for i in range(self.n_bipartitions)]:
                minDist = np.min(self.DIST[bipar, indices])
                TD[bipar] += min (minDist / (self.BipartitionP[bipar]-1), 1)
        return TD
            
    
    def __init__(self, input_trees : TreeList_with_support):
        super().__init__(input_trees)
        
        self.DIST = self._ComputeDIST(method="simple")
        
        self.TD = self._ComputeTransferDissimilarityCost()
        
        
        # Initialize self.MATCH and self.SECOND_MATCH
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        # self.n_bipartitions is used to represent the match with leaf bipartitions.
        
    def specify_initial_tree(self, initial_tree : Tree_with_support, *args, **kwargs):
        # initial tree should only have those bipartitions in the input trees
        assert (initial_tree.taxon_namespace == self.taxon_namespace)
        internal_branches = initial_tree.internal_edges(exclude_seed_edge = True)
        internal_bipar_keys = [item.bipartition.split_as_int() for item in internal_branches]
        if not np.all([item in self.BipartitionList for item in internal_bipar_keys]):
            raise ValueError("Initial tree should only have bipartitions present in the input trees.")
        
        # renew self.current_tree_included
        self.current_tree_included = np.zeros(self.n_bipartitions)
        internal_indices = np.array([self.Bipartition2Index[item] for item in internal_bipar_keys])
        for item in internal_indices:
            self.current_tree_included[item] = 1
            
        # renew MATCH and SECOND_MATCH
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        
        if len(internal_indices) == 0:
            return
        
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            min_ind = internal_indices[np.argmin(self.DIST[bipar_ind, internal_indices])]
            if (self.DIST[bipar_ind, min_ind]) / (self.BipartitionP[bipar_ind]-1) >= 1:
                continue
            self.MATCH[bipar_ind] = min_ind
            candidate_indices = internal_indices[internal_indices != min_ind]
            if len(candidate_indices) == 0:
                break
            min2_ind = candidate_indices[np.argmin(self.DIST[bipar_ind, candidate_indices])]
            if (self.DIST[bipar_ind, min2_ind]) / (self.BipartitionP[bipar_ind]-1) >= 1:
                continue
            self.SECOND_MATCH[bipar_ind] = min2_ind
    
    def reset_initial_tree(self, *args, **kwargs):
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.current_tree_included = np.zeros(self.n_bipartitions)
        
    
    def is_compatible(self, BiparBit):
        
        included_indices = np.nonzero(self.current_tree_included)[0]
        
        if len(included_indices) == 0:
            return True
        elif len(included_indices) == self.n_taxa - 3:
            return False
        else:
            # check compatibility with each edge.
            compatible = True
            for index in included_indices:
                current = self.BipartitionBitsList[index]
                if not _compatible(BiparBit, current):
                    return False
        return compatible
                
    def _remove_benefit(self, branch_index):
        benefit_fp = self.TD[branch_index]
        # compute benefit_fn
        benefit_fn = 0
        branch = self.BipartitionList[branch_index]
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            # check if match
            if self.MATCH[bipar_ind] == branch_index:
                benefit_fn += ( (self.DIST[bipar_ind, branch_index] - self.DIST[bipar_ind,self.SECOND_MATCH[bipar_ind]]) / (self.BipartitionP[bipar_ind]-1) ) * self.BipartitionCounts[bipar_ind]
        
        return benefit_fp + benefit_fn 

    
    def _add_benefit(self, branch_index):
        benefit_fp = - self.TD[branch_index]
        # compute benefit_fn
        benefit_fn = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            dist_diff = self.DIST[self.MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
            if dist_diff > 0:
                # risk reduction occurs for this bipartition
                benefit_fn += dist_diff / (self.BipartitionP[bipar_ind] - 1) * self.BipartitionCounts[bipar_ind]
        
        return benefit_fn + benefit_fp

    
    def benefit(self, branch_index, method = "add"):
        if method == "add":
            return self._add_benefit(branch_index)
        elif method == "remove":
            return self._remove_benefit(branch_index)
        else:
            raise ValueError(f"Invalid method specified: {method}")
    
    
    def add(self, branch_index, *args, **kwargs):
        self.current_tree_included[branch_index] = 1
        # renew match and second match
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            first_dist_diff = self.DIST[self.MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
            if first_dist_diff > 0:
                # renew match and second match
                self.SECOND_MATCH[bipar_ind] = self.MATCH[bipar_ind]
                self.MATCH[bipar_ind] = branch_index
            else:
                second_dist_diff = self.DIST[self.SECOND_MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
                if second_dist_diff > 0:
                    # renew second match
                    self.SECOND_MATCH[bipar_ind] = branch_index
        
    def remove(self, branch_index, *args, **kwargs):
        self.current_tree_included[branch_index] = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            # check if match
            if self.MATCH[bipar_ind] == branch_index:
                # renew match and search for the new second match
                self.MATCH[bipar_ind] = self.SECOND_MATCH[bipar_ind]
                # Search for second match
                if self.MATCH[bipar_ind] == self.n_bipartitions:
                    self.SECOND_MATCH[bipar_ind] == self.n_bipartitions # equivalent to doing nothing
                else:
                    current_indices = np.nonzero(self.current_tree_included)[0]
                    candidate_indices = current_indices[current_indices!=self.MATCH[bipar_ind]]
                    cand = candidate_indices[np.argmin(self.DIST[bipar_ind,candidate_indices])]
                    if self.DIST[bipar_ind, cand]/(self.BipartitionP[bipar_ind] - 1) >= 1:
                        self.SECOND_MATCH[bipar_ind] = self.n_bipartitions
                    else:
                        self.SECOND_MATCH[bipar_ind] = cand
            if self.SECOND_MATCH[bipar_ind] == branch_index:
                # just renew the second match
                current_indices = np.nonzero(self.current_tree_included)[0]
                candidate_indices = current_indices[current_indices!=self.MATCH[bipar_ind]]
                cand = candidate_indices[np.argmin(self.DIST[bipar_ind,candidate_indices])]
                if self.DIST[bipar_ind, cand]/(self.BipartitionP[bipar_ind] - 1) >= 1:
                    self.SECOND_MATCH[bipar_ind] = self.n_bipartitions
                else:
                    self.SECOND_MATCH[bipar_ind] = cand

    def _first_greedy(self, sorted_index):
        improving = True
        while improving:
            improving = False
            for bipar_ind in sorted_index:
                if self.current_tree_included[bipar_ind] == 1:
                    if self.benefit(bipar_ind, method="remove") > 0:
                        self.remove(bipar_ind)
                        improving =True
                        break
                elif self.is_compatible(self.BipartitionBitsList[bipar_ind]):
                    if self.benefit(bipar_ind, method="add") > 0:
                        self.add(bipar_ind)
                        improving=True
                        break
        
    def _risk_greedy(self, sorted_index):
        while True:
            benefit = np.zeros(self.n_bipartitions)
            for bipar_ind in sorted_index:
                if self.current_tree_included[bipar_ind] == 1:
                    benefit[bipar_ind] = self.benefit(bipar_ind, method="remove")
                elif self.is_compatible(self.BipartitionBitsList[bipar_ind]):
                    benefit[bipar_ind] = self.benefit(bipar_ind, method="add")
            most_benefit_arg = np.argmax(benefit)
            if benefit[most_benefit_arg] > 0:
                if self.current_tree_included[most_benefit_arg] == 1:
                    self.remove(most_benefit_arg)
                else:
                    self.add(most_benefit_arg)
            else:
                break
        
    def greedy(self, method = "first", order="BS", refine_majority=False):
        
        if order == "BS":
            # order branch by branch support
            if refine_majority:
                # specify initial tree
                maj = self.input_trees.majority_consensus()
                self.specify_initial_tree(maj)
                
                # edit sorted index
                maj_ints = [branch.bipartition.split_as_int() for branch in maj.internal_edges(exclude_seed_edge=True)]
                mask = []
                for index, bipar_int in enumerate(self.BipartitionList):
                    if self._Compatible_with_bipar_ints(bipar_int, maj_ints):
                        mask.append(index)
                # now mask contains indices of those branches that are compatible with majority
                
                sorted_index = np.argsort(self.BipartitionCounts[mask])[::-1]
                sorted_index = [mask[item] for item in sorted_index]
            else:
                sorted_index = np.argsort(self.BipartitionCounts)[::-1]
               
            if method == "first":
                self._first_greedy(sorted_index)
            if method == "most":
                self._risk_greedy(sorted_index)
            
        
    def return_current_tree(self):
        bipartition_keys = self.BipartitionList[np.nonzero(self.current_tree_included)[0]]
        bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.n_taxa - 1) for item in bipartition_keys]
        dendropy_tree = dendropy.Tree.from_bipartition_encoding(bipartitions, self.taxon_namespace)
        return Tree_with_support(dendropy_tree)  
    
    def current_loss(self):
        nonzero = np.nonzero(self.current_tree_included)[0]
        fp = np.sum(self.TD[nonzero])
        fn = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            fn += ((self.DIST[bipar_ind, self.MATCH[bipar_ind]]) / (self.BipartitionP[bipar_ind] - 1)) * self.BipartitionCounts[bipar_ind]
        return (fp + fn) / len(self.input_trees)
                                
                   
class SQDGreedyConsensus(GreedyConsensusBase):
              
    def __init__(self, input_trees : TreeList_with_support):
        super().__init__(input_trees)
        self.input_trees_string = input_trees.as_string("newick", suppress_rooting=True)
        
        
     
    def specify_initial_tree(self, initial_tree : Tree_with_support, *args, **kwargs):
        # initial tree should only have those bipartitions in the input trees
        assert (initial_tree.taxon_namespace == self.taxon_namespace)
        internal_branches = initial_tree.internal_edges(exclude_seed_edge = True)
        internal_bipar_keys = [item.bipartition.split_as_int() for item in internal_branches]
        if not np.all([item in self.BipartitionList for item in internal_bipar_keys]):
            raise ValueError("Initial tree should only have bipartitions present in the input trees.")
        
        # renew self.current_tree_included
        self.current_tree_included = np.zeros(self.n_bipartitions)
        internal_indices = np.array([self.Bipartition2Index[item] for item in internal_bipar_keys])
        for item in internal_indices:
            self.current_tree_included[item] = 1
    
    def reset_initial_tree(self, *args, **kwargs):
        self.current_tree_included = np.zeros(self.n_bipartitions)
    
    def SQD_cost(self, bipar_keys, exec_dir=None):    
        bipartition_list = [self.BipartitionDict[key] for key in bipar_keys]
        tree = dendropy.Tree.from_bipartition_encoding(bipartition_list, taxon_namespace = self.taxon_namespace)
        return SQD_loss(tree, self.input_trees_string, len(self.input_trees), False, exec_dir)
     
    def add(self, branch_index, *args, **kwargs):
        self.current_tree_included[branch_index] = 1
        
    def remove(self, branch_index, *args, **kwargs):
        self.current_tree_included[branch_index] = 0
    
    def _first_greedy(self, sorted_index, exec_dir = None):
        
        mask = self.current_tree_included == 1
        bipar_keys = self.BipartitionList[mask]
        current_cost = self.SQD_cost(bipar_keys)
        
        improving = True
        while improving:
            improving = False
            for bipar_ind in sorted_index:
                if self.current_tree_included[bipar_ind] == 1:
                    new_bipar_keys = bipar_keys[bipar_keys!=bipar_ind]
                    newcost = self.SQD_cost(new_bipar_keys, exec_dir=exec_dir)
                    if newcost < current_cost:
                        self.remove(bipar_ind)
                        improving =True
                        current_cost = newcost
                        bipar_keys = new_bipar_keys
                        break
                elif self.is_compatible(self.BipartitionBitsList[bipar_ind]):
                    new_bipar_keys = np.append(bipar_keys, self.BipartitionList[bipar_ind])
                    newcost = self.SQD_cost(new_bipar_keys, exec_dir=exec_dir)
                    if newcost < current_cost:
                        self.add(bipar_ind)
                        improving=True
                        current_cost = newcost
                        bipar_keys = new_bipar_keys
                        break
    
    def greedy(self, method="first", order="BS", exec_dir = None, refine_majority=False):
        if order == "BS":
            # order branch by branch support
            if refine_majority:
                # specify initial tree
                maj = self.input_trees.majority_consensus()
                self.specify_initial_tree(maj)
                
                # edit sorted index
                maj_ints = [branch.bipartition.split_as_int() for branch in maj.internal_edges(exclude_seed_edge=True)]
                mask = []
                for index, bipar_int in enumerate(self.BipartitionList):
                    if self._Compatible_with_bipar_ints(bipar_int, maj_ints):
                        mask.append(index)
                # now mask contains indices of those branches that are compatible with majority
                
                sorted_index = np.argsort(self.BipartitionCounts[mask])[::-1]
                sorted_index = [mask[item] for item in sorted_index]
            else:
                sorted_index = np.argsort(self.BipartitionCounts)[::-1]
            if method == "first":
                self._first_greedy(sorted_index, exec_dir = exec_dir)
            if method == "most":
                NotImplementedError("Method 'most' is not implemented for SQDGreedyConsensus.")
               

class SUTDGreedyConsensus(GreedyConsensusBase):
     
    def _ComputeDIST(self, method="simple"):
        DIST = np.zeros((self.n_bipartitions+1, self.n_bipartitions+1))
        if method == "simple":
            for i, j in combinations([i for i in range(self.n_bipartitions)], 2):
                bit_i = self.BipartitionBitsList[i]
                bit_j = self.BipartitionBitsList[j]
                DIST[i,j] = DIST[j,i] = _MinHammingDist(bit_i, bit_j)
        elif method == "efficient":
            refinfos = _create_refinfo_splitkey(self.BipartitionList, self.n_taxa)  
            for ind, bipar_key in enumerate(self.BipartitionList):
                refinfo_b = refinfos[bipar_key]
                p = refinfo_b[1]
                b_bitstr = refinfo_b[0]
                d = p - 1
                for tree in self.input_trees:
                    # WE ASSUME THAT bipartition and tree has the exact same TAXON_NAMESPACE 
                    taxon_namespace = tree.taxon_namespace
                    taxon_labels = [taxon.label for taxon in taxon_namespace]
                    n_taxa = len(taxon_labels)
                    taxon_labels_to_bitstr_digit = dict(zip(taxon_labels, [i for i in range(n_taxa)]))
                    node_biparint_to_postorder_index = dict()

                    if tree.bipartition_encoding is None:
                        tree.encode_bipartitions()
                    
                    m = len(tree.bipartition_encoding)
                    CountOnesSubtree = np.zeros(m)
                    CountZerosSubtree = np.zeros(m)

                    for i, node in enumerate(tree.postorder_node_iter()):
                        node_biparint_to_postorder_index[node.bipartition.split_as_int()] = i
                        if node.is_leaf():
                            digit = taxon_labels_to_bitstr_digit[node.taxon.label]
                            CountOnesSubtree[i] =  int(b_bitstr[- digit - 1])
                            CountZerosSubtree[i] = 1 - CountOnesSubtree[i]
                        else:
                            for child in node.child_node_iter():
                                CountOnesSubtree[i] += CountOnesSubtree[ node_biparint_to_postorder_index[child.bipartition.split_as_int()] ]
                                CountZerosSubtree[i] += CountZerosSubtree[ node_biparint_to_postorder_index[child.bipartition.split_as_int()] ]
                            actDist = p - CountZerosSubtree[i] + CountOnesSubtree[i]
                            if actDist > n_taxa / 2:
                                actDist = n_taxa - actDist
                            
                            if not node == tree.seed_node:
                                node_ind = self.Bipartition2Index[node.bipartition.split_as_int()]
                                DIST[ind, node_ind] = DIST[node_ind, ind] = actDist                            
                            if actDist < d:
                                d = actDist  
        else:
            raise ValueError(f"Invalid method: {method}")
        DIST[self.n_bipartitions, :-1] = self.BipartitionP - 1
        DIST[:-1, self.n_bipartitions] = self.BipartitionP - 1
        return DIST

    
    def _ComputeUnnormalizedTransferDissimilarityCost(self):
        TD = np.zeros(self.n_bipartitions)
        for tree in self.input_trees:
            internal_bipartitions = [item.bipartition.split_as_int() for item in tree.internal_edges(exclude_seed_edge=True)]
            indices = [self.Bipartition2Index[item] for item in internal_bipartitions]
            for bipar in [i for i in range(self.n_bipartitions)]:
                minDist = np.min(self.DIST[bipar, indices])
                TD[bipar] += min (minDist, self.BipartitionP[bipar]-1)
        return TD
    
    def __init__(self, input_trees : TreeList_with_support):
        super().__init__(input_trees)
        
        self.DIST = self._ComputeDIST(method="simple")
        
        self.UTD = self._ComputeUnnormalizedTransferDissimilarityCost()
        
        
        # Initialize self.MATCH and self.SECOND_MATCH
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        # self.n_bipartitions is used to represent the match with leaf bipartitions.
        
    def specify_initial_tree(self, initial_tree : Tree_with_support, *args, **kwargs):
        # initial tree should only have those bipartitions in the input trees
        assert (initial_tree.taxon_namespace == self.taxon_namespace)
        internal_branches = initial_tree.internal_edges(exclude_seed_edge = True)
        internal_bipar_keys = [item.bipartition.split_as_int() for item in internal_branches]
        if not np.all([item in self.BipartitionList for item in internal_bipar_keys]):
            raise ValueError("Initial tree should only have bipartitions present in the input trees.")
        
        # renew self.current_tree_included
        self.current_tree_included = np.zeros(self.n_bipartitions)
        internal_indices = np.array([self.Bipartition2Index[item] for item in internal_bipar_keys])
        for item in internal_indices:
            self.current_tree_included[item] = 1
            
        # renew MATCH and SECOND_MATCH
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        
        if len(internal_indices) == 0:
            return
        
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            min_ind = internal_indices[np.argmin(self.DIST[bipar_ind, internal_indices])]
            if (self.DIST[bipar_ind, min_ind]) / (self.BipartitionP[bipar_ind]-1) >= 1:
                continue
            self.MATCH[bipar_ind] = min_ind
            candidate_indices = internal_indices[internal_indices != min_ind]
            if len(candidate_indices) == 0:
                break
            min2_ind = candidate_indices[np.argmin(self.DIST[bipar_ind, candidate_indices])]
            if (self.DIST[bipar_ind, min2_ind]) / (self.BipartitionP[bipar_ind]-1) >= 1:
                continue
            self.SECOND_MATCH[bipar_ind] = min2_ind
    
    def reset_initial_tree(self, *args, **kwargs):
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.current_tree_included = np.zeros(self.n_bipartitions)
        
    
    def is_compatible(self, BiparBit):
        
        included_indices = np.nonzero(self.current_tree_included)[0]
        
        if len(included_indices) == 0:
            return True
        elif len(included_indices) == self.n_taxa - 3:
            return False
        else:
            # check compatibility with each edge.
            compatible = True
            for index in included_indices:
                current = self.BipartitionBitsList[index]
                if not _compatible(BiparBit, current):
                    return False
        return compatible
                
    def _remove_benefit(self, branch_index):
        benefit_fp = self.UTD[branch_index]
        # compute benefit_fn
        benefit_fn = 0
        branch = self.BipartitionList[branch_index]
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            # check if match
            if self.MATCH[bipar_ind] == branch_index:
                benefit_fn += (self.DIST[bipar_ind, branch_index] - self.DIST[bipar_ind,self.SECOND_MATCH[bipar_ind]])  * self.BipartitionCounts[bipar_ind]
        
        return benefit_fp + benefit_fn 

    
    def _add_benefit(self, branch_index):
        benefit_fp = - self.UTD[branch_index]
        # compute benefit_fn
        benefit_fn = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            dist_diff = self.DIST[self.MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
            if dist_diff > 0:
                # risk reduction occurs for this bipartition
                benefit_fn += dist_diff * self.BipartitionCounts[bipar_ind]
        
        return benefit_fn + benefit_fp

    
    def benefit(self, branch_index, method = "add"):
        if method == "add":
            return self._add_benefit(branch_index)
        elif method == "remove":
            return self._remove_benefit(branch_index)
        else:
            raise ValueError(f"Invalid method specified: {method}")
    
    
    def add(self, branch_index, *args, **kwargs):
        self.current_tree_included[branch_index] = 1
        # renew match and second match
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            first_dist_diff = self.DIST[self.MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
            if first_dist_diff > 0:
                # renew match and second match
                self.SECOND_MATCH[bipar_ind] = self.MATCH[bipar_ind]
                self.MATCH[bipar_ind] = branch_index
            else:
                second_dist_diff = self.DIST[self.SECOND_MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
                if second_dist_diff > 0:
                    # renew second match
                    self.SECOND_MATCH[bipar_ind] = branch_index
        
    def remove(self, branch_index, *args, **kwargs):
        self.current_tree_included[branch_index] = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            # check if match
            if self.MATCH[bipar_ind] == branch_index:
                # renew match and search for the new second match
                self.MATCH[bipar_ind] = self.SECOND_MATCH[bipar_ind]
                # Search for second match
                if self.MATCH[bipar_ind] == self.n_bipartitions:
                    self.SECOND_MATCH[bipar_ind] == self.n_bipartitions # equivalent to doing nothing
                else:
                    current_indices = np.nonzero(self.current_tree_included)[0]
                    candidate_indices = current_indices[current_indices!=self.MATCH[bipar_ind]]
                    cand = candidate_indices[np.argmin(self.DIST[bipar_ind,candidate_indices])]
                    if self.DIST[bipar_ind, cand]/(self.BipartitionP[bipar_ind] - 1) >= 1:
                        self.SECOND_MATCH[bipar_ind] = self.n_bipartitions
                    else:
                        self.SECOND_MATCH[bipar_ind] = cand
            if self.SECOND_MATCH[bipar_ind] == branch_index:
                # just renew the second match
                current_indices = np.nonzero(self.current_tree_included)[0]
                candidate_indices = current_indices[current_indices!=self.MATCH[bipar_ind]]
                cand = candidate_indices[np.argmin(self.DIST[bipar_ind,candidate_indices])]
                if self.DIST[bipar_ind, cand]/(self.BipartitionP[bipar_ind] - 1) >= 1:
                    self.SECOND_MATCH[bipar_ind] = self.n_bipartitions
                else:
                    self.SECOND_MATCH[bipar_ind] = cand

    def _first_greedy(self, sorted_index):
        improving = True
        while improving:
            improving = False
            for bipar_ind in sorted_index:
                if self.current_tree_included[bipar_ind] == 1:
                    if self.benefit(bipar_ind, method="remove") > 0:
                        self.remove(bipar_ind)
                        improving =True
                        break
                elif self.is_compatible(self.BipartitionBitsList[bipar_ind]):
                    if self.benefit(bipar_ind, method="add") > 0:
                        self.add(bipar_ind)
                        improving=True
                        break
        
    def _risk_greedy(self, sorted_index):
        while True:
            benefit = np.zeros(self.n_bipartitions)
            for bipar_ind in sorted_index:
                if self.current_tree_included[bipar_ind] == 1:
                    benefit[bipar_ind] = self.benefit(bipar_ind, method="remove")
                elif self.is_compatible(self.BipartitionBitsList[bipar_ind]):
                    benefit[bipar_ind] = self.benefit(bipar_ind, method="add")
            most_benefit_arg = np.argmax(benefit)
            if benefit[most_benefit_arg] > 0:
                if self.current_tree_included[most_benefit_arg] == 1:
                    self.remove(most_benefit_arg)
                else:
                    self.add(most_benefit_arg)
            else:
                break
        
    def greedy(self, method = "first", order="BS", refine_majority=False):
        
        if order == "BS":
            if refine_majority:
                # specify initial tree
                maj = self.input_trees.majority_consensus()
                self.specify_initial_tree(maj)
                
                # edit sorted index
                maj_ints = [branch.bipartition.split_as_int() for branch in maj.internal_edges(exclude_seed_edge=True)]
                mask = []
                for index, bipar_int in enumerate(self.BipartitionList):
                    if self._Compatible_with_bipar_ints(bipar_int, maj_ints):
                        mask.append(index)
                # now mask contains indices of those branches that are compatible with majority
                
                sorted_index = np.argsort(self.BipartitionCounts[mask])[::-1]
                sorted_index = [mask[item] for item in sorted_index]
            else:
                sorted_index = np.argsort(self.BipartitionCounts)[::-1]
            # order branch by branch support
            if method == "first":
                self._first_greedy(sorted_index)
            if method == "most":
                self._risk_greedy(sorted_index)
            
        
    def return_current_tree(self):
        bipartition_keys = self.BipartitionList[np.nonzero(self.current_tree_included)[0]]
        bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.n_taxa - 1) for item in bipartition_keys]
        dendropy_tree = dendropy.Tree.from_bipartition_encoding(bipartitions, self.taxon_namespace)
        return Tree_with_support(dendropy_tree)  
    
    def current_loss(self):
        nonzero = np.nonzero(self.current_tree_included)[0]
        fp = np.sum(self.UTD[nonzero])
        fn = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            fn += (self.DIST[bipar_ind, self.MATCH[bipar_ind]])  * self.BipartitionCounts[bipar_ind]
        return (fp + fn) / len(self.input_trees)
 
