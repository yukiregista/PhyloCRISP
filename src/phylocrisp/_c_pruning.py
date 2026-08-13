from ._booster import load_booster, prepare_tbe_support_and_match, prepare_tbe_support_and_match_args, organize_support_and_match, OrganizedTransferResult, _generate_bitmask, prepare_recompute, transferDistanceAll, prepare_prune_and_return_newick, prepare_prune_and_return_newick_args
import dendropy
import numpy as np
import ctypes
from bitstring import Bits

import re


class PruningVars():
    """
        Class to hold all variables used in pruning.
    """
    def __init__(self, FP: np.ndarray, FN_diff: np.ndarray, E: set, R1: np.ndarray, R2: np.ndarray, W1: list[set], W2: list[set], OSR: OrganizedTransferResult, booster_lib: ctypes.CDLL):
        self.FP: np.ndarray = FP # False Positive Loss caused by each internal bipartitoin int the ref tree.
        self.FN_diff: np.ndarray = FN_diff # False Negaive Loss increase that will be induced by pruning each branch in the ref tree.
        self.E: set = E # Set of indices of ref_tree to keep
        self.pruned: set = set() # Set of indices that have been pruned.
        self.R1: np.ndarray = R1 # Rank of the first match (starts with index 0). 
        self.R2: np.ndarray = R2 # Rank of the second match (initialized with 1 at first).
        self.W1: list[set] = W1 # List of size maximum n-3, containing the set of branches matched to each internal bipar in ref tree.
        self.W2: list[set] = W2 # List of size maximum n-3, containing the set of branches second matched to each internal bipar in ref tree.
        self.L: np.ndarray = self.FP - self.FN_diff # Utility of pruning branches
        self.OSR: OrganizedTransferResult = OSR # OrganizedTransferResult
        self.use_K: list[bool] = [True for _ in range(len(self.R1))] # Whether to use the originally returned K matches, or the result of recomputation.
        self.booster_lib = booster_lib
        
        # some useful variables
        self.n_internal_ref = len(self.FP) # Number of internal edges in the initial tree
        self.m = len(self.R1) # Number of unique internal branches in the input trees
        self.K = self.OSR.M.shape[1] # K
        
        # For storing recomputation results
        self.RECOMP: dict[int, RecomputationResults] = dict()


class RecomputationResults():
    
    def __init__(self, M: np.ndarray, S: np.ndarray):
        self.M = M
        self.S = S
    
    
def transfer_distance(a: tuple, b: tuple, bits: int, n_taxa: int):
    assert(len(a)==len(b))
    spl = len(a)
    mask = _generate_bitmask(bits)
    mask_last = _generate_bitmask(n_taxa % bits)
    hamdist1 = 0
    for i in range(spl):
        hamdist1 += bin(a[i] ^ b[i]).count('1')
    
    hamdist2 = 0
    for i in range(spl-1):
        hamdist2 += bin(a[i] ^ (~b[i] & mask)).count('1')
    hamdist2 += bin(a[spl-1] ^  (~b[spl-1] & mask_last)).count('1')
    
    return min(hamdist1, hamdist2)
        

def recompute_td(osr: OrganizedTransferResult, bipar_id: int, booster_lib: ctypes.CDLL, scaled=True):
    h1 = osr.H1[bipar_id]
    h2 = osr.H2[bipar_id]
    tda: transferDistanceAll = booster_lib.recompute(ctypes.c_uint(h1), ctypes.c_uint(h2)).contents
    
    # We need to reconcile the order between the original OSR result and tda result.
    node_id_list = []
    tda_dict = dict()
    use_external=False
    for i in range(tda.n_elements):
        if (tda.td[i].contents.dist >= tda.topo_depth-1):
            continue
        else:  
            bipar_loc = osr.REF_BIPAR_ID_to_location[ tda.td[i].contents.node_id ]
            node_id_list.append( bipar_loc )
            if scaled:
                tda_dict[bipar_loc] = (tda.td[i].contents.dist / (tda.topo_depth-1))
            else:
                tda_dict[bipar_loc] = tda.td[i].contents.dist 


    if (len(node_id_list)) == 0:
        # Others are all filled with externals
        new_M = np.concatenate([osr.M[bipar_id], [-1 for _ in range(tda.n_elements - len(osr.M[bipar_id]))]]).astype(int)
        if scaled:
            new_S = np.concatenate([osr.S[bipar_id], [1 for _ in range(tda.n_elements - len(osr.M[bipar_id]))]])
        else:
            new_S = np.concatenate([osr.S[bipar_id], [tda.topo_depth - 1 for _ in range(tda.n_elements - len(osr.M[bipar_id]))]])
    else:
        first_K_IDS = set(osr.M[bipar_id])
        IDs_to_add = [item for item in node_id_list if item not in first_K_IDS]
        Values = [tda_dict[item] for item in IDs_to_add]
        # sort the values
        order = np.argsort(Values)
        M_mid = [IDs_to_add[x] for x in order]
        S_mid = [Values[x] for x in order]
        rem_elements = tda.n_elements - len(osr.M[bipar_id]) - len(M_mid)
        new_M = np.concatenate([ osr.M[bipar_id], M_mid, [-1 for _ in range(rem_elements)] ]).astype(int)
        if scaled:
            new_S = np.concatenate([ osr.S[bipar_id], S_mid, [1 for _ in range(rem_elements)] ])
        else:
            new_S = np.concatenate([ osr.S[bipar_id], S_mid, [tda.topo_depth-1 for _ in range(rem_elements)] ])
    
    return RecomputationResults(new_M, new_S)
            
    
     
def compute_all_normalized_td(bipar, matched, K, reverse_id_ref, bipar_val, topo_depth, bits, n_taxa):
    processed_set = set(matched)
    remaining_branches = [(id,val) for id, val in reverse_id_ref.items() if id not in processed_set]
    td = [transfer_distance(bipar_val, val[1], bits, n_taxa) for val in remaining_branches]
    order = np.argsort(td)
    ns = np.array([td[item]/(topo_depth - 1.0) for item in order])
    ids = [remaining_branches[item][0] for item in order]
    
    ext_index = 0
    for i in range(len(ns)):
        if ns[i] > topo_depth-2:
            ext_index=i
            break
    
    ids = [remaining_branches[item][0] for item in order[:ext_index]] + [-1 for _ in range(len(ns)-ext_index)]
    ns[ext_index:] = 1.0
    return ids, ns
    
    
def renew_match(pv: PruningVars, b: int, scaled=True):
    """
        renew match by pruning b.
    """
    pv.L[b] = -1 # Set to -1 so that it will not be picked again.
    for bipar_id in pv.W1[b]: # bipartitions in the first match
        # Promote second match to the fist match
        second_match = pv.OSR.M[bipar_id, pv.R2[bipar_id]] if pv.use_K[bipar_id] else pv.RECOMP[bipar_id].M[pv.R2[bipar_id]]
        if (second_match!=-1): # If the original second match of bipar_id is not external one.
            pv.W1[second_match].add(bipar_id) # Add bipar_id to the original second match's W1
            pv.W2[second_match].remove(bipar_id) # Remove bipar_id from the original second match's W2
            pv.R1[bipar_id] = pv.R2[bipar_id]
        
            # Find new second match
            now_rank = pv.R1[bipar_id]+1
            while True:
                if pv.use_K[bipar_id]:
                    if now_rank > pv.K-1: # Then do recomputation
                        # recomputation
                        
                        pv.RECOMP[bipar_id] = recompute_td(pv.OSR, bipar_id, pv.booster_lib, scaled=scaled)
                        pv.use_K[bipar_id] = False
                    elif pv.OSR.M[bipar_id, now_rank] in pv.pruned: # Then proceed to next rank
                        # it is already pruned
                        now_rank += 1
                    else: # Then store the new second match and modify losses
                        if pv.OSR.M[bipar_id, now_rank]!=-1: # if matched to one of the internal branches
                            pv.W2[pv.OSR.M[bipar_id, now_rank]].add(bipar_id)
                        pv.R2[bipar_id] = now_rank
                        # Adjust the FN_diff and L of the first match
                        diff = (pv.OSR.S[bipar_id, pv.R2[bipar_id]] - pv.OSR.S[bipar_id, pv.R1[bipar_id]]) * pv.OSR.C[bipar_id]
                        pv.FN_diff[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] += diff
                        pv.L[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] -= diff
                        break
                        
                else:
                    if pv.RECOMP[bipar_id].M[now_rank] in pv.pruned: # Then proceed to next rank
                        now_rank += 1
                    else: # Then store the new second match and modify losses
                        if pv.RECOMP[bipar_id].M[now_rank] != -1:
                            pv.W2[pv.RECOMP[bipar_id].M[now_rank]].add(bipar_id)
                        pv.R2[bipar_id] = now_rank
                        # Adjust the FN_diff and L of the first match
                        diff = (pv.RECOMP[bipar_id].S[pv.R2[bipar_id]] - pv.RECOMP[bipar_id].S[pv.R1[bipar_id]]) * pv.OSR.C[bipar_id]
                        pv.FN_diff[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] += diff
                        pv.L[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] -= diff
                        break
     
    for bipar_id in pv.W2[b]: # Bipartitions that had b as the second match -> renew second match
        now_rank = pv.R2[bipar_id] + 1
        while True:
            if pv.use_K[bipar_id]:
                if now_rank > pv.K-1: # Then do recomputation
                    # recomputation
                    pv.RECOMP[bipar_id] = recompute_td(pv.OSR, bipar_id, pv.booster_lib, scaled=scaled)
                    pv.use_K[bipar_id] = False
                elif pv.OSR.M[bipar_id, now_rank] in pv.pruned: # Then proceed to next rank
                    now_rank += 1
                else: # Then store the new second match and modify losses
                    if pv.OSR.M[bipar_id, now_rank]!=-1: # if matched to one of the internal branches
                        pv.W2[pv.OSR.M[bipar_id, now_rank]].add(bipar_id)
                    old_rank = pv.R2[bipar_id]
                    pv.R2[bipar_id] = now_rank
                    # Adjust the FN_diff and L of the first match
                    diff = (pv.OSR.S[bipar_id, pv.R2[bipar_id]] - pv.OSR.S[bipar_id, old_rank]) * pv.OSR.C[bipar_id]
                    pv.FN_diff[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] += diff
                    pv.L[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] -= diff
                    break
            else:
                if pv.RECOMP[bipar_id].M[now_rank] in pv.pruned: # Then proceed to next rank
                    now_rank += 1
                else: # Then store the new second match and modify losses
                    if pv.RECOMP[bipar_id].M[now_rank] != -1:
                        pv.W2[pv.RECOMP[bipar_id].M[now_rank]].add(bipar_id)
                    old_rank = pv.R2[bipar_id]
                    pv.R2[bipar_id] = now_rank
                    # Adjust the FN_diff and L of the first match
                    diff = (pv.RECOMP[bipar_id].S[pv.R2[bipar_id]] - pv.RECOMP[bipar_id].S[old_rank]) * pv.OSR.C[bipar_id]
                    pv.FN_diff[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] += diff
                    pv.L[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] -= diff
                    break      
                        
        
          
        
        
        
def greedy_pruning(pv: PruningVars, scaled=True):
    while True:
        best_index = np.argmax(pv.L)
        if (pv.L[best_index] <= 0): break # No edges induce loss reduction. Exiting the while loop
        # Otherwise, prune best_index
        pv.E.remove(best_index)
        pv.pruned.add(best_index)
        renew_match(pv, best_index, scaled=scaled)



def _prune_first_comp(res: OrganizedTransferResult, booster_lib: ctypes.CDLL, scaled=True):
    # False Positive Loss
    if scaled:
        FP_loss = (1-np.array(res.TS)) * res.num_alt_trees
    else:
        FP_loss = (1-np.array(res.TS)) * (np.array(res.TD)-1) * res.num_alt_trees
        
    # False Negative Losses and W1, W2
    W1 = [set() for _ in range(len(res.TS))]
    W2 = [set() for _ in range(len(res.TS))]
    FN_diff = np.zeros_like(FP_loss)
    for i in range(res.S.shape[0]):
        if res.M[i,0]!= -1:
            FN_diff[res.M[i,0]] += (res.S[i,1] - res.S[i,0]) * res.C[i]
            W1[res.M[i,0]].add(i)
            if res.M[i,1]!=-1:
                W2[res.M[i,1]].add(i)   
    # Initialize Other variables
    E = [i for i in range(len(res.TS))]
    R1 = np.zeros(res.S.shape[0], dtype=int)
    R2 = np.full(res.S.shape[0], 1, dtype=int)
    
    return PruningVars(FP_loss, FN_diff, E, R1, R2, W1, W2, res, booster_lib)
           
def split_set_to_bipar_encodings(bipar_set, n_taxa, bits):
    split_length = (n_taxa-1) // bits + 1
    last_bits = n_taxa % bits
    
    # first take unique list of all unsigned integers to convert.
    all_integers_pre = [item for integer_list in bipar_set for item in integer_list[:split_length]]
    all_integers_pre = np.unique(all_integers_pre)
    
    all_integers_last = [integer_list[-1] for integer_list in bipar_set]
    all_integers_last = np.unique(all_integers_last)
    
    # split to bipar
    pre_bits = [''.join('1' if item & (1 << i) else '0' for i in range(bits)) for item in all_integers_pre]
    pre_dict = dict(zip(all_integers_pre, pre_bits))
    last_bits = [''.join('1' if item & (1 << i) else '0' for i in range(last_bits)) for item in all_integers_last]
    last_dict = dict(zip(all_integers_last, last_bits))
    
    bipartition_list = []
    for bipar in bipar_set:
        bipartition_bin = ''.join([pre_dict[item] for item in bipar[:split_length-1]] + [last_dict[bipar[split_length-1]]])
        bipartition_list.append(dendropy.Bipartition(leafset_bitmask = Bits(bin=bipartition_bin).uint, tree_leafset_bitmask=2**n_taxa-1))
    return bipartition_list
 
def remove_branch_lengths(newick_string):
    # Use regex to remove branch lengths (e.g., ":0.123")
    return re.sub(r':-?\d+(\.\d+)?([eE][-+]?\d+)?', '', newick_string)

def c_prune(inittree_file: str, inputtrees_file: str, K=30, scaled=True):
    try:
        booster_lib = load_booster()
    except ImportError as e:
        raise ImportError(
            "The optional 'booster' native extension is not available. "
            "Fast prune-only mode for Scaled/Unscaled Transfer methods requires booster. "
            "On Windows, install a build that includes booster support, or use Add and Prune / Quartet methods."
        ) from e
    booster_lib = prepare_tbe_support_and_match(booster_lib)
    booster_lib = prepare_recompute(booster_lib)
    booster_lib = prepare_prune_and_return_newick(booster_lib)
    
    
    with open(inputtrees_file, "r") as f:
        inputtrees = f.read()
    with open(inittree_file, "r") as f:
        inittree = f.read()
    args = prepare_tbe_support_and_match_args(inittree, inputtrees, K)
    res_ptr = None
    try:
        res_ptr = booster_lib.tbe_support_and_match(*args)

        organized_results: OrganizedTransferResult = organize_support_and_match(res_ptr, scaled=scaled)
        
        pruning_vars: PruningVars = _prune_first_comp(organized_results, booster_lib, scaled=scaled)
        
        greedy_pruning(pruning_vars, scaled=scaled)
        
        loc_to_id = dict()
        for id, loc in pruning_vars.OSR.REF_BIPAR_ID_to_location.items():
            loc_to_id[loc] = id
        id_list = [loc_to_id[item] for item in pruning_vars.pruned]
        args = prepare_prune_and_return_newick_args(id_list)
        newick_string = booster_lib.prune_and_return_newick(*args)
        newick_string_decoded = ctypes.cast(newick_string, ctypes.c_char_p).value.decode('utf-8')
        
        booster_lib.free_buffer(newick_string)
    finally:
        if res_ptr is not None:
            booster_lib.free_tbe_support_and_match(res_ptr)
        if hasattr(booster_lib, "free_TDA_LIST"):
            booster_lib.free_TDA_LIST()

    
    
    return remove_branch_lengths(newick_string_decoded)
