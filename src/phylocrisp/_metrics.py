import os
import subprocess
import shutil
import tempfile

import numpy as np
from bitstring import Bits
import dendropy

from ._consensus import Tree_with_support, TreeList_with_support, transfer_support, unnormalized_transfer_support, _tqdist_fp_fn, _tqdist_fp_fn_trees_str, _randomname


def _compute_all_depths_dfs(tree):
    """
    Compute depths for all internal nodes in O(n) time using DFS.
    
    The depth of a bipartition is the minimum of the number of leaves on each side.
    
    Parameters
    ----------
    tree : dendropy.Tree
        Tree with encoded bipartitions
        
    Returns
    -------
    dict
        Mapping from node to depth value
    """
    depths = {}
    n_taxa = len(tree.taxon_namespace)
    
    def count_leaves(node):
        """DFS to count leaves in subtree and compute depth."""
        if node.is_leaf():
            return 1
        
        # Count leaves in left and right subtrees
        leaf_count = sum(count_leaves(child) for child in node.child_nodes())
        
        # For internal nodes (excluding root), compute depth
        if not node.is_seed_node():
            # depth = min(leaves on one side, leaves on other side)
            depths[node] = min(leaf_count, n_taxa - leaf_count)
        
        return leaf_count
    
    count_leaves(tree.seed_node)
    return depths


def _compute_metric_from_result(result, metric_type):
    """
    Compute metric from OrganizedTransferResult directly.
    
    Parameters
    ----------
    result : OrganizedTransferResult
        Result from organize_support_and_match in _booster.py
    metric_type : str
        'std', 'sutd', or 'sbd'
    
    Returns
    -------
    float
        Metric value
    """
    metric_val = 0.0
    
    if metric_type == 'std':
        # Sum of (1 - transfer_support)
        for ts in result. TS: 
            metric_val += 1.0 - ts
    
    elif metric_type == 'sutd':
        # Sum of (1 - transfer_support) * (depth - 1)
        for ts, td in zip(result.TS, result. TD):
            metric_val += (1.0 - ts) * (td - 1)
    
    elif metric_type == 'sbd':
        # Same as std for bipartition distance
        for ts in result. TS:
            metric_val += 1.0 - ts
    
    return metric_val


def _compute_fp_fn_with_booster_direct(estimate, inputs, metric_type, aggregate=False):
    """
    Compute FP and FN using direct booster library calls (no subprocess).
    
    Parameters
    ----------
    estimate : Tree_with_support
        Estimate tree
    inputs : Tree_with_support or TreeList_with_support or list
        Input tree(s)
    metric_type : str
        Metric type ('std' or 'sutd')
    aggregate : bool, optional
        If True, compute aggregated metrics (n+1 booster calls, returns scalars).
        If False (default), compute per-tree metrics (2n booster calls, returns vectors).
        
    Returns
    -------
    tuple
        (fp, fn) - scalars if aggregate=True or inputs is single tree, 
        arrays if aggregate=False and inputs is list/TreeList
    """
    from ._booster import (
        load_booster, 
        prepare_tbe_support_and_match,
        prepare_tbe_support_and_match_args,
        organize_support_and_match
    )
    
    # Load library once
    booster_lib = load_booster()
    booster_lib = prepare_tbe_support_and_match(booster_lib)
    
    # Convert trees to newick strings
    estimate_str = estimate.as_string(schema="newick", suppress_rooting=True).strip()
    
    inputs_list = list(inputs) if not isinstance(inputs, Tree_with_support) else [inputs]
    is_single_input = isinstance(inputs, Tree_with_support)
    
    # Determine if we use scaled (normalized) transfer distance
    scaled = (metric_type == 'std')
    
    try:
        if aggregate and not is_single_input:
            # Aggregate mode: 1 FP call (all inputs), n FN calls
            # FP: estimate supported by all inputs
            all_inputs_str = "\n".join(
                tree.as_string(schema="newick", suppress_rooting=True).strip() 
                for tree in inputs_list
            )
            
            ref_tree_cstr, input_trees_array, num_trees, K = \
                prepare_tbe_support_and_match_args(estimate_str, all_inputs_str, 1)
            
            res_ptr = booster_lib. tbe_support_and_match(
                ref_tree_cstr, input_trees_array, num_trees, K
            )
            
            try:
                result = organize_support_and_match(res_ptr, scaled=scaled)
                fp = _compute_metric_from_result(result, metric_type)
            finally:
                # Free memory
                booster_lib.free_tbe_support_and_match(res_ptr)
            
            # FN: each input supported by estimate (average)
            fn_list = []
            for input_tree in inputs_list:
                input_str = input_tree.as_string(schema="newick", suppress_rooting=True).strip()
                ref_tree_cstr, input_trees_array, num_trees, K = \
                    prepare_tbe_support_and_match_args(input_str, estimate_str, 1)
                
                res_ptr = booster_lib.tbe_support_and_match(
                    ref_tree_cstr, input_trees_array, num_trees, K
                )
                
                try:
                    result = organize_support_and_match(res_ptr, scaled=scaled)
                    fn_i = _compute_metric_from_result(result, metric_type)
                    fn_list.append(fn_i)
                finally:
                    # Free memory
                    booster_lib.free_tbe_support_and_match(res_ptr)
            
            return fp, np.mean(fn_list)
        
        else:
            # Per-tree mode: n FP calls, n FN calls (or single if is_single_input)
            fp_list = []
            fn_list = []
            
            for input_tree in inputs_list: 
                input_str = input_tree.as_string(schema="newick", suppress_rooting=True).strip()
                # FP: estimate supported by this input
                ref_tree_cstr, input_trees_array, num_trees, K = \
                    prepare_tbe_support_and_match_args(estimate_str, input_str, 1)
                
                res_ptr = booster_lib.tbe_support_and_match(
                    ref_tree_cstr, input_trees_array, num_trees, K
                )
                
                try:
                    result = organize_support_and_match(res_ptr, scaled=scaled)
                    fp_i = _compute_metric_from_result(result, metric_type)
                finally:
                    # Free memory
                    booster_lib.free_tbe_support_and_match(res_ptr)
                
                # FN: input supported by estimate
                ref_tree_cstr, input_trees_array, num_trees, K = \
                    prepare_tbe_support_and_match_args(input_str, estimate_str, 1)
                
                res_ptr = booster_lib.tbe_support_and_match(
                    ref_tree_cstr, input_trees_array, num_trees, K
                )
                
                try:
                    result = organize_support_and_match(res_ptr, scaled=scaled)
                    fn_i = _compute_metric_from_result(result, metric_type)
                finally:
                    # Free memory
                    booster_lib.free_tbe_support_and_match(res_ptr)
                
                fp_list.append(fp_i)
                fn_list.append(fn_i)
            
            if is_single_input:
                return fp_list[0], fn_list[0]
            else: 
                return np.array(fp_list), np.array(fn_list)
    
    finally:
        # Clean up any remaining TDA list memory
        if hasattr(booster_lib, 'free_TDA_LIST'):
            booster_lib. free_TDA_LIST()


def _SUTD1_fp_fn(true_tree, tree2:  Tree_with_support):
    tree1_int_bipars = [node. bipartition for node in true_tree.postorder_internal_node_iter(exclude_seed_node=True)]
    tree2_int_bipars = [node.bipartition for node in tree2.postorder_internal_node_iter(exclude_seed_node=True)]
    fp = np.sum(unnormalized_transfer_support(tree2_int_bipars, [true_tree]))
    fn = np.sum(unnormalized_transfer_support(tree1_int_bipars, [tree2]))
    return fp, fn


def _STD1_fp_fn(true_tree, tree2: Tree_with_support):
    tree1_int_bipars = [node.bipartition for node in true_tree.postorder_internal_node_iter(exclude_seed_node=True)]
    tree2_int_bipars = [node.bipartition for node in tree2.postorder_internal_node_iter(exclude_seed_node=True)]
    fp = np.sum((1-transfer_support(tree2_int_bipars, [true_tree])))
    fn = np.sum((1-transfer_support(tree1_int_bipars, [tree2])))
    return fp, fn


def SUTD_fp_fn(estimate, inputs, parent_dir=None, aggregate=False, use_booster=True):
    """Computes false positive and false negative of unnormalized STD. 
        
    Parameters
    ----------
    estimate : Tree_with_support
        Estimates. 
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate. 
    parent_dir : str, optional
        Path to place intermediate files, by default None (place files in the folder of execution).
    aggregate : bool, optional
        If True, compute aggregated metrics (n+1 booster calls, returns scalars).
        If False (default), compute per-tree metrics (2n booster calls, returns vectors).
    use_booster : bool, optional
        If True (default), attempt to use booster library for faster computation.
        If False, use pure Python implementation. 

    Returns
    -------
    (float, float) or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives.
    """
    
    # Try using booster library directly if enabled
    if use_booster: 
        try:
            fp, fn = _compute_fp_fn_with_booster_direct(estimate, inputs, 'sutd', aggregate=aggregate)
            return fp, fn
        except Exception as e:
            # If booster fails, fall through to pure Python implementation
            pass
    
    # Pure Python implementation
    if isinstance(inputs, Tree_with_support):
        return _SUTD1_fp_fn(inputs, estimate)
    elif isinstance(inputs, TreeList_with_support):
        fp_list = []
        fn_list = []
        for i in range(len(inputs)):
            fp_tmp, fn_tmp = _SUTD1_fp_fn(inputs[i], estimate)
            fp_list.append(fp_tmp)
            fn_list.append(fn_tmp)
        if aggregate:
            return np.mean(fp_list), np.mean(fn_list)
        return np.array(fp_list), np.array(fn_list)
    else:
        raise TypeError("input trees must be either Tree_with_support or TreeList_with_support.")


def STD_fp_fn(estimate, inputs, parent_dir=None, aggregate=False, use_booster=True):
    """Computes false positive and false negative of STD.
        
    Parameters
    ----------
    estimate : Tree_with_support
        Estimates. 
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate.
    parent_dir : str, optional
        Path to place intermediate files, by default None (place files in the folder of execution).
    aggregate : bool, optional
        If True, compute aggregated metrics (n+1 booster calls, returns scalars).
        If False (default), compute per-tree metrics (2n booster calls, returns vectors).
    use_booster : bool, optional
        If True (default), attempt to use booster library for faster computation.
        If False, use pure Python implementation.

    Returns
    -------
    (float, float) or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives.
    """
    
    # Try using booster library directly if enabled
    if use_booster:
        try:
            fp, fn = _compute_fp_fn_with_booster_direct(estimate, inputs, 'std', aggregate=aggregate)
            return fp, fn
        except Exception as e:
            # If booster fails, fall through to pure Python implementation
            pass
    
    # Pure Python implementation
    if isinstance(inputs, Tree_with_support):
        return _STD1_fp_fn(inputs, estimate)
    elif isinstance(inputs, TreeList_with_support):
        fp_list = []
        fn_list = []
        for i in range(len(inputs)):
            fp_tmp, fn_tmp = _STD1_fp_fn(inputs[i], estimate)
            fp_list.append(fp_tmp)
            fn_list.append(fn_tmp)
        if aggregate:
            return np.mean(fp_list), np.mean(fn_list)
        return np.array(fp_list), np.array(fn_list)
    else:
        raise TypeError("input trees must be either Tree_with_support or TreeList_with_support.")


def SQD_fp_fn(estimate, inputs, parent_dir=None, aggregate=False):
    """Compute false positives and false negatives of symmetric quartet distance.
    

    Parameters
    ----------
    estimate : Tree_with_support
        Estimates.
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate.
    parent_dir : str, optional
        Path to place intermediate files, by default None (place files in the folder of execution).
        
    Returns
    -------
    (float, float) or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives respectively.
    """
    if isinstance(inputs, Tree_with_support):
        fp, fn = _tqdist_fp_fn(estimate, inputs, parent_dir=parent_dir)
        return fp, fn
    elif isinstance(inputs, TreeList_with_support):
        inputs_string = inputs.as_string("newick", suppress_rooting=True)
        fp, fn = _tqdist_fp_fn_trees_str(estimate, inputs_string, len(inputs), parent_dir=parent_dir)
        if aggregate:
            return np.mean(fp), np.mean(fn)
        return fp, fn
    else:
        raise TypeError("inputs must be an instance of Tree_with_support or TreeList_with_support.")


def SBD_fp_fn(estimate, inputs, parent_dir=None, aggregate=False):
    """Compute false positives and false negatives of symmetric bipartition distance.
    

    Parameters
    ----------
    estimate : Tree_with_support
        Estimates. 
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate.
    parent_dir : str, optional
        Path to place intermediate files, by default None (place files in the folder of execution).
    aggregate : bool, optional
        If True, compute aggregated metrics (returns mean of fp and fn).
        If False (default), compute per-tree metrics (returns vectors).
        
    Returns
    -------
    (float, float) or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives respectively. 
    """
    # Use dendropy implementation only
    if isinstance(inputs, Tree_with_support):
        fp, fn = dendropy.calculate.treecompare.false_positives_and_negatives(inputs, estimate)
        return fp, fn
    elif isinstance(inputs, TreeList_with_support):
        fps = []
        fns = []
        for tree in inputs:
            fp, fn = dendropy.calculate.treecompare.false_positives_and_negatives(tree, estimate)
            fps.append(fp)
            fns.append(fn)
        if aggregate: 
            return np.mean(fps), np.mean(fns)
        return np.array(fps), np.array(fns)
    else:
        raise TypeError("inputs must be an instance of Tree_with_support or TreeList_with_support.")
