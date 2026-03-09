import dendropy
import numpy as np
from bitstring import Bits
from io import StringIO
import warnings
from collections import OrderedDict
import sys
from copy import deepcopy
import itertools
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import to_tree
from scipy.cluster.hierarchy import linkage
import shutil
import subprocess
import os
import tempfile
import random, string
import time
import importlib
import platform
from pathlib import Path
from functools import lru_cache


def _run_subprocess_hidden(cmd, **kwargs):
    """Run subprocess while suppressing console windows on Windows."""
    if os.name == "nt":
        creationflags = kwargs.get("creationflags", 0)
        creationflags |= getattr(subprocess, "CREATE_NO_WINDOW", 0)
        startupinfo = kwargs.get("startupinfo")
        if startupinfo is None:
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        kwargs["creationflags"] = creationflags
        kwargs["startupinfo"] = startupinfo
    return subprocess.run(cmd, **kwargs)


def _import_biophylo():
    """Import Bio.Phylo lazily so optional plotting deps are not forced at startup."""
    try:
        return importlib.import_module("Bio.Phylo")
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Bio.Phylo is required for this operation. Install optional dependency with "
            "`pip install biopython`."
        ) from e


def _resolve_work_dir(parent_dir=None):
    """Return a writable directory for temporary quartet/tqdist files."""
    candidates = []
    if parent_dir:
        candidates.append(parent_dir)
    candidates.append(os.path.join(tempfile.gettempdir(), "Consensus"))
    candidates.append(tempfile.gettempdir())

    for directory in candidates:
        try:
            os.makedirs(directory, exist_ok=True)
            return directory
        except OSError:
            continue
    return "."


def _tmp_nwk_path(parent_dir=None):
    work_dir = _resolve_work_dir(parent_dir)
    return os.path.join(work_dir, "__" + _randomname(10) + ".nwk")


def _platform_folder():
    if sys.platform == "win32":
        return "win"
    if sys.platform == "darwin":
        return "mac"
    return "linux"


def _arch_folder():
    machine = platform.machine().lower()
    if machine in ("x86_64", "amd64"):
        return "x64"
    if machine in ("arm64", "aarch64"):
        return "arm64"
    if machine in ("i386", "i686", "x86"):
        return "x86"
    return machine or "unknown"


def _candidate_tqdist_dirs():
    candidates = []

    env_dir = os.environ.get("CONSENSUS_TQDIST_DIR")
    if env_dir:
        candidates.append(env_dir)

    # PyInstaller bundle locations.
    if getattr(sys, "frozen", False):
        meipass = getattr(sys, "_MEIPASS", "")
        if meipass:
            candidates.append(meipass)
        exe_dir = os.path.dirname(os.path.abspath(sys.executable))
        candidates.append(exe_dir)
        candidates.append(os.path.join(exe_dir, "_internal"))
        if sys.platform == "darwin":
            candidates.append(os.path.normpath(os.path.join(exe_dir, "..", "Frameworks")))
            candidates.append(os.path.normpath(os.path.join(exe_dir, "..", "Resources")))

    # Repository/package layout for development and bundling.
    package_dir = Path(__file__).resolve().parent
    platform_bin_root = package_dir / "bin" / _platform_folder()
    if sys.platform == "darwin":
        candidates.append(os.fspath(platform_bin_root / "universal2"))
    candidates.append(os.fspath(platform_bin_root / _arch_folder()))

    deduped = []
    seen = set()
    for directory in candidates:
        if not directory:
            continue
        directory = os.path.abspath(os.fspath(directory))
        if directory in seen:
            continue
        seen.add(directory)
        deduped.append(directory)
    return deduped


@lru_cache(maxsize=16)
def _resolve_tqdist_executable(executable_name):
    # 1) PATH lookup first.
    found = shutil.which(executable_name)
    if found:
        return found
    if os.name == "nt" and not executable_name.lower().endswith(".exe"):
        found = shutil.which(executable_name + ".exe")
        if found:
            return found

    # 2) Bundled/local lookup.
    name_candidates = [executable_name]
    if os.name == "nt" and not executable_name.lower().endswith(".exe"):
        name_candidates.append(executable_name + ".exe")

    for directory in _candidate_tqdist_dirs():
        for name in name_candidates:
            candidate = os.path.join(directory, name)
            if os.path.isfile(candidate):
                return candidate
    return None


def _require_tqdist_executable(executable_name):
    executable_path = _resolve_tqdist_executable(executable_name)
    if executable_path is None:
        raise FileNotFoundError(
            f"Error: '{executable_name}' not found. Install tqDist on PATH or bundle it under "
            f"'src/Consensus/bin/{_platform_folder()}/{_arch_folder()}'."
        )
    return executable_path



def _get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
    
    Code copied from https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format.

    Parameters
    ----------
    node : _type_
        output of sciply.cluster.hierarchy.to_tree()
    parent_dist : _type_
        output of sciply.cluster.hierarchy.to_tree().dist
    leaf_names : List
        list of leaf names
    newick : str, optional
        leave empty, this variable is used in recursion., by default ''

    Returns
    -------
    str
        tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = _get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = _get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick


def _quartet_resolution(tree, parent_dir = None, normalized=True):
    quartet_dist_exe = _require_tqdist_executable("quartet_dist")
    conname = _tmp_nwk_path(parent_dir)
    conname2 = _tmp_nwk_path(parent_dir)
    tree.write(path = conname, schema="newick", suppress_rooting=True)
    tree.write(path = conname2, schema="newick", suppress_rooting=True)
    p = _run_subprocess_hidden([quartet_dist_exe, "-v", conname, conname2], capture_output=True, text=True)
    keys = ["number_of_leaves", "number_of_all_quartets", "quartet_dist", "normalized_quartet_dist", "number_of_resolved_quartets_agreed", 
            "normalized_number_of_resolved_quartets_agreed", "number_of_unresolved_quartets_agreed", "normalized_number_of_unresolved_quaretets_agreed"]
    items = [float(item) for item in p.stdout.split()]
    q_dict = dict(zip(keys, items))
    for item in [conname,conname2]:
        if os.path.exists(item):
            os.remove(item)
    if normalized:
        return (q_dict["number_of_all_quartets"] - q_dict["number_of_unresolved_quartets_agreed"])/q_dict["number_of_all_quartets"]
    else:
        return (q_dict["number_of_all_quartets"] - q_dict["number_of_unresolved_quartets_agreed"])

def _quartet_resolution_str(tree_string, parent_dir = None, normalized=True):
    quartet_dist_exe = _require_tqdist_executable("quartet_dist")
    conname = _tmp_nwk_path(parent_dir)
    conname2 = _tmp_nwk_path(parent_dir)
    with open(conname, "w") as f:
        f.write(tree_string)
    with open(conname2, "w") as f:
        f.write(tree_string)
    p = _run_subprocess_hidden([quartet_dist_exe, "-v", conname, conname2], capture_output=True, text=True)
    keys = ["number_of_leaves", "number_of_all_quartets", "quartet_dist", "normalized_quartet_dist", "number_of_resolved_quartets_agreed", 
            "normalized_number_of_resolved_quartets_agreed", "number_of_unresolved_quartets_agreed", "normalized_number_of_unresolved_quaretets_agreed"]
    items = [float(item) for item in p.stdout.split()]
    for item in [conname,conname2]:
        if os.path.exists(item):
            os.remove(item)
    num_all_quartets = items[1]
    num_unresolved_quartets_agreed = np.array(items[6::8])
    if normalized:
        return (num_all_quartets - num_unresolved_quartets_agreed)/num_all_quartets
    else:
        return (num_all_quartets - num_unresolved_quartets_agreed)


def _tqdist_fp_fn(estimate, true, parent_dir = None):
    executable_path = _require_tqdist_executable("quartet_dist")
    estimatename = _tmp_nwk_path(parent_dir)
    truename = _tmp_nwk_path(parent_dir)
    
    estimate.write(path = estimatename, schema="newick", suppress_rooting=True)
    true.write(path = truename, schema="newick", suppress_rooting=True)   
    
    p = _run_subprocess_hidden([executable_path, "-v", estimatename, truename], capture_output=True, text=True)
    keys = ["number_of_leaves", "number_of_all_quartets", "quartet_dist", "normalized_quartet_dist", "number_of_resolved_quartets_agreed", 
            "normalized_number_of_resolved_quartets_agreed", "number_of_unresolved_quartets_agreed", "normalized_number_of_unresolved_quaretets_agreed"]
    items = [float(item) for item in p.stdout.split()]
    tqdist_dict = dict(zip(keys, items))
    
    for item in [estimatename,truename]:
        if os.path.exists(item):
            os.remove(item)

    
    # check resolution
    estimate_num_resolved = _quartet_resolution(estimate, parent_dir=parent_dir, normalized=False)
    estimate_num_unresolved = tqdist_dict["number_of_all_quartets"] - estimate_num_resolved
    true_num_resolved = _quartet_resolution(true, parent_dir=parent_dir, normalized=False)
    true_num_unresolved = tqdist_dict["number_of_all_quartets"] - true_num_resolved
    
    # compute number of each component
    unresolved_resolved_fn = estimate_num_unresolved -  tqdist_dict["number_of_unresolved_quartets_agreed"]
    resolved_unresolved_fp = true_num_unresolved -  tqdist_dict["number_of_unresolved_quartets_agreed"]
    resolved_resolved_disagree = estimate_num_resolved - resolved_unresolved_fp - tqdist_dict['number_of_resolved_quartets_agreed']
    assert resolved_resolved_disagree == (true_num_resolved - unresolved_resolved_fn - tqdist_dict['number_of_resolved_quartets_agreed'])
    
    # compute fn and fp
    fn = unresolved_resolved_fn + resolved_resolved_disagree
    fp = resolved_unresolved_fp + resolved_resolved_disagree
    return fp, fn

def _tqdist_fp_fn_trees_str(estimate, input_trees_string, n_trees, parent_dir = None):
    executable_path = _require_tqdist_executable("pairs_quartet_dist")
    estimatename = _tmp_nwk_path(parent_dir)
    truename = _tmp_nwk_path(parent_dir)
    estimate_string = estimate.as_string(schema="newick", suppress_rooting = True)
    with open(estimatename, "w") as f:
        f.write(estimate_string*n_trees)
    #estimate.write(path = estimatename, schema="newick", suppress_rooting=True)
    with open(truename, "w") as f:
        f.write(input_trees_string)
    #input_trees.write(path = truename, schema="newick", suppress_rooting=True)   
    
    p = _run_subprocess_hidden([executable_path, "-v", estimatename, truename], capture_output=True, text=True)
    # keys = ["number_of_leaves", "number_of_all_quartets", "quartet_dist", "normalized_quartet_dist", "number_of_resolved_quartets_agreed", 
    #         "normalized_number_of_resolved_quartets_agreed", "number_of_unresolved_quartets_agreed", "normalized_number_of_unresolved_quaretets_agreed"]
    items = [float(item) for item in p.stdout.split()]
    num_unresolved_quartets_agreed = np.array(items[6::8])
    num_resolved_quartets_agreed = np.array(items[4::8])
    
    assert len(num_unresolved_quartets_agreed) == n_trees
    assert len(num_resolved_quartets_agreed) == n_trees
    
    num_all_quartets = items[1]
    
    
    for item in [estimatename,truename]:
        if os.path.exists(item):
            os.remove(item)
    
    # check resolution
    estimate_num_resolved = _quartet_resolution(estimate, parent_dir=parent_dir, normalized=False)
    estimate_num_unresolved = num_all_quartets - estimate_num_resolved
    true_num_resolved = _quartet_resolution_str(input_trees_string, parent_dir=parent_dir, normalized=False)
    true_num_unresolved = num_all_quartets - true_num_resolved
    
    # compute number of each component
    unresolved_resolved_fn = estimate_num_unresolved -  num_unresolved_quartets_agreed
    resolved_unresolved_fp = true_num_unresolved -  num_unresolved_quartets_agreed
    resolved_resolved_disagree = estimate_num_resolved - resolved_unresolved_fp - num_resolved_quartets_agreed
    assert (resolved_resolved_disagree == (true_num_resolved - unresolved_resolved_fn - num_resolved_quartets_agreed)).all()
    
    # compute fn and fp
    fn = unresolved_resolved_fn + resolved_resolved_disagree
    fp = resolved_unresolved_fp + resolved_resolved_disagree
    return fp, fn
     

def SQD_loss(consensus_tree, input_trees_string, n_trees, normalized=True, parent_dir = None):
    # use tqdist's `pairs_quartet_dist` function
    fp, fn = _tqdist_fp_fn_trees_str(consensus_tree, input_trees_string, n_trees, parent_dir) # fp, fn are ndarrays of length n=len(input_trees)
    loss = np.sum(fp+fn)
    if normalized:
        loss = loss/n_trees
    return loss

def SQD_pruning(consensus_tree, input_trees, parent_dir=None):
    bipartitions = np.array(consensus_tree.encode_bipartitions())
    bipartition_ints = np.array([bipartition.split_as_int() for bipartition in bipartitions])
    internal_edges = consensus_tree.internal_edges(exclude_seed_edge=True)
    internal_bipartitions = np.array([edge.bipartition for edge in internal_edges])
    internal_bipartition_ints = np.array([bipartition.split_as_int() for bipartition in internal_bipartitions])
    internal_edge_dict = {edge.bipartition.split_as_int():edge for edge in internal_edges}


    taxon_namespace = consensus_tree.taxon_namespace
    trees_string = input_trees.as_string("newick", suppress_rooting=True)
    n_trees = len(input_trees)
    current_loss = SQD_loss(consensus_tree, trees_string, n_trees, False, parent_dir)
    next_updates = internal_bipartition_ints
    loss_reduction_dict = {bipar_int: 0 for bipar_int in internal_bipartition_ints}
    iteration = 0
    reduction_list = []
    while True:
        iteration += 1
        st = time.time()
        for bipar_int in next_updates:
            if bipar_int in internal_bipartition_ints:
                mask = bipar_int!=bipartition_ints
                pruned_tree = dendropy.Tree.from_bipartition_encoding(bipartitions[mask], taxon_namespace=taxon_namespace)
                loss = SQD_loss(pruned_tree, trees_string, n_trees, False, parent_dir)
                loss_reduction_dict[bipar_int] = current_loss - loss # if risk reduction happends, this is a positive value
        
        # look for edge that induces maximum risk reduction
        max_key = max(loss_reduction_dict, key=loss_reduction_dict.get)
        max_value = loss_reduction_dict[max_key]
        renew = (max_value > 0) # if renew = True, prune. otherwise, don't prune.
        if not renew:
            break # the current `consensus_tree` is the best

        next_updates = [edge.bipartition.split_as_int() for edge in internal_edge_dict[max_key].adjacent_edges]
        best_mask = max_key!=bipartition_ints
        consensus_tree = dendropy.Tree.from_bipartition_encoding(bipartitions[best_mask], taxon_namespace=taxon_namespace)
        bipartitions = np.array(consensus_tree.encode_bipartitions())
        bipartition_ints = np.array([bipartition.split_as_int() for bipartition in bipartitions])
        internal_edges = consensus_tree.internal_edges(exclude_seed_edge=True)
        internal_bipartitions = np.array([edge.bipartition for edge in internal_edges])
        internal_bipartition_ints = np.array([bipartition.split_as_int() for bipartition in internal_bipartitions])
        internal_edge_dict = {edge.bipartition.split_as_int():edge for edge in internal_edges}
        current_loss = current_loss - max_value
        loss_reduction_dict.pop(max_key)
        ed = time.time()
        reduction_list.append(max_value)
    return consensus_tree, reduction_list


class Tree_with_support(dendropy.Tree):
    """Child class of `dendropy.datamodel.treemodel.Tree`, with support for computing and storing support values.
    
    Attributes:
        x : 1
    """
    def __init__(self, *args, **kwargs):
        branch_support = kwargs.pop("branch_support", None)
        transfer_support = kwargs.pop("transfer_support", None)
        original_taxon_namespace = kwargs.pop("original_taxon_namespace", None)
        
    
        super().__init__(*args, **kwargs)
        
        if original_taxon_namespace == "not provided":
            taxon_labels = [item.label for item in self.taxon_namespace]
            taxon_labels.sort() # sort in alphabetical order.
            new_taxon_namespace = dendropy.TaxonNamespace(taxon_labels) # new taxonnamespaces that is sorted.
            self.migrate_taxon_namespace(new_taxon_namespace)
        
        '''
        Functions to sort taxon_namespace and migrate.
        Note: C implementation uses strcmp() to compare. Same behaviour should be implemented here.ps
        '''
        
        self.branch_support = branch_support
        self.transfer_support = transfer_support
        self.n_taxa = len(self.taxon_namespace)
        self.encode_bipartitions() #some methods (from_bipartitions etc) seem to produce a tree with bad bipartition_encoding, so renew it 

    def compute_branch_support(self, treelist):
        """Compute branch support values of self against input ``treelist`` and update `self.branch_support`
        If Bayesian posterior draws are given as ``treelist``, `branch_support` corresponds to the posterior support.
        If Bootstrap estimates are given as ``treelist``, `branch_support` corresponds to the traditional bootstrap support. 

        Parameters
        ----------
        treelist : TreeList_with_support
            List of trees to evaluate bootstrap support against.

        Returns
        -------
        None
        """
        assert self.taxon_namespace == treelist.taxon_namespace
        branch_support = dict()
        if self.bipartition_encoding is None:
            self.encode_bipartitions()
        for bipar in self.bipartition_encoding:
            key = bipar.split_as_int()
            if key in treelist.edge_dict.keys():
                branch_support[key] = treelist.edge_dict[key]/len(treelist)
            else:
                branch_support[key] = 0
        self.branch_support = branch_support
        return None
    
    
    def compute_transfer_support(self, treelist):
        """Compute transfer_support support of self against input ``treelist`` and update `self.transfer_support`

        Parameters
        ----------
        treelist : TreeList_with_support
            List of trees to evaluate transfer_support support against.

        Returns
        -------
        None
        """
        assert self.taxon_namespace == treelist.taxon_namespace
        if self.bipartition_encoding  is None:
            self.encode_bipartitions()
        bipartitions = self.bipartition_encoding
        bipartitions_ints = [bipar.split_as_int() for bipar in bipartitions]
        transfer_support_list = transfer_support(bipartitions, treelist)
        self.transfer_support = dict(zip(bipartitions_ints, transfer_support_list))
        return None

    def clone(self, depth=1):
        """clone

        Parameters
        ----------
        depth : int, optional
            _description_, by default 1

        Returns
        -------
        _type_
            _description_
        """
        cloned_tree = super().clone(depth=1)
        return Tree_with_support(cloned_tree, branch_support=self.branch_support, transfer_support = self.transfer_support)
        

    def _make_Bio_tree(self, support = None):
        """Helper function to convert self to Bio tree.

        Parameters
        ----------
        support : dict or None, optional
            Support values to provide to Bio tree.   
            If dict, it should have integer representing each bipartition as key and support values as value.
            If None, no support values are provided to Bio tree.
            By default None

        Returns
        -------
        Bio.Phylo.Newick.Tree
            Newly created Bio tree.
        """
        phylo = _import_biophylo()
        tree_Bio = phylo.read(file=StringIO(self.as_string(schema="newick")), format="newick")
        if support is not None:
            for clade in tree_Bio.find_clades(order="postorder"):
                terminals = [item.name for item in clade.get_terminals()]
                taxonnames_array = np.array([item.label for item in self.taxon_namespace])
                digits = [ np.where( item == taxonnames_array )[0][0] for item in terminals ]
                clade_bool = [False for i in range(self.n_taxa)]
                for digit in digits:
                    clade_bool[-(digit+1)] = True
                clade_bit = Bits(clade_bool) # bitstring corresponding to the clade
                if int(clade_bit.bin[-1]) == 1:
                    clade_bit = (~clade_bit) # convert to least significant 0 format
                if clade_bit.uint in support.keys():
                    clade_support = support[clade_bit.uint]
                    clade.confidence = clade_support
        return tree_Bio      
    def branch_resolution(self):
        """Returns branch resolution of self.
        
        The branch resolution is defined by:
        
        .. math::

            \\frac{\\text{# internal edges in self}}{\\text{self.n_taxa} - 3}

            

        Returns
        -------
        float
            Branch resolution
        """
        possible_num = self.n_taxa - 3 # dendropy Tree stores n external edges + 1 external edge connecting to seed_node. There are n-3 internal edges.
        internal_num = len(self.internal_edges(exclude_seed_edge=True))
        return internal_num / possible_num
    
    def quartet_resolution(self):
        """Returns quartet resolution of self.
        
        The quartet resolution is defined by:
        
        .. math::

            \\frac{\\text{# resolved quartets in self}}{\\binom{n}{4}}

            

        Returns
        -------
        float
            Quartet resolution
        """
        return _quartet_resolution(self)

    def STD_greedy_pruning(self, treelist, normalized=True):
        """Apply greedy pruning algorithm w.r.t. STD loss.

        Parameters
        ----------
        treelist : TreeList_with_support
            Input trees to evaluate STD loss.
        normalized : bool, optional
            Whether to normalized STD loss, by default True

        Returns
        -------
        Tree_with_support
            Consensus tree after applying greedy pruning
        """
        self_copy = self.clone(depth=1)
        srp = std_risk_prune(self_copy, treelist, normalized)
        srp.greedy_pruning()
        return srp.current_tree
    
    def SQD_greedy_pruning(self, treelist, parent_dir=None):
        """Apply greedy pruning algorithm w.r.t. SQD loss.

        Parameters
        ----------
        treelist : TreeList_with_support
            Input trees to evaluate SQD loss.

        Returns
        -------
        Tree_with_support
            Consensus tree after applying greedy pruning
        """
        self_copy = self.clone(depth=1)
        res, reduction_list = SQD_pruning(self_copy, treelist, parent_dir)
        return res
    
    def BS_prune(self, treelist, threshold=0.5):
        """Prune all edges that has branch support <= threshold.

        Parameters
        ----------
        treelist : TreeList_with_support
            Input trees to compute branch support
        threshold : float, optional
            Threshold of branch support, by default 0.5

        Returns
        -------
        Tree_with_support
            Thresholded consensus
        """
        if self.branch_support is None:
            self.compute_branch_support(treelist)
            
        assert(self.taxon_namespace == treelist.taxon_namespace)
        
        bipars = []; branch_support = dict(); transfer_support = dict()
        for node in self.postorder_node_iter():
            splitint = node.bipartition.split_as_int()
            if self.branch_support[splitint] > threshold:
                bipars.append(node.bipartition)
                branch_support[splitint] = self.branch_support[splitint]
                if self.transfer_support is not None:
                    transfer_support[splitint] = self.transfer_support[splitint]
        thresholded = Tree_with_support(dendropy.Tree.from_bipartition_encoding(bipars, self.taxon_namespace), branch_support = branch_support, transfer_support = transfer_support)
        return thresholded
         
    def TS_prune(self, treelist, threshold=0.7):
        if self.transfer_support is None:
            self.compute_transfer_support(treelist)
        
        assert(self.taxon_namespace == treelist.taxon_namespace)
        
        bipars = []; branch_support = dict(); transfer_support = dict()
        for node in self.postorder_node_iter():
            splitint = node.bipartition.split_as_int()
            if self.transfer_support[splitint] > threshold:
                bipars.append(node.bipartition)
                if self.branch_support is not None:
                    branch_support[splitint] = self.branch_support[splitint]
                transfer_support[splitint] = self.transfer_support[splitint]
        thresholded = Tree_with_support(dendropy.Tree.from_bipartition_encoding(bipars, self.taxon_namespace), branch_support = branch_support, transfer_support = transfer_support)
        return thresholded
        
        
    @classmethod
    def get(cls, **kwargs):
        """get

        Returns
        -------
        _type_
            _description_
        """
        tree = dendropy.Tree.get(**kwargs)
        branch_support = kwargs.pop("branch_support", None)
        transfer_support = kwargs.pop("transfer_support", None)
        original_taxon_namespace = kwargs.get("taxon_namespace", "not provided") # check if taxon_namespace is passed
        return cls(tree, branch_support = branch_support, transfer_support = transfer_support, original_taxon_namespace = original_taxon_namespace)


def _create_refinfo(bipartitions, n_taxa):
    """Helper function for transfer_support

    Parameters
    ----------
    bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        _description_
    n_taxa : int
        _description_

    Returns
    -------
    dict
        bipar_int as key an tuple of (bitstr, p) as value.
    """
    
    refinfo = dict()
    for bipartition in bipartitions:
        edge_bitstr = Bits(uint = bipartition.split_as_int(), length=n_taxa)
        counts = (edge_bitstr.count(0), edge_bitstr.count(1))
        if counts[0] < counts[1]:
            # zero lready assigned to light side
            bitstr = edge_bitstr
            p = counts[0]
        else:
            # reverse bitstring
            bitstr = (~edge_bitstr)
            p = counts[1]
        refinfo[bipartition.split_as_int()] = (bitstr, p)
    # O( len(bipartitions) *  n_taxa )
    return refinfo

def _minDist(refinfo_b, tree):
    """Helper minDist function for computing transfer_support.

    Parameters
    ----------
    refinfo_b : dict
        refinfo of bipartition b created by _create_refinfo function.
    tree : Tree_with_support or dendropy.datamodel.treemodel.Tree
        Tree to compute support of the bipartition
        
    Notes
    -----
    This function expects input bipartition (``refinfo_b``) to be NOT included in ``tree``.

    Returns
    -------
    int
        Hamming distance.
    """
    # bipartition: Bipartition object
    # tree: Tree_with_support or Tree object
    p = refinfo_b[1]
    b_bitstr = refinfo_b[0]

    d = p - 1

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
            if actDist < d:
                d = actDist
                if d == 1:
                    return d
    
    return d



def _minDist_and_match2(refinfo_b, tree):
    # bipartition: Bipartition object
    # tree: Tree_with_support or Tree object
    p = refinfo_b[1]
    b_bitstr = refinfo_b[-0]

    d = p - 1
    d2 = p - 1
    matched=-1
    second_match = -1

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
            if actDist < d:
                # renew second match
                d2 = d
                second_match = matched
                # renew first match
                d = actDist
                matched = node.bipartition.split_as_int()
            elif actDist < d2:
                # renew only second match
                d2 = actDist
                second_match = node.bipartition.split_as_int()
    return d, matched, d2, second_match



def _transfer_support_and_match2(bipartitions, tree:Tree_with_support):
    n_taxa = len(tree.taxon_namespace)
    refinfo = _create_refinfo(bipartitions, n_taxa)
    n_bipartitions = len(bipartitions)
    totalSupport = np.zeros(n_bipartitions)
    secondSupport = np.zeros(n_bipartitions)

    if tree.bipartition_encoding is None:
        tree.encode_bipartitions()
    tree_bipartition_ints = [bipartition.split_as_int() for bipartition in tree.bipartition_encoding]
    matched_list = []
    second_match_list = []
    for i, bipartition in enumerate(bipartitions):
        bipar_int = bipartition.split_as_int()
        if refinfo[bipar_int][1]<=1:
            support=1
            matched=-1
            second_support=1
            second_match=-1
        elif bipar_int not in tree_bipartition_ints and refinfo[bipar_int][1] == 2:
            support = 0
            matched = -1
            second_support = 0
            second_match = -1
        else:
            dis, matched, dis2, second_match = _minDist_and_match2(refinfo[bipar_int], tree)
            support = 1 - dis / (refinfo[bipar_int][1] - 1)
            second_support = 1 - dis2 / (refinfo[bipar_int][1] - 1)
        totalSupport[i] = support
        secondSupport[i] = second_support
        matched_list.append(matched)
        second_match_list.append(second_match)
    
    return totalSupport, matched_list, secondSupport, second_match_list

def _unnormalized_transfer_support_and_match2(bipartitions, tree:Tree_with_support):
    n_taxa = len(tree.taxon_namespace)
    refinfo = _create_refinfo(bipartitions, n_taxa)
    n_bipartitions = len(bipartitions)
    totalDistance = np.zeros(n_bipartitions)
    secondDistance = np.zeros(n_bipartitions)

    if tree.bipartition_encoding is None:
        tree.encode_bipartitions()
    tree_bipartition_ints = [bipartition.split_as_int() for bipartition in tree.bipartition_encoding]
    matched_list = []
    second_match_list = []
    for i, bipartition in enumerate(bipartitions):
        bipar_int = bipartition.split_as_int()
        if refinfo[bipar_int][1]<=1:
            distance=0
            matched=-1
            second_distance=1
            second_match=-1
        elif bipar_int not in tree_bipartition_ints and refinfo[bipar_int][1] == 2:
            distance = 1
            matched = -1
            second_distance = 0
            second_match = -1
        else:
            distance, matched, second_distance, second_match = _minDist_and_match2(refinfo[bipar_int], tree)
            #support = 1 - dis / (refinfo[bipar_int][1] - 1)
            #second_support = 1 - dis2 / (refinfo[bipar_int][1] - 1)
        totalDistance[i] = distance
        secondDistance[i] = second_distance
        matched_list.append(matched)
        second_match_list.append(second_match)
    
    return totalDistance, matched_list, secondDistance, second_match_list




class std_risk_prune:

    def __init__(self, initial_tree, trees, normalized=True):
        self.current_tree = initial_tree
        self.trees = trees
        self.n_trees = len(trees)
        self.normalized=normalized
        assert self.current_tree.taxon_namespace == self.trees.taxon_namespace
        if self.current_tree.bipartition_encoding  is None:
            self.current_tree.encode_bipartitions()
        bipartitions = self.current_tree.bipartition_encoding
        bipartition_ints = [bipar.split_as_int() for bipar in bipartitions]
        if normalized:
            transfer_support_list = transfer_support(bipartitions, self.trees)
        else:
            transfer_support_list = unnormalized_transfer_support(bipartitions, self.trees)
        self.transfer_support = dict(zip(bipartition_ints, transfer_support_list))
        if normalized:
            self.fp = np.sum([1-support for key, support in self.transfer_support.items()]) * self.n_trees
        else:
            self.fp = np.sum([support for key, support in self.transfer_support.items()]) * self.n_trees

        # use edge dict
        self.edge_dict = self.trees.edge_dict
        self.trees_bipartition_ints = list(self.edge_dict.keys())
        self.trees_bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.current_tree.n_taxa - 1) for item in self.trees_bipartition_ints]
        if normalized:
            supports, matched_list, second_supports, second_matched_list = _transfer_support_and_match2(self.trees_bipartitions, self.current_tree)
        else:
            supports, matched_list, second_supports, second_matched_list = _unnormalized_transfer_support_and_match2(self.trees_bipartitions, self.current_tree)

        # initialize match_dict, reverse_match_dict, 
        self.match_dict = dict(zip(self.transfer_support.keys(), [list() for i in range(len(self.transfer_support))]))
        self.second_match_dict = dict(zip(self.transfer_support.keys(), [list() for i in range(len(self.transfer_support))]))
        self.reverse_match_dict = dict(zip(self.trees_bipartition_ints, [list() for i in range(len(self.trees_bipartition_ints))]))
        self.match_dict[-1] = list()
        self.second_match_dict[-1] = list()
        for i in range(len(self.trees_bipartition_ints)):
            bipar_int = self.trees_bipartition_ints[i]
            matched = matched_list[i]
            second_matched = second_matched_list[i]
            self.match_dict[matched].append(bipar_int)
            self.second_match_dict[second_matched].append(bipar_int)
            self.reverse_match_dict[bipar_int] = [matched, second_matched]

        # intialize rsupp_dict, second_supp_dict
        self.rsupp_dict = dict(zip(self.trees_bipartition_ints, supports))
        self.second_supp_dict = dict(zip(self.trees_bipartition_ints, second_supports))
        # intialize fn_dict
        if normalized:
            self.fn = np.sum([ (1-item) * self.edge_dict[key] for key, item in self.rsupp_dict.items() ])
        else:
            self.fn = np.sum([ item * self.edge_dict[key] for key, item in self.rsupp_dict.items() ])

    def find_risk_reduction(self, prune_bipar_int):
        if self.normalized:
            fp_reduction = (1-self.transfer_support[prune_bipar_int]) * self.n_trees #positive value
        else:
            fp_reduction = self.transfer_support[prune_bipar_int] * self.n_trees
        # compute fn_increase
        matched_bipars = self.match_dict[prune_bipar_int]
        
        # compute support difference
        if self.normalized:
            fn_diff = np.sum([ (self.rsupp_dict[matched_bipars[i]] - self.second_supp_dict[matched_bipars[i]]) * self.edge_dict[matched_bipars[i]] 
                for i in range(len(matched_bipars))]) # positive value
        else:
            fn_diff = np.sum([ (self.second_supp_dict[matched_bipars[i]] - self.rsupp_dict[matched_bipars[i]]) * self.edge_dict[matched_bipars[i]] 
                for i in range(len(matched_bipars))]) # positive value

        return fp_reduction - fn_diff
    
    def prune(self, prune_bipar_int):
        if self.normalized:
            fp_reduction = (1-self.transfer_support[prune_bipar_int]) * self.n_trees
        else:
            fp_reduction = self.transfer_support[prune_bipar_int] * self.n_trees
        matched_bipars = self.match_dict[prune_bipar_int]
        rematch_bipartitions = {item.split_as_int():item for item in self.trees_bipartitions if item.split_as_int() in matched_bipars}
        rematch_bipartitions = [rematch_bipartitions[item] for item in matched_bipars]
        second_matched_bipars = self.second_match_dict[prune_bipar_int]
        second_rematch_bipartitions = {item.split_as_int():item for item in self.trees_bipartitions if item.split_as_int() in second_matched_bipars}
        second_rematch_bipartitions = [second_rematch_bipartitions[item] for item in second_matched_bipars]
        

        if self.normalized:
            fn_diff = np.sum([ (self.rsupp_dict[matched_bipars[i]] - self.second_supp_dict[matched_bipars[i]]) * self.edge_dict[matched_bipars[i]] 
                for i in range(len(matched_bipars))]) # positive value
        else:
            fn_diff = np.sum([ (self.second_supp_dict[matched_bipars[i]] - self.rsupp_dict[matched_bipars[i]]) * self.edge_dict[matched_bipars[i]] 
                for i in range(len(matched_bipars))]) # positive value
        new_bipartitions = [deepcopy(bipar) for bipar in self.current_tree.bipartition_encoding if bipar.split_as_int()!=prune_bipar_int]
        new_tree = Tree_with_support(dendropy.Tree.from_bipartition_encoding(new_bipartitions, self.current_tree.taxon_namespace),
                                        transfer_support = self.transfer_support)
        self.current_tree = new_tree
        self.current_tree.encode_bipartitions()
        self.transfer_support.pop(prune_bipar_int)
        self.fp = self.fp - fp_reduction

        # run matching for rematch_bipartitions.
        if self.normalized:
            supports, matched_list, second_supports, second_matched_list = _transfer_support_and_match2(rematch_bipartitions, self.current_tree)
        else:
            supports, matched_list, second_supports, second_matched_list = _unnormalized_transfer_support_and_match2(rematch_bipartitions, self.current_tree)
    
        # renew match dict, second_match_dict, and reverse_match_dict and second_match_dict
        self.match_dict.pop(prune_bipar_int)
        for i in range(len(matched_bipars)):
            bipar_int = matched_bipars[i]
            original_first_match = self.reverse_match_dict[bipar_int][0]
            original_second_match = self.reverse_match_dict[bipar_int][1]
            self.reverse_match_dict[bipar_int] = [matched_list[i], second_matched_list[i]]
            self.match_dict[matched_list[i]].append(bipar_int) # set the newly matched one
            self.second_match_dict[original_second_match].remove(bipar_int) # remove from original second match bipar_int
            self.second_match_dict[second_matched_list[i]].append(bipar_int)
            # renew supports
            if np.abs(self.second_supp_dict[bipar_int] - supports[i]) > 1e-5:
                print("bad second support", self.second_supp_dict[bipar_int], supports[i])
                print(bipar_int)
                print(matched_list[i], second_matched_list[i], original_first_match, original_second_match)
            self.rsupp_dict[bipar_int] = supports[i]
            self.second_supp_dict[bipar_int] = second_supports[i]

        if self.normalized:
            supports, matched_list, second_supports, second_matched_list = _transfer_support_and_match2(second_rematch_bipartitions, self.current_tree)
        else:
            supports, matched_list, second_supports, second_matched_list = _unnormalized_transfer_support_and_match2(second_rematch_bipartitions, self.current_tree)
        self.second_match_dict.pop(prune_bipar_int)
        for i in range(len(second_matched_bipars)):
            bipar_int = second_matched_bipars[i]
            first_match = matched_list[i]; second_match = second_matched_list[i]
            second_support = second_supports[i]
            if first_match != self.reverse_match_dict[bipar_int][0]:
                print("different first match")
                second_match = first_match # this way we don't have to modify the first match
                first_match = self.reverse_match_dict[bipar_int][0]
                second_support = self.rsupp_dict[bipar_int]
            self.reverse_match_dict[bipar_int][1] = second_match
            self.second_match_dict[second_match].append(bipar_int) # add new second match
            self.second_supp_dict[bipar_int] = second_support

        self.fn = self.fn + fn_diff
    

    def greedy_pruning(self):
        while True:
            bipartition_ints = list(self.transfer_support.keys())
            # find risk reduction
            risk_reductions = [self.find_risk_reduction(bipar_int) for bipar_int in bipartition_ints]
            max_risk_reduction_index = np.argmax(risk_reductions)
            if risk_reductions[max_risk_reduction_index] <= 0:
                # risk reduction does not happen, so break
                break
            # edge with max risk reduction will be removed.
            prune_bipar = bipartition_ints[max_risk_reduction_index]
            self.prune(prune_bipar)

def _randomname(n):
   return ''.join(random.choices(string.ascii_letters + string.digits, k=n))

def _compatible(bits_a : Bits, bits_b : Bits):
    """Check if two bitsting represented clades are compatible.

    Parameters
    ----------
    bits_a : Bits
        Bitstring representation of first bipartition.
    bits_b : Bits
        Bitstring representation of second bipartition.

    Returns
    -------
    int
        1 indicates compatible, 0 indicates incompatible.
    """
    count1 = ( bits_a & bits_b ).count(1)
    count2 = ((~bits_a) & bits_b).count(1)
    count3 = (bits_a & (~bits_b)).count(1)
    count4 = ((~bits_a) & (~bits_b)).count(1)
    if np.min([count1, count2, count3, count4]) == 0:
        return 1 # compatible
    else:
        return 0 # incompatible

def _Rstar(trees,root_name : str):
    """Computes R* consensus tree using Bio.Phylo package.

    Parameters
    ----------
    trees : Iterable of Bio.Phylo.Newick.Tree
        Input trees to compute R* consensus.
    root_name : str
        Name of the root node picked.

    Returns
    -------
    Bio.Phylo.Newick.Tree
        R* consensus tree
    """
    # root: clade instance of the leaf we make a root.
    # currently very inefficient and messy
    rooted_trees = []
    for tree in trees:
        rooted_tree = deepcopy(tree)
        terminals = rooted_tree.get_terminals()
        root_ind = np.where(np.array([item.name for item in terminals]) == root_name)[0][0]
        root = terminals[root_ind]
        rooted_tree.root_with_outgroup(root)
        rooted_tree.prune(root)
        rooted_trees.append(rooted_tree)
    
    all_taxa_list = [item.name for item in rooted_trees[0].get_terminals()]
    all_taxa_set = set(all_taxa_list)

    
    # Initialization: dict_RT_count O(n^3)
    dict_RT_count = {}
    for combi in itertools.combinations(all_taxa_list, 3):
        dict_RT_count[ ( frozenset((combi[0], combi[1])) ,combi[2] ) ] = 0
        dict_RT_count[ ( frozenset((combi[1], combi[2])) ,combi[0] ) ] = 0
        dict_RT_count[ ( frozenset((combi[2], combi[0])) ,combi[1] ) ] = 0
    # For each tree, we traverse by DFS and add counts of rooted triples O(Kn^3)
    for tree in rooted_trees:
        for clade in tree.find_clades(order="postorder"):
            nonterminals = clade.get_nonterminals()
            if len(nonterminals) == 0:
                continue # leaf clade
            
            # obtain the list of subtrees
            subtree_list = clade.clades
            # obtain list of list of taxon in each subtree
            subtree_taxa_list = [ [item.name for item in  tr.get_terminals()] for tr in subtree_list]
            # set of taxon in this current clade and the other side of bipartition
            taxa_set = {item.name for item in clade.get_terminals()}
            taxa_other = all_taxa_set - taxa_set
            # if this corresponds to root, skip
            if len(taxa_other)==0:
                continue

            for st_a, st_b in itertools.combinations(subtree_taxa_list, 2):
                for taxon_a in st_a:
                    for taxon_b in st_b:
                        for taxon_c in taxa_other:
                            # add count of (a,b)|c
                            dict_RT_count[(frozenset({taxon_a, taxon_b}), taxon_c)] += 1
    # Initialization of Rmaj
    Rmaj = {}
    for combi in itertools.combinations(all_taxa_set, 2):
        Rmaj[frozenset({combi[0], combi[1]})] = list()
    # again iterate over all rooted triples
    for combi in itertools.combinations(all_taxa_set, 3):
        count_0 = dict_RT_count[(frozenset({combi[0], combi[1]}) , combi[2])]
        count_1 = dict_RT_count[(frozenset({combi[1], combi[2]}) , combi[0])]
        count_2 = dict_RT_count[(frozenset({combi[2], combi[0]}) , combi[1])]
        counts = [count_0, count_1, count_2]
        argmax_ind = np.where(counts == np.max(counts))[0]
        if len(argmax_ind) > 1:
            # multiple candidates, so add nothing to Rmaj
            continue
        # else, add to Rmaj
        Rmaj[ frozenset( { combi[argmax_ind[0]], combi[(argmax_ind[0] + 1)%3 ] } ) ].append(combi[(argmax_ind[0] + 2)%3])
    n_taxa = len(all_taxa_list)
    s_maj = np.diag([-1 for i in range(n_taxa)])# initialization
    for i in range (n_taxa):
        for j in range(i):
            s_maj[i,j] = s_maj[j,i] = len( Rmaj[ frozenset( { all_taxa_list[i], all_taxa_list[j] } ) ] )
    condensed_similarity = squareform(s_maj, checks=False) # condensed SIMILARITY matrix
    condensed_dissimilarity = - condensed_similarity + n_taxa # transforming into nonnegative value as well
    linkage_matrix = linkage(condensed_dissimilarity, method = "single", metric="IGNORED")
    slink_tree = to_tree(linkage_matrix)
    slink_newick = _get_newick(slink_tree, slink_tree.dist, leaf_names=all_taxa_list)
    phylo = _import_biophylo()
    slink = phylo.read(StringIO(slink_newick), "newick")
    #unit branch length version
    slink_unit = deepcopy(slink)
    for clade in slink_unit.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = 1
    collapse_clades = []
    for clade in slink_unit.find_clades(order = "postorder"):
        if clade.is_terminal():
            # leaf node
            continue
        terminals = {item.name for item in clade.get_terminals()}
        outsiders = all_taxa_set - terminals
        if len(outsiders) ==0:
            # trivial clade
            continue
        end = False
        for a, b in itertools.combinations(terminals, 2):
            for c in outsiders:
                if c not in Rmaj[ frozenset( {a, b} ) ]:
                    collapse_clades.append(clade)
                    end = True
                    break
            if end:
                break
    slink_unit.collapse_all(lambda c: c in collapse_clades)
    return slink_unit


class TreeList_with_support(dendropy.TreeList):
    """Child class of dendropy.datamodel.treecollectionmodel.TreeList, with support for several consensus constructions and support values.
    """
    # EXPECTS TreeList with the same leaves. 
    # EXPECTS TreeList of unrooted trees

    def _make_edge_dict(self):
        edge_dict = OrderedDict()
        for tree in self:
            dict_keys = edge_dict.keys() # no duplicate keys in one tree
            for bipartition in tree.encode_bipartitions():
                key = bipartition.split_as_int()
                if key in dict_keys:
                    edge_dict[key] += 1
                else:
                    edge_dict[key] = 1
        return edge_dict
    
    def _make_all_transfer_support_dict(self):
        keys = list(self.edge_dict.keys())
        bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.n_taxa - 1, compile_bipartition=True) for item in keys]
        transfer_support_list = self.transfer_support(bipartitions)
        bipartition_ints = [bipar.split_as_int() for bipar in bipartitions]
        self.all_transfer_support_dict =  dict(zip(bipartition_ints, transfer_support_list))
        return None
    
    def __init__(self, *args, **kwargs):
        edge_dict = kwargs.pop("edge_dict", None)
        self.all_transfer_support_dict =  kwargs.pop("all_transfer_support_dict", None)
        if len(args) ==  1 and isinstance(args[0], dendropy.TreeList):
            # just reference it
            original_taxon_namespace = kwargs.pop("original_taxon_namespace", None)
            self.__dict__.update(args[0].__dict__)
        else:
            original_taxon_namespace = kwargs.pop("original_taxon_namespace", None)
            super().__init__(*args, **kwargs)
        if original_taxon_namespace is None:
            taxon_labels = [item.label for item in self.taxon_namespace]
            taxon_labels.sort() # sort in alphabetical order.
            new_taxon_namespace = dendropy.TaxonNamespace(taxon_labels) # new taxonnamespaces that is sorted.
            self.migrate_taxon_namespace(new_taxon_namespace)
        trees_init = True
        if len(args) == 1:
            try:
                if len(args[0]) > 0 and all(isinstance(x, Tree_with_support) for x in args[0]):
                    # the given trees are Tree_with_support, so initialize it with this (not to lose additional info)
                    self._trees = list(args[0])
                    trees_init = False
            except:
                pass
        if trees_init:
            self._trees = [Tree_with_support(tree) for tree in self._trees]
        self.edge_dict = edge_dict
        self.n_trees = len(self)
        self.n_taxa = len(self[0].leaf_nodes())
        if self.edge_dict is None:
            self.edge_dict = self._make_edge_dict()
    
    @classmethod
    def get(cls, **kwargs):
        """get

        Returns
        -------
        _type_
            _description_
        """
        edge_dict = kwargs.pop("edge_dict", None)
        tree = dendropy.TreeList.get(**kwargs)
        taxon_namespace = kwargs.get("taxon_namespace", None)
        return cls(tree, edge_dict = edge_dict, original_taxon_namespace = taxon_namespace)



    def branch_support(self, bipartition):
        """Compute bootstrap support of a bipartition against self.

        Parameters
        ----------
        bipartition : `dendropy.datamodel.treemodel.Bipartition`
            Bipartition to evaluate bootstrap support

        Returns
        -------
        float
            Bootstrap support.
        """
        # computes bootstrap support of an edge.
        # bipartition: Bipartition object

        key = bipartition.split_as_int()
        if key in self.edge_dict.keys():
            return self.edge_dict[key]/self.n_trees
        return 0
    
    def compute_all_transfer_support(self):
        """Compute all transfer_support of all bipartitions present in the treelist against self and save it in self.all_transfer_support_dict.

        Returns
        -------
        dict
            Dictionary that has integers representing bipartitions as keys and transfer_support support as values. 
        """
        self._make_all_transfer_support_dict()
        return self.all_transfer_support_dict
    
    def transfer_support(self, bipartitions):
        """Compute transfer_support support of bipartitions against self.

        Parameters
        ----------
        bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        Bipartitions to evaluate transfer_support support.

        Returns
        -------
        numpy.ndarray
            Array of transfer_support support
        """
        return transfer_support(bipartitions, self)
    
    def majority_rule_consensus(self):
        """Computes majority rule consensus.

        Returns
        -------
        Tree_with_support
            Majority rule consensus tree.
        """
        tree_leafset_bitmask = 2**self.n_taxa - 1 # corresponds to "1111......1"
        bipartitions = []; bipar_ints = []; bipar_bs = []
        for k, v in self.edge_dict.items():
            if v > self.n_trees/2:
                bipartitions.append(dendropy.Bipartition(bitmask = k, leafset_bitmask = k, tree_leafset_bitmask = tree_leafset_bitmask, is_rooted=False))
                bipar_ints.append(k); bipar_bs.append(v/self.n_trees)
        majority_tree = dendropy.Tree.from_bipartition_encoding(bipartitions, self.taxon_namespace)
        branch_support = dict(zip(bipar_ints, bipar_bs))
        # majority_tree.branch_support = branch_support
        return Tree_with_support(majority_tree, branch_support = branch_support, taxon_namespace = self.taxon_namespace)

    def MCC(self):
        """Computes Maximum Clade Credibility (MCC) tree.

        Returns
        -------
        Tree_with_support
            Maximum Clade Credibility (MCC) tree.
        """
        max_ind = 0
        maxcredi = 0
        for index, tree in enumerate(self):
            internals = tree.internal_nodes(exclude_seed_node=True)
            credi = 1
            for node in internals:
                credi *= self.edge_dict[node.bipartition.split_as_int()]
            if credi > maxcredi:
                # renew max_ind and maxcredi
                max_ind = index
                maxcredi = credi
        return self[max_ind].clone() # return cloned tree with depth 1
            
    
    def MAP(self):
        """Returns the most appearing tree topology.

        Returns
        -------
        List
            List of the most appearing tree topologies (has length > 1 when more than two topologies have the same maximum count).
        
        int
            The number of occurrences of MAP trees.
        """
        dict_topologies = dict()
        for tree_idx in range(len(self)):
            tree = self[tree_idx]
            bipar_set = frozenset({node.bipartition.split_as_int() for node in tree.postorder_internal_node_iter(exclude_seed_node=True)})
            if bipar_set in dict_topologies.keys():
                dict_topologies[bipar_set].append(tree_idx)
            else:
                dict_topologies[bipar_set] = [tree_idx]
        len_list = [len(item) for item in dict_topologies.values()]
        max_len = np.max(len_list)
        max_len_topology_indices = [item for item in dict_topologies.values() if len(item)==max_len]
        repr_trees = []
        for item in max_len_topology_indices:
            repr_tree = self[item[0]].clone(depth=1)
            len_dict = dict()
            for edge in self[item[0]].postorder_internal_edge_iter(exclude_seed_edge=True):
                len_dict[edge._bipartition.split_as_int()] = edge.length/max_len
            for treeind in item[1:]:
                tree = self[treeind]
                for edge in tree.postorder_internal_edge_iter(exclude_seed_edge=True):
                    len_dict[edge._bipartition.split_as_int()] += edge.length/max_len
            for edge in repr_tree.postorder_internal_edge_iter(exclude_seed_edge=True):
                edge.length = len_dict[edge._bipartition.split_as_int()]
            repr_trees.append(repr_tree)
        return repr_trees, max_len

    def Rstar_consensus(self, root_index = None, root_name = None):
        """Compute R* consensus.

        Parameters
        ----------
        root_index : int, optional
            The index of root in self.taxon_namespace, by default None.
            You need to specify either ``root_index`` or ``root_name``.
        root_name : _type_, optional
            The name of root taxa, by default None
            You need to specify either ``root_index`` or ``root_name``.

        Returns
        -------
        Tree_with_support
            R* consensus tree.
        """

        # specify either root_index or root_name
        # currently very inefficient and messy

        if root_index is not None:
            if root_name is not None:
                raise ValueError(
                    "Both root_index and root_name are specified. Specify only one of them."
                )
            root_name = self.taxon_namespace[root_index].label
        elif root_name is None:
            raise ValueError("You need to specify either root_index or root_name.")
        
        taxon_idx = np.where(np.array( [item.label for item in self.taxon_namespace] ) == root_name)[0][0]
        taxon = self.taxon_namespace[int(taxon_idx)]

        phylo = _import_biophylo()
        trees = [phylo.read(file=StringIO(tree.as_string(schema="newick")), format="newick") for tree in self]
        
        tree = Tree_with_support.get(data = _Rstar(trees, root_name).__format__("newick"), schema="newick", taxon_namespace = self.taxon_namespace, rooting="force-rooted")
        node = dendropy.Node(taxon = taxon, edge_length = 1)
        tree.seed_node.add_child(node)
        tree.is_rooted = False
        return tree

def transfer_support(bipartitions, treelist):
    """Compute transfer_support of given set of bipartitions against input trees.

    Parameters
    ----------
    bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        Bipartitions to evaluate transfer support.
    treelist : TreeList_with_support or Iterable of either Tree_with_support or 
        Input trees to evaluate transfer support against.
    
    Notes
    -----
    Trees that produce input bipartitions and treelist need to have the same taxon_namespace.
    
    Otherwise, the produced transfer support will not be a valid one.
    
    Please check this manually, as we don't check it inside the function.

    Returns
    -------
    numpy.ndarray
        Array of transfer support. 
    """
    
    n_taxa = len(treelist[0].taxon_namespace)
    refinfo = _create_refinfo(bipartitions, n_taxa)
    n_bipartitions = len(bipartitions)
    totalSupport = np.zeros(n_bipartitions)

    for tree in treelist:
        if tree.bipartition_encoding is None:
            tree.encode_bipartitions()
        tree_bipartition_ints = [bipartition.split_as_int() for bipartition in tree.bipartition_encoding]
        for i, bipartition in enumerate(bipartitions):
            bipar_int = bipartition.split_as_int()
            if bipar_int in tree_bipartition_ints:
                support = 1
            elif refinfo[bipar_int][1] == 2:
                support  = 0
            else:
                support = 1 - _minDist(refinfo[bipar_int], tree) / (refinfo[bipar_int][1] - 1)
            totalSupport[i] += support
    
    totalSupport = totalSupport / len(treelist)
    return totalSupport

def unnormalized_transfer_support(bipartitions, treelist):
    """Unnormalized version of transfer  expectation.
    
    Due to the lack of normalization, contrary to the usual transfer_support, the lower value (close to zero) indicates well supported edges.

    Parameters
    ----------
    bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        Bipartitions to evaluate unnormalized transfer_support.
    treelist : TreeList_with_support or Iterable of either Tree_with_support or 
        Input trees to evaluate unnormalized transfer_support against.

    Returns
    -------
    numpy.ndarray
        Array of unnormalized transfer_support. 
    """
    n_taxa = len(treelist[0].taxon_namespace)
    refinfo = _create_refinfo(bipartitions, n_taxa)
    n_bipartitions = len(bipartitions)
    totalDistance = np.zeros(n_bipartitions)

    for tree in treelist:
        if tree.bipartition_encoding is None:
            tree.encode_bipartitions()
        tree_bipartition_ints = [bipartition.split_as_int() for bipartition in tree.bipartition_encoding]
        for i, bipartition in enumerate(bipartitions):
            bipar_int = bipartition.split_as_int()
            if bipar_int in tree_bipartition_ints:
                distance = 0
            elif refinfo[bipar_int][1] == 2:
                distance  = 1
            else:
                distance = _minDist(refinfo[bipar_int], tree)
            totalDistance[i] += distance
    
    totalDistance = totalDistance/ len(treelist)
    return totalDistance
