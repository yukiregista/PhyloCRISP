import dendropy
import random
import numpy as np

import sys
sys.setrecursionlimit(10000)


def birthdeath_sampling(n_taxa=100, birth_rate = 0.1, death_rate = 0.1, gsa_prop = 10, seed = None) -> dendropy.Tree:
    """Implements the general sampling approach using birth-death model and returns unrooted tree.

    Parameters
    ----------
    n_taxa : int, optional
        number of tips (taxa), by default 100
    birth_rate : float, optional
        birth rate, by default 0.1
    death_rate : float, optional
        death rate, by default 0.1
    gsa_prop : int, optional
        Birth-death model will be simulated until `gsa_prop` * `n_taxa` tips are obtained, by default 10
    seed : int, optional
        random seed, by default None

    Returns
    -------
    dendropy.Tree
        Sampled tree.
    """
    
    # create birth-death tree
    t = dendropy.simulate.treesim.birth_death_tree(birth_rate = birth_rate, death_rate = death_rate, gsa_ntax = int(n_taxa * gsa_prop), 
                                                    num_extant_tips=n_taxa, is_retain_extinct_tips=True,is_assign_extinct_taxa = False, 
                                                    rng = random.Random(seed))
    
    # only retain extant taxa
    ll = []
    count2 = 0
    for node in t.nodes():
        if node.is_extinct is False:
            count2 += 1
            node.taxon.label = f"S{count2}"
            ll.append(node.taxon.label)
    t = t.extract_tree_with_taxa_labels(ll)
    t.purge_taxon_namespace()
    t.seed_node.edge_length=None
    
    return t

def normalize_tree_with_lognormal(tree, height = 0.05, lognormal_mean=1, lognormal_std=0.5) -> dendropy.Tree:
    """Normalize edge length to have the specified height, then each edge length is multiplied by lognormal variable.
    If `tree.is_rooted` is `True`, the height is defined to be the (maximum) distance from the root to the tips.
    If `tree.is_rooted` is `False`, then the height is defined to be (the maximum patristic distance)/2. 

    Parameters
    ----------
    tree : `dendropy.Tree` or `Tree_with_support`
        The tree to apply normalization.
    height : float, optional
        see above, by default 0.05
    lognormal_mean : int, optional
        The mean of lognormal variable to multiply for each edge, by default 1
    lognormal_std : float, optional
        The standard deviation of lognormal variable to multiply for each edge, by default 0.5

    Returns
    -------
    dendropy.Tree
        Normalized and perturbed tree.
    """
    
    # compute lognormal parameter
    lognormal_var = lognormal_std**2
    scale2 = np.log(lognormal_var / lognormal_mean + 1)
    scale = np.sqrt(scale2)
    mu = np.log(lognormal_mean) - scale2 / 2
    
    tt =tree.clone(depth=1)

    # adjust length
    maxdist = tt.max_distance_from_root()
    for edge in tt.postorder_edge_iter():
        if edge.length is not None:
            edge.length = edge.length/(maxdist) * height
            edge.length = edge.length * np.random.lognormal(mean=mu, sigma=scale)
    
    return tt
    
