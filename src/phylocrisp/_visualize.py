from ete4 import Tree
from ete4.smartview import TextFace
import dendropy
from bitstring import Bits
import numpy as np
from ._consensus import Tree_with_support


def build_tree_with_support_faces(
    newick: str,
    namespace,
    branch_support: dict | None = None,
    transfer_support: dict | None = None,
    leaf_support: bool = False,
    attach_faces: bool = True,
):
    """
    Returns (t, support_maps)

    support_maps:
      { "branch": {node: value}, "transfer": {node: value} }
    """
    t = Tree(newick)

    support_maps = {"branch": {}, "transfer": {}}

    if branch_support is not None:
        _ = get_support(
            t, namespace, branch_support,
            pos=0, leaf_support=leaf_support,
            attach_faces=attach_faces,
            out_map=support_maps["branch"],
        )
    if transfer_support is not None:
        _ = get_support(
            t, namespace, transfer_support,
            pos=1, leaf_support=leaf_support,
            attach_faces=attach_faces,
            out_map=support_maps["transfer"],
        )
    return t, support_maps

def try_make_treestyle_with_legend():
    """
    GUI-only. If not available, returns None.
    """
    try:
        from ete4.treeview import TreeStyle
        ts = TreeStyle()
        return ts
    except Exception:
        return None

def plot_example_func(tree: Tree_with_support):
    # If GUI is missing, avoid calling add_face().
    ts = try_make_treestyle_with_legend()
    attach_faces = (ts is not None)

    t, support_maps = build_tree_with_support_faces(
        newick=tree.as_string(schema="newick", suppress_rooting=True),
        namespace=tree.taxon_namespace,
        branch_support=tree.branch_support,
        transfer_support=tree.transfer_support,
        leaf_support=False,
        attach_faces=attach_faces,
    )

    # --- second part (legend, etc.) rewritten ---
    if ts is not None:
        color = ["#006BA4", "#FF800E", "#ABABAB", "#595959",
                 "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]

        # NOTE: In ete4, adding faces/legend requires treeview(Qt).
        # We only do this when ts is available.
        if tree.branch_support is not None:
            f = TextFace("branch_support")
            setattr(f, "color", color[0])
            ts.legend.add_face(f, column=0)

        if tree.transfer_support is not None:
            f = TextFace("transfer_support")
            setattr(f, "color", color[1])
            ts.legend.add_face(f, column=0)

    # If you want to use values programmatically in headless mode,
    # support_maps contains {node: value} mappings.
    return t, ts, support_maps
    

def get_support(Node, namespace, support_hashtable, pos=0, leaf_support=True, attach_faces: bool = True, out_map: dict | None = None):
    """A function to add supports(branch support, transfer support) from The Tree of Tree_with_support　class to Nodes in ete4 for visualization.
    
        Parameters
        ---------
        Node: ete4.TreeNode
        namespace: A list of species of `dendropy`
        support_hashtable: dict(), key: clade_bit value: support_score
        pos:int, Position information for adding support to a Node using ete4.TreeNode.add_face
        leaf_support: bool
    
        Returns
        -------
        clade_bit: A bit string representing the bipartition of leaf species.
    """
    color = ["#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]
    taxonnames_array = np.array([item.label for item in namespace])
    clade_bool = [False for i in range(len(taxonnames_array))]
    if(Node.is_leaf):
        digits = [np.where( Node.name == taxonnames_array )][0][0][0]
        clade_bool[-(digits+1)] = True
        clade_bit = Bits(clade_bool)
        if int(clade_bit.bin[-1]) == 1:
            clade_bit = (~clade_bit)
    else:
        clade_bit = Bits(clade_bool)
        if int(clade_bit.bin[-1]) == 1:
            clade_bit = (~clade_bit)
        for child in Node.children:
            clade_bit = clade_bit | get_support(
                child,
                namespace,
                support_hashtable,
                pos=pos,
                leaf_support=leaf_support,
                attach_faces=attach_faces,
                out_map=out_map,
            )
    if Node.is_leaf == False or leaf_support == True:
        if clade_bit.uint in support_hashtable.keys():
            val = support_hashtable[clade_bit.uint]

            # store externally (ete4 nodes don't allow arbitrary attributes)
            if out_map is not None:
                out_map[Node] = val

            if attach_faces:
                color = ["#006BA4", "#FF800E", "#ABABAB", "#595959",
                         "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]
                textface = TextFace("{:.3f}".format(val))
                setattr(textface, "color", color[pos])

                if (pos % 2 == 0):
                    Node.add_face(textface, pos // 2, position="branch-top")
                else:
                    Node.add_face(textface, pos // 2, position="branch-bottom")

    return clade_bit

def read_consensus_NeXML(NeXML_path):
    newtree = dendropy.Tree.get(path=NeXML_path,schema="nexml")
    consensus = Tree_with_support(newtree)
    consensus.branch_support = get_support_from_NeXML(newtree,"branch_support")
    consensus.transfer_support = get_support_from_NeXML(newtree,"transfer_support")
    return consensus


def write_consensus_NeXML(Tree_with_support,NexML_path):
    if(Tree_with_support.branch_support != None):
        for edge in Tree_with_support.postorder_edge_iter():
            edge.branch_support = None
            edge.annotations.add_bound_attribute("branch_support")
            edge.branch_support = Tree_with_support.branch_support[int(edge.bipartition)]

    if(Tree_with_support.transfer_support != None):
        for edge in Tree_with_support.postorder_edge_iter():
            edge.transfer_support = None
            edge.annotations.add_bound_attribute("transfer_support")
            edge.transfer_support = Tree_with_support.transfer_support[int(edge.bipartition)]

    Tree_with_support.write(
        path= NexML_path,
        schema='nexml',
        ignore_unrecognized_keyword_arguments=False,
        )

def get_support_from_NeXML(dendropy_Tree_from_NeXML,support_name):
    if(support_name != "branch_support" and support_name != "transfer_support"):
        raise ValueError(
            "support_name is invalid; select either 'branch_support' or 'transfer_support'."
        )

    support_added = dict()
    dendropy_Tree_from_NeXML.encode_bipartitions()
    if(support_name == "branch_support"): 
        for edge in dendropy_Tree_from_NeXML.postorder_edge_iter():
            if edge.annotations.find(name="branch_support") == None: continue
            support_added[int(edge.bipartition)] = float(str(edge.annotations.find(name="branch_support")).strip("branch_support=").strip("'"))
    if(support_name == "transfer_support"):
        for edge in dendropy_Tree_from_NeXML.postorder_edge_iter():
            if edge.annotations.find(name="transfer_support") == None: continue
            support_added[int(edge.bipartition)] = float(str(edge.annotations.find(name="transfer_support")).strip("transfer_support=").strip("'"))

    return support_added
    
    
