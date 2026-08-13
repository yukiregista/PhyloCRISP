import ctypes
import importlib.util
import os
import platform
import sys
from pathlib import Path

from bitstring import Bits
import numpy as np

class Split(ctypes.Structure):
    _fields_ = [("data", ctypes.POINTER(ctypes.c_uint))]

class BipartitionMatches(ctypes.Structure):
    _fields_ = [
        ("bipartition", ctypes.POINTER(ctypes.c_uint)),
        ("topo_depth", ctypes.c_int),
        ("matches", ctypes.POINTER(ctypes.POINTER(ctypes.c_uint))),
        ("td", ctypes.POINTER(ctypes.c_int)),
        ("is_external", ctypes.c_bool)
    ]

class BipartitionDict(ctypes.Structure):
    _fields_ = [
        ("entries", ctypes.POINTER(ctypes.POINTER(BipartitionMatches))),
        ("num_entries", ctypes.c_int),
        ("bipartition_size", ctypes.c_size_t), # number of unsigned integers used to represent the bipartition
        ("num_matches", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
        ("uint_bits", ctypes.c_int) # number of bits in one unsigned integer
    ]


class BipartitionDictArray(ctypes.Structure):
    _fields_ = [
        ("bdict", ctypes.POINTER(BipartitionDict)),
        ("num_trees", ctypes.c_int)
    ]



class BipartitionSupport(ctypes.Structure):
    _fields_ = [
        ("bipartition", ctypes.POINTER(ctypes.c_uint)),
        ("is_external", ctypes.c_bool),
        ("support", ctypes.c_double)
    ]
    
class BipartitionSupportArray(ctypes.Structure):
    _fields_ = [
        ("num_entries", ctypes.c_int),
        ("bipartition_size", ctypes.c_size_t),
        ("n_taxa", ctypes.c_int),
        ("uint_bits", ctypes.c_int),
        ("bipartition_supports", ctypes.POINTER(ctypes.POINTER(BipartitionSupport)))
    ]


class BipartitionMatchesHash(ctypes.Structure):
    _fields_ = [
        ("h1", ctypes.c_uint),
        ("h2", ctypes.c_uint),
        ("topo_depth", ctypes.c_int),
        ("matched_ids", ctypes.POINTER(ctypes.c_int)),
        ("td", ctypes.POINTER(ctypes.c_int)),
        ("input_count", ctypes.c_int)
    ]

class BipartitionDictHash(ctypes.Structure):
    _fields_ = [
        ("entries", ctypes.POINTER(ctypes.POINTER(BipartitionMatchesHash))),
        ("num_entries", ctypes.c_int),
        ("num_matches", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
    ]
    

class BipartitionSupportHash(ctypes.Structure):
    _fields_ = [
        ("h1", ctypes.c_uint),
        ("h2", ctypes.c_uint),
        ("is_external", ctypes.c_bool),
        ("support", ctypes.c_double),
        ("node_id", ctypes.c_int),
        ("topo_depth", ctypes.c_int)
        ]
    
class BipartitionSupportArrayHash(ctypes.Structure):
    _fields_ = [
        ("num_entries", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
        ("bipartition_supports", ctypes.POINTER(ctypes.POINTER(BipartitionSupportHash))),
        ("num_alt_trees", ctypes.c_int)
    ]
    
    
class TbeSupportMatchResult(ctypes.Structure):
    _fields_ = [
        ("bdict", ctypes.POINTER(BipartitionDictHash)),
        ("bsupp_arr", ctypes.POINTER(BipartitionSupportArrayHash))
    ]
    
    
class transferDistance(ctypes.Structure):
    _fields_ = [
        ("dist", ctypes.c_int),
        ("h1", ctypes.c_uint),
        ("h2", ctypes.c_uint),
        ("node_id", ctypes.c_int)
    ]

class transferDistanceAll(ctypes.Structure):
    _fields_ = [
        ("td", ctypes.POINTER(ctypes.POINTER(transferDistance))),
        ("n_elements", ctypes.c_int),
        ("node_h1", ctypes.c_uint),
        ("node_h2", ctypes.c_uint),
        ("topo_depth", ctypes.c_int)
    ]


class OrganizedTransferResult():
    def __init__(self, TS: list, REF_BIPAR_ID_to_location: dict, TD: list, H1: list, H2: list, S: np.ndarray, M: np.ndarray, C: list, num_alt_trees: int):
        self.TS = TS # Tranfer Supports
        # Dictionary containing bipar_id of reference tree's internal bipartition as keys, and its location in the TS array as values
        self.REF_BIPAR_ID_to_location = REF_BIPAR_ID_to_location
        self.TD = TD # Topological Depth.
        self.H1 = H1 # stores hash h1 of input bipartitions
        self.H2 = H2 # stores hash h2 of input bipartitions
        self.S = S # stores transfer distance (either scaled or unscaled) to the first K best matches: size m \times K.
        self.M = M # stores first K best matched IDs : size m \times K.
        self.C = C # stores the number of appearances for each bipartition in the input trees
        self.num_alt_trees = num_alt_trees # Number of alt trees
        

_DLL_DIR_HANDLES = []


def _register_windows_dll_dir(path: str):
    """Register a DLL search directory on Windows (Python 3.8+)."""
    if os.name != "nt" or not path:
        return
    add_dll_directory = getattr(os, "add_dll_directory", None)
    if add_dll_directory is None:
        return
    try:
        handle = add_dll_directory(path)
    except OSError:
        return
    _DLL_DIR_HANDLES.append(handle)


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


def _iter_booster_candidate_paths():
    """Yield likely paths for a ctypes-loadable booster library."""
    seen = set()

    def _yield(path_like):
        if not path_like:
            return
        path = os.path.abspath(os.fspath(path_like))
        norm = os.path.normcase(path)
        if norm in seen:
            return
        seen.add(norm)
        if os.path.exists(path):
            yield path

    env_path = os.environ.get("CONSENSUS_BOOSTER_DLL")
    if env_path:
        yield from _yield(env_path)

    module_dir = Path(__file__).resolve().parent
    platform_dir = Path(_platform_folder())
    arch_dir = Path(_arch_folder())
    bin_dir = module_dir / "bin"
    parent_bin_dir = module_dir.parent / "bin"
    candidate_dirs = [
        bin_dir,
        bin_dir / platform_dir / arch_dir,
        bin_dir / platform_dir / "universal2",
        module_dir,
        parent_bin_dir,
        parent_bin_dir / platform_dir / arch_dir,
        parent_bin_dir / platform_dir / "universal2",
        module_dir.parent,
    ]

    if hasattr(sys, "_MEIPASS"):
        meipass = Path(sys._MEIPASS)
        candidate_dirs.insert(0, meipass)
        candidate_dirs.insert(1, meipass / "_internal")

    exe_dir = Path(sys.executable).resolve().parent if sys.executable else None
    if exe_dir:
        candidate_dirs.append(exe_dir)
        candidate_dirs.append(exe_dir / "_internal")

    for base in candidate_dirs:
        for name in ("booster.dll", "booster.pyd", "libbooster.so", "libbooster.dylib"):
            yield from _yield(base / name)


def _load_ctypes_library(path: str) -> ctypes.CDLL:
    _register_windows_dll_dir(os.path.dirname(path))
    return ctypes.CDLL(path)


def load_booster() -> ctypes.CDLL:
    load_errors = []

    for candidate in _iter_booster_candidate_paths():
        try:
            return _load_ctypes_library(candidate)
        except OSError as e:
            load_errors.append(f"{candidate}: {e}")

    spec = importlib.util.find_spec("booster")

    if spec is not None and spec.origin is not None:
        try:
            return _load_ctypes_library(spec.origin)
        except OSError as e:
            load_errors.append(f"{spec.origin}: {e}")

    message = "Cannot find/load the 'booster' shared library"
    if load_errors:
        message += "\nTried:\n- " + "\n- ".join(load_errors)
    raise ImportError(message)

def int_to_reversed_bin(x, bit_length):
    '''
    This conversion is inefficient!: 
    It would be better if we also use the similar structure (array of uints)
    in Python as well.
    '''
    return ''.join('1' if x & (1 << i) else '0' for i in range(bit_length))


def to_bitvector(bipartition: ctypes.POINTER(ctypes.c_uint), n_tips: int=100, bit_length: int=32) -> Bits:
    n_uints = (n_tips-1) // bit_length + 1
    bipartition_array = [bipartition[i] for i in range(n_uints)]
    bit_size = [bit_length for _ in range(n_uints-1)]
    bit_size.append(n_tips - bit_length * (n_uints-1))
    bit_length = 32
    binary_strings = [int_to_reversed_bin(x, y) for x, y in zip(bipartition_array, bit_size)]
    concatenated_bits = ''.join(binary_strings)
    return Bits(bin=concatenated_bits[::-1])


def prepare_tbe_match(booster_lib: ctypes.CDLL) -> ctypes.CDLL:
    """Prepare `tbe_match` and `after_tbe_match`

    Args:
        booster_lib (ctypes.CDLL): booster library.
    """
    booster_lib.tbe_match.restype = ctypes.c_void_p # To get pointer accurately
    booster_lib.after_tbe_match.argtypes = ctypes.c_void_p,
    return booster_lib

def prepare_tbe_support(booster_lib: ctypes.CDLL) -> ctypes.CDLL:
    """Prepare `tbe_match` and `after_tbe_match`

    Args:
        booster_lib (ctypes.CDLL): booster library.
    """
    booster_lib.tbe_support.restype = ctypes.c_void_p # To get pointer accurately
    booster_lib.after_tbe_support.argtypes = ctypes.c_void_p,
    return booster_lib


def prepare_tbe_support_and_match(booster_lib: ctypes.CDLL) -> ctypes.CDLL:
    booster_lib.tbe_support_and_match.restype = ctypes.POINTER(TbeSupportMatchResult)

    # Specify the argument types of the function
    booster_lib.tbe_support_and_match.argtypes = [
        ctypes.c_char_p,                    # char* nh_initial_tree
        ctypes.POINTER(ctypes.c_char_p),    # char** nh_input_trees
        ctypes.c_int,                       # int num_trees
        ctypes.c_int                        # int k
    ]
    return booster_lib


def prepare_recompute(booster_lib: ctypes.CDLL):
    booster_lib.recompute.restype = ctypes.POINTER(transferDistanceAll)

    # Specify the argument types of the function
    booster_lib.recompute.argtypes = [
        ctypes.c_uint,                       # unsigned int h1
        ctypes.c_uint                        # unsigned int h2
    ]
    return booster_lib

def prepare_prune_and_return_newick(booster_lib: ctypes.CDLL):
    booster_lib.prune_and_return_newick.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int
    ]
    booster_lib.prune_and_return_newick.restype = ctypes.POINTER(ctypes.c_char)
    
    booster_lib.free_buffer.argtypes = [ctypes.c_char_p]
    booster_lib.free_buffer.restype = None
    return booster_lib



def prepare_tbe_match_args(reftrees_str:str, alttree_str: str, K: int):
    tbe_match_args = ['program_name', '-k', f'{K}']
    argc = len(tbe_match_args)
    argv_ctypes = (ctypes.POINTER(ctypes.c_char) * argc)()
    for i, str in enumerate(tbe_match_args):
        enc_str = str.encode('utf-8')
        argv_ctypes[i] = ctypes.create_string_buffer(enc_str)  
        
    tree1_lines = reftrees_str.split("\n")
    tree1_lines = [item for item in tree1_lines if len(item)>0]
    num_trees = len(tree1_lines)
    tree1_ctypes = (ctypes.POINTER(ctypes.c_char) * num_trees)()
    tree2_ctypes = ctypes.create_string_buffer(alttree_str.encode('utf-8'))
    for i in range(num_trees):
        buffer = ctypes.create_string_buffer(tree1_lines[i].encode('utf-8'))
        tree1_ctypes[i] = buffer
    
    num_trees_ctypes = ctypes.c_int(num_trees)    
    return argc, argv_ctypes, tree1_ctypes, tree2_ctypes, num_trees_ctypes


def prepare_tbe_support_args(reftree_str:str, alttrees_str: str):
    tree2_lines = alttrees_str.split("\n")
    tree2_lines = [item for item in tree2_lines if len(item)>0]
    num_trees = len(tree2_lines)
    tree2_ctypes = (ctypes.POINTER(ctypes.c_char) * num_trees)()
    tree1_ctypes = ctypes.create_string_buffer(reftree_str.encode('utf-8'))
    for i in range(num_trees):
        buffer = ctypes.create_string_buffer(tree2_lines[i].encode('utf-8'))
        tree2_ctypes[i] = buffer
    
    num_trees_ctypes = ctypes.c_int(num_trees)    
    return tree1_ctypes, tree2_ctypes, num_trees_ctypes



def prepare_tbe_support_and_match_args(reftree_str: str, input_trees_str: str, K: int):
    reference_tree_cstr = ctypes.c_char_p(reftree_str.encode('utf-8'))
    input_trees_lines = input_trees_str.split("\n")
    input_trees_lines = [item for item in input_trees_lines if len(item)>0] # ignore empty line.
    input_trees_array = (ctypes.c_char_p * len(input_trees_lines))(*[line.encode('utf-8') for line in input_trees_lines])
    num_trees = len(input_trees_lines)
    return reference_tree_cstr, input_trees_array, num_trees, K


def prepare_prune_and_return_newick_args(node_ids: list[int]):
    n_nodes = len(node_ids)
    c_array = (ctypes.c_int * n_nodes)(*node_ids)
    return c_array, n_nodes


def _generate_bitmask(bits: int)-> int:
    return (1<<bits) - 1

def organize_support_and_match(res_ptr, scaled=True) -> OrganizedTransferResult:
    bdict = res_ptr.contents.bdict.contents
    bsupp_arr = res_ptr.contents.bsupp_arr.contents
    num_trees = bsupp_arr.num_alt_trees
    
    # Create TS, H1_REF, H2_REF array. 
    TS = []; REF_BIPAR_IDS=[]; REF_BIPAR_ID_to_location = dict()
    TD = []
    count = 0
    for i in range(bsupp_arr.num_entries):
        entry = bsupp_arr.bipartition_supports[i].contents
        if not entry.is_external:
            TS.append(entry.support)
            TD.append(entry.topo_depth)
            REF_BIPAR_IDS.append(entry.node_id)
            REF_BIPAR_ID_to_location[entry.node_id] = count
            count += 1
    
    
    # Create Match, H1, H2, Score, Count
    K = bdict.num_matches; m = bdict.num_entries
    H1=[]; H2 =[]; S = np.zeros((m, K)); M = np.zeros((m,K), dtype=int); C = []
    count = 0
    for i in range(m):
        entry = bdict.entries[i].contents
        H1.append(entry.h1)
        H2.append(entry.h2)
        C.append(entry.input_count)
        
        # Use -1 if it reaches the maximum td.
        use_external=False
        if scaled:
            for j in range(K): 
                if use_external:
                    S[count][j] = 1 # maximum normalized td
                    M[count][j] = -1 # match to external
                elif (entry.td[j]>=entry.topo_depth-1):
                    use_external=True
                    S[count][j] = 1 # maximum normalized td
                    M[count][j] = -1 # match to external
                else:
                    S[count][j] = entry.td[j]/(entry.topo_depth-1)
                    M[count][j] = REF_BIPAR_ID_to_location[entry.matched_ids[j]]
        else:
            for j in range(K): 
                if use_external:
                    S[count][j] = entry.topo_depth-1 # maximum normalized td
                    M[count][j] = -1 # match to external
                elif (entry.td[j]>=entry.topo_depth-1):
                    use_external=True
                    S[count][j] = entry.topo_depth-1 # maximum normalized td
                    M[count][j] = -1 # match to external
                else:
                    S[count][j] = entry.td[j]
                    M[count][j] = REF_BIPAR_ID_to_location[entry.matched_ids[j]]
        count += 1
    
    
    
    return OrganizedTransferResult(TS, REF_BIPAR_ID_to_location, TD, H1, H2, S, M, C, num_trees)
            
            
        
    


def match_score_list(res_ptr: int, id_dict_ref: dict, K: int):
    bdict_arr: BipartitionDictArray = ctypes.cast(res_ptr, ctypes.POINTER(BipartitionDictArray)).contents  
    num_trees: int = bdict_arr.num_trees
    
    id = 0
    
    match_list = []
    score_list = []
    id_dict = dict()
    count_list = []
    topo_depth_list = []
    n_bits = bdict_arr.bdict[0].uint_bits
    mask = _generate_bitmask(n_bits)
    mask_last = _generate_bitmask(bdict_arr.bdict[0].n_taxa % n_bits)
    for i_tree in range(num_trees):
        bipartition_dict = bdict_arr.bdict[i_tree]
        for i in range(bipartition_dict.num_entries): 
            entry = bipartition_dict.entries[i]
            size = bipartition_dict.bipartition_size
            if(entry.contents.is_external): 
                continue
            
            bipartition = entry.contents.bipartition
            if bipartition[0] %2:
                bipartition = tuple(entry.contents.bipartition[item] for item in range(size))
            else:
                bipartition = tuple([~entry.contents.bipartition[item] & mask for item in range(size-1)] + [~entry.contents.bipartition[size-1] & mask_last])
            
            if bipartition in id_dict.keys():
                count_list[id_dict[bipartition]] += 1
                continue
            
            id_dict[bipartition] = id
            topo_depth = entry.contents.topo_depth
            better_than_external = True
            matches = []; scores =[]
            for i in range(K):
                if not better_than_external:
                    matches.append(-1)
                    scores.append(topo_depth-1)
                elif topo_depth-2 < entry.contents.td[i]:
                    matches.append(-1)
                    scores.append(topo_depth-1)
                    better_than_external = False
                else:
                    match_bipar = entry.contents.matches[i]
                    if match_bipar[0]%2:
                        match_bipar = id_dict_ref[tuple(match_bipar[i] for i in range(size))]
                    else:
                        complemented_key = tuple(
                            [~match_bipar[i] & mask for i in range(size - 1)]
                            + [~match_bipar[size - 1] & mask_last]
                        )
                        try:
                            match_bipar = id_dict_ref[complemented_key]
                        except KeyError as e:
                            raise RuntimeError(
                                "Failed to map complemented matched bipartition to reference ids in "
                                "match_score_list(). This usually indicates a mismatch between "
                                "booster output and reference bipartition encoding."
                            ) from e
                    matches.append(match_bipar)
                    scores.append((entry.contents.td[i]))
            match_list.append(matches)
            score_list.append(scores)
            count_list.append(1)
            topo_depth_list.append(entry.contents.topo_depth)
            id += 1
    
    return id_dict, count_list, match_list, score_list, topo_depth_list, num_trees

def tbe_support(res_ptr: int):
    bdict_supp: BipartitionSupportArray = ctypes.cast(res_ptr, ctypes.POINTER(BipartitionSupportArray)).contents
    num_edges = bdict_supp.num_entries
    size = bdict_supp.bipartition_size
    support_list = []; id_dict = dict()
    mask = _generate_bitmask(bdict_supp.uint_bits)
    mask_last = _generate_bitmask(bdict_supp.n_taxa % bdict_supp.uint_bits)
    id = 0
    for i in range(num_edges):
        if bdict_supp.bipartition_supports[i].contents.is_external:
            continue
        bipartition = bdict_supp.bipartition_supports[i].contents.bipartition
        if bipartition[0] %2:
            bipartition = tuple(bipartition[item] for item in range(size))
        else:
            bipartition = tuple([~bipartition[item] & mask for item in range(size-1)] + [~bipartition[size-1] & mask_last])
        support_list.append(bdict_supp.bipartition_supports[i].contents.support)
        id_dict[bipartition] = id
        id += 1
    return id_dict, support_list, bdict_supp.uint_bits, bdict_supp.n_taxa
        
    
    
    
