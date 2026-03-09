'''
- Bipartitionは  0,1のbitvectorとなっています。
- Cではすべてrooted treeとして扱っており、枝のsubtreeに含まれるtaxaに対応するbitが1になっています。
- taxaはアルファベット順に並び替えられています。bitvectorのk番目の要素はアルファベット順k番目のtaxonに対応しています。
'''


import ctypes

from bitstring import Bits


# Number of matches
K = 10

# Load the shared library
lib = ctypes.CDLL('../src/libbooster.so')

class Split(ctypes.Structure):
    _fields_ = [("data", ctypes.POINTER(ctypes.c_uint))]

class BipartitionMatches(ctypes.Structure):
    _fields_ = [
        ("bipartition", ctypes.POINTER(ctypes.c_uint)),
        ("matches", ctypes.POINTER(ctypes.POINTER(ctypes.c_uint))),
        ("td", ctypes.POINTER(ctypes.c_int)),
        ("is_external", ctypes.c_bool)
    ]

class BipartitionDict(ctypes.Structure):
    _fields_ = [
        ("entries", ctypes.POINTER(ctypes.POINTER(BipartitionMatches))),
        ("num_entries", ctypes.c_int),
        ("bipartition_size", ctypes.c_size_t),
        ("num_matches", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
        ("bits_uint", ctypes.c_int)
    ]

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
    return Bits(bin=concatenated_bits)


# Define the argument and return types of the C function
# lib.tbe_match.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.c_char_p, ctypes.POINTER(ctypes.c_char_p)]
#lib.tbe_match.restype = ctypes.c_int
lib.tbe_match.restype = ctypes.c_void_p # To get pointer accurately
lib.print_bdict.argtypes = ctypes.c_void_p,
lib.after_tbe_match.argtypes = ctypes.c_void_p,

# # Create a list of strings
custom_args = ['program_name', '-k', f'{K}']

# # Convert list of strings to ctypes array of char*
argc = len(custom_args)
#argv_ctypes = (ctypes.c_char_p * argc)(*map(ctypes.create_string_buffer, custom_args))

argv_ctypes = (ctypes.POINTER(ctypes.c_char) * argc)()
for i, str in enumerate(custom_args):
    enc_str = str.encode('utf-8')
    argv_ctypes[i] = ctypes.create_string_buffer(enc_str)

# Call the C function
## newick strings
nh_tree1 = "((S1,(((((((((((S69,(S71,S70)0.950267:0.050038)1.000000:0.561803,S72)1.000000:6.503290,S73)1.000000:1.071165,(S74,S75)1.000000:5.456159)1.000000:4.485022,(S76,S77)1.000000:5.116995)1.000000:4.629517,((((((S60,S59)1.000000:6.503290,S58)1.000000:6.503290,(S56,S57)1.000000:6.503290)1.000000:0.404785,S61)1.000000:0.506943,(((S62,S63)0.727403:0.047200,S64)1.000000:1.816539,S65)1.000000:6.503290)1.000000:0.216140,((S68,S67)1.000000:6.503290,S66)1.000000:6.503290)1.000000:1.160477)1.000000:2.540835,(((((S47,(S45,S46)1.000000:2.892372)1.000000:3.665186,S44)1.000000:2.244026,S55)1.000000:0.266510,(((S49,S50)1.000000:0.166464,S48)1.000000:6.503290,((S51,S52)1.000000:6.503290,S53)1.000000:0.149549)1.000000:0.135687)1.000000:0.297050,S54)1.000000:6.503290)1.000000:0.931136,(((((S94,S95)1.000000:5.411872,((((S97,S96)1.000000:0.243235,S99)0.991032:0.064405,S100)1.000000:0.295985,S98)1.000000:0.714942)1.000000:0.263404,S93)1.000000:6.503290,((((S85,(S86,S87)1.000000:6.503290)1.000000:6.503290,S88)1.000000:0.591315,(S89,(S92,(S90,S91)1.000000:6.503290)1.000000:2.072473)1.000000:6.503290)0.998856:0.079637,(S83,S84)1.000000:6.503290)1.000000:6.503290)1.000000:0.288682,((S78,(S80,S79)1.000000:2.866246)1.000000:2.212830,(S81,S82)1.000000:6.503290)1.000000:4.423848)1.000000:6.503290)1.000000:6.503290,((((S41,S42)1.000000:6.503290,S43)1.000000:4.423848,(S39,S40)1.000000:5.404677)1.000000:6.503290,(S38,S37)1.000000:6.503290)1.000000:3.938340)1.000000:6.503290,((((S15,(S13,S14)1.000000:6.503290)1.000000:6.503290,(S11,S12)1.000000:6.503290)1.000000:5.499988,((S9,S10)1.000000:1.022651,(((S6,S5)1.000000:2.742569,S4)1.000000:3.795239,(S8,S7)1.000000:6.503290)1.000000:6.503290)1.000000:6.503290)1.000000:1.952685,S3)1.000000:0.377830)1.000000:0.762649,(S36,(((S23,S22)1.000000:3.367795,(((S21,S16)0.883454:0.043581,(S17,S18)1.000000:5.418664)0.999788:0.097235,(S19,S20)1.000000:3.135994)1.000000:6.503290)1.000000:4.128384,((S30,(S24,((S29,S28)1.000000:6.503290,((S25,S26)1.000000:6.503290,S27)1.000000:1.272463)1.000000:0.533730)1.000000:6.166817)1.000000:6.503290,(S35,(S31,((S33,S32)1.000000:6.503290,S34)1.000000:2.789718)1.000000:6.503290)1.000000:6.503290)1.000000:0.383954)1.000000:3.466735)1.000000:0.355515)1.000000:5.898335),S2);"
nh_tree2 = "((S1,(((((((((((S61,S60)1.000000:6.503290,(S59,S58)1.000000:0.155718)1.000000:4.151914,((S65,((S64,S63)1.000000:0.356960,S62)1.000000:6.097825)1.000000:3.458767,S66)1.000000:6.503290)1.000000:6.503290,S67)1.000000:6.503290,((((((S51,S52)1.000000:6.503290,S50)1.000000:0.448899,S53)1.000000:0.662403,S54)1.000000:5.133434,((S55,S56)1.000000:6.503290,S57)1.000000:2.479881)1.000000:0.522718,((S48,S49)1.000000:4.711530,S47)1.000000:6.503290)1.000000:6.503290)1.000000:4.893852,((((S42,S38)0.901605:0.122038,((S40,S41)1.000000:4.306065,S39)1.000000:1.811942)1.000000:6.503290,(S44,S43)1.000000:0.518244)1.000000:1.261543,(S45,S46)1.000000:6.503290)1.000000:6.503290)1.000000:0.793446,((((((((S73,(S75,S74)0.999999:0.197014)1.000000:6.503290,S76)1.000000:3.006782,S77)1.000000:5.810142,((S79,S80)1.000000:6.503290,S78)1.000000:1.506077)1.000000:3.730701,(S82,S81)1.000000:6.503290)1.000000:6.503290,((((S69,S70)1.000000:6.503290,S71)1.000000:2.122436,S72)1.000000:0.444946,S68)1.000000:6.503290)1.000000:3.207453,((S84,(S86,S85)1.000000:2.048942)1.000000:6.503290,S83)1.000000:6.503290)1.000000:5.810142,((S100,(S96,((S99,S98)1.000000:6.503290,S97)1.000000:6.503290)1.000000:2.495956)1.000000:6.503290,(((S88,S87)1.000000:2.719100,((S89,S90)1.000000:6.503290,(S93,(S92,S91)1.000000:6.482881)1.000000:2.839728)1.000000:4.711530)1.000000:6.503290,(S94,S95)1.000000:0.942608)1.000000:6.503290)1.000000:6.503290)1.000000:2.838205)1.000000:0.926660,((((((((S22,S23)1.000000:6.503290,((S25,S24)1.000000:6.503290,((S27,S26)1.000000:6.503290,(S29,S28)1.000000:6.503290)1.000000:6.503290)1.000000:0.882456)1.000000:0.243390,(S31,S30)1.000000:1.035230)1.000000:6.503290,((S34,S35)1.000000:6.503290,(S33,S32)1.000000:4.694318)1.000000:1.878317)1.000000:6.503290,(S36,S37)1.000000:6.503290)1.000000:0.352739,(S20,S21)1.000000:6.503290)1.000000:3.265949,((S17,((S15,S16)1.000000:6.503290,S18)0.999963:0.117305)1.000000:0.814057,S19)1.000000:6.503290)1.000000:2.108841,((S14,(S12,S13)1.000000:1.928579)1.000000:6.503290,(((S8,S9)1.000000:1.414358,S11)1.000000:0.490798,S10)1.000000:6.503290)1.000000:6.503290)1.000000:5.810142)1.000000:0.225176,S7)1.000000:6.503290,(S5,S6)1.000000:0.421071)1.000000:6.503290,(S3,S4)1.000000:4.423848)1.000000:6.503290),S2);"

nh_tree2_ctypes = (ctypes.POINTER(ctypes.c_char) * 1)()
nh_tree2_ctypes[0] = ctypes.create_string_buffer(nh_tree2.encode('utf-8'))

print(type(argv_ctypes), type(ctypes.create_string_buffer(nh_tree1.encode('utf-8'))), type(nh_tree2_ctypes))

# これはポインタ
res = lib.tbe_match(argc, argv_ctypes, ctypes.create_string_buffer(nh_tree1.encode('utf-8')) ,nh_tree2_ctypes)

# これはBipartitionDict
bipartition_dict = ctypes.cast(res, ctypes.POINTER(BipartitionDict)).contents
print(bipartition_dict) # 

for i in range(bipartition_dict.num_entries): 
    entry = bipartition_dict.entries[i] # 各枝に対応: BipartitionMatches インスタンスへのポインタ
    if(entry.contents.is_external): 
        continue # 外部枝にはmatchの情報がありません。
    matches = entry.contents.matches # matchのポインタ配列
    bipartition = to_bitvector(entry.contents.bipartition) # 枝のbipartitionのBitsインスタンス
    #print(bipartition.bin) # Bipartitionのプリント
    matches = [] # First K match bipartition をストアするリスト
    scores = [] # First K match の transfer distance をストアするリスト
    for i in range(K):
        matches.append(to_bitvector(entry.contents.matches[i]))
        scores.append(entry.contents.td[i])


lib.after_tbe_match(res) # メモリをFreeするための関数なので、tbe_match()をして中身を見終わったらこれを必ずしてください！
