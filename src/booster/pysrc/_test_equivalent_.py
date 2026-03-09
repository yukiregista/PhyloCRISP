import dendropy
import numpy as np


# from the edited code.
t1_data = "(S100:0,(S98:0,((S94:0,S93:0)1:6.50329,S97:0)1:2.86291,(S96:0,S95:0)1:6.50329,((((S69:0,((S72:0,S70:0)0.945896:0.056823,((S74:0,S75:0)1:2.31364,S73:0)1:4.10824,S71:0)1:0.435157,(S80:0,(S77:0,S76:0)1:6.50329,S79:0,S78:0)1:0.8078,(((S89:0,S90:0)1:6.50329,S88:0,S87:0)0.904872:0.047802,S86:0)1:3.0454,((S84:0,S85:0)1:0.203309,(S83:0,S81:0,S82:0)1:3.97756)1:3.58552)1:0.925048,S36:0,S41:0,(S46:0,S47:0)1:6.50329,S40:0,((S42:0,S44:0)1:0.14434,S43:0)0.919764:0.054413,(S39:0,S38:0)1:3.41656,S37:0,S48:0,((S54:0,S56:0,S57:0,(S59:0,S58:0)0.854019:0.039115,S60:0,S55:0)1:2.02583,((S49:0,S50:0)0.958314:0.051241,S51:0)1:6.50329,(S53:0,S52:0)1:4.14343)1:0.775676,(S62:0,S65:0,S66:0,S61:0,S63:0,S64:0,S68:0,S67:0)1:6.50329,S45:0)1:6.50329,(S91:0,S92:0)1:6.50329)1:1.33427,((((S3:0,(S1:0,S2:0)1:6.50329)1:3.55885,(S7:0,((S4:0,S6:0)0.422885:0.008834,S5:0)1:2.66782,((S9:0,S10:0)1:1.79837,(S11:0,S12:0)1:4.29972)1:1.85003,S8:0)1:6.50329)1:3.3678,((S15:0,(S14:0,S19:0)0.885403:0.044103,S13:0,S18:0,S17:0,S16:0)1:2.43944,(S34:0,S32:0,S31:0,S33:0,S30:0,S28:0,S29:0)1:6.50329,(S20:0,S21:0)1:1.6034,S22:0,(S25:0,S23:0,S26:0,S24:0)1:0.919886,S27:0)1:6.50329)1:2.2634,S35:0)1:0.272108)1:6.50329)1:2.60769,S99:0);"
# from the original booster code
t2_data = "(S1,S2,(S3,(((S35,(((S36,S37,S40,S41,S45,S48,(S43,(S42,S44)),(S38,S39),(S46,S47),((S54,S55,S56,S57,S60,(S58,S59)),(S51,(S49,S50)),(S52,S53)),(S61,S62,S63,S64,S65,S66,S67,S68),(S69,(S71,(S73,(S74,S75)),(S70,S72)),(S78,S79,S80,(S76,S77)),(S86,(S87,S88,(S89,S90))),((S84,S85),(S81,S82,S83)))),(S91,S92)),(S98,(S95,S96),(S100,S99),(S97,(S93,S94))))),(S22,S27,(S13,S15,S16,S17,S18,(S14,S19)),(S28,S29,S30,S31,S32,S33,S34),(S23,S24,S25,S26),(S20,S21))),(S7,S8,(S5,(S4,S6)),((S10,S9),(S11,S12))))));"

t1 = dendropy.Tree.get(schema="newick", data = t1_data)
t1.encode_bipartitions()
t2 = dendropy.Tree.get(schema="newick", data = t2_data, taxon_namespace=t1.taxon_namespace)
t2.encode_bipartitions()

print("topological difference:", dendropy.calculate.treecompare.false_positives_and_negatives(t1, t2))


t1_nodes = t1.internal_nodes(exclude_seed_node=True)
t2_nodes = t2.internal_nodes(exclude_seed_node=True)

t1_node_splits = [item.bipartition.split_as_int() for item in t1_nodes]
t2_node_splits = [item.bipartition.split_as_int() for item in t2_nodes]

t1_node_order = np.argsort(t1_node_splits)
t2_node_order = np.argsort(t2_node_splits)

t1_sorted = [t1_nodes[i] for i in t1_node_order]
t2_sorted = [t2_nodes[i] for i in t2_node_order]

print(np.array(t1_node_splits)[t1_node_order] - np.array(t2_node_splits)[t2_node_order])

t1_to_t2 = dict(zip(t1_sorted, t2_sorted))


agree = 0
disagree = 0

for k, v in t1_to_t2.items():

    if float(k.label) == float(v.label):
        agree += 1
    else:
        disagree += 1
        print("disagreement", float(k.label) - float(v.label))

print(agree, disagree) # disagree should be equal to zero
