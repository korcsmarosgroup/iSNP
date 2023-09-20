from pandas import DataFrame
from scipy import sparse
import networkx as nx
import pandas as pd
import numpy as np
import argparse
import random
import re
import sys


def parse_args(argv):
    """ Command line interface for the module """
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--graph',
                        help="<path to the ncol network file [mandatory]",
                        dest="ncol_file",
                        action="store",
                        required=True)

    parser.add_argument("-pf", "--patient_file",
                        help="<path to a tab separated patient file where the first column is the same ID as in the network nodes. The column headers are the patient Ids, the values are initial wheights for propagation. In the test example we have them as binary values (0-1) [mandatory]",
                        dest="patient_file",
                        action="store",
                        required=True)

    parser.add_argument("-of", "--output_file",
                        help="<path to the outputfile [mandatory]",
                        dest="output_file",
                        action="store",
                        required=True)

    parser.add_argument("-tf", "--tf_target_network",
                        help="transcription factor target gene network as an ncol file. If given the program calculates the TF-TG step. [optional]",
                        dest="tf_target_network",
                        action="store",
                        required=False)

    parser.add_argument("-rr", "--random_runs",
                        help="Random run number default 1000 [optional]",
                        dest="random_runs",
                        default=1000,
                        type=int,
                        action="store",
                        required=False)

    results = parser.parse_args(argv)
    return results


def calculate_alpha(edgenum, m = -0.02935302, b = 0.74842057):
    log_edge_count = np.log10(edgenum)
    alpha_val = round(m * log_edge_count + b, 3)

    if alpha_val <= 0:
        raise ValueError('Alpha <= 0 - Network Edge Count is too high')

    else:
        return alpha_val


def get_edge_weight(G):
    """
    Graph weight calculation based on outdegree. It overwrites the original weight attribute.
    :param G: networkx graph object
    :return: networkx graph object with degree based weight attribute.
    """
    for n in G:
        G.nodes[n]["out_degree"] = G.out_degree[n]
    for e in G.edges:
        G.edges[e]["weight"] = 1 / G.nodes[e[0]]["out_degree"]
    return G


def create_vector_from_node_ids(node_id_set, nx_graph, patient_series=False):
    """
    This function creates a row vector (list) which contain the nodes with SNP effected weights.
    Initially the SNP affected nodes has a weight of one.
    :param node_id_set:
    :param nx_graph:
    :param
    :return: Output file where you store the data.
    """
    start_vector = []
    for node in nx_graph.nodes:
        if node in node_id_set:
            if type(patient_series) == bool:
                start_vector.append(1)
            else:
                values = patient_series[node]
                try:
                    value = max(values)
                except:
                    value = values
                start_vector.append(value)
        else:
            start_vector.append(0)
    return start_vector


def calcaulate_network_propagation(ppi_graph, patient_df, outfile, random_runs, tf_tg_graph=None):
    """
    :param ppi_graph: networkx graph object
    :param patient_df: Dataframe of patients in the format of columns patient ids rows gene #IDs which are in the
    graph node IDs
    :param outfile: Output file where you store the data.
    :return:
    """
    number_of_nodes = len(set(ppi_graph.nodes))
    print("The number of nodes in the giant component are:", number_of_nodes)
    number_of_edges = len(ppi_graph.edges)
    print("The number of edges in the giant component are:", number_of_edges)
    ppi_names_set = set(ppi_graph.nodes)
    ppi_names_list = list(ppi_graph.nodes)

    # Create kernel
    ppi_graph = get_edge_weight(ppi_graph)
    A = np.array(nx.convert_matrix.to_numpy_array(ppi_graph, weight = 'weight'))
    alpha = calculate_alpha(number_of_edges)
    term2 = np.identity(A.shape[0]) - alpha * A
    term2_inv = np.linalg.inv(term2)
    term1 = (1.0 - alpha) * np.identity(A.shape[0])
    kernel = np.dot(term1, term2_inv)
    kernel_sparse = sparse.csr_matrix(kernel)
    print("Kernel calculated")

    # Patient specific datasets
    number_of_affected_proteins_set = set()
    patient_out_df_dic = {}

    for patient_id in patient_df.columns.tolist():
        print("This is patient:", patient_id)
        SNP_affected_nodes = set(patient_df[patient_df[patient_id] > 0].index.tolist())
        print("It has", len(SNP_affected_nodes), "SNPs.")
        source_nodes_in_graph = SNP_affected_nodes & ppi_names_set
        print("From them there are", len(source_nodes_in_graph), "in the graph")
        patient_out_df_dic[patient_id] = {}
        number_of_affected_proteins_set.add(len(source_nodes_in_graph))
        start_vector = create_vector_from_node_ids(source_nodes_in_graph, ppi_graph, patient_series=patient_df[patient_id])
        start_vector = np.array([start_vector])
        start_vector_sparse = sparse.csr_matrix(start_vector)
        patient_out_df_dic[patient_id]["propagation_vector"] = np.dot(start_vector_sparse, kernel_sparse)
        patient_out_df_dic[patient_id]["source_node_in_the_graph"] = len(source_nodes_in_graph)

    #  Random datasets
    print("Running random tests")
    random_out_dictionarry = {}

    for random_node_number in number_of_affected_proteins_set:
        random_out_dictionarry[random_node_number] = {}
        print("Running a case when", random_node_number, "nodes are in the graph.")

        for i in range(random_runs):
            # This sampling chooses from random nodes which are in the network.
            random_nodes = random.sample(ppi_names_list, random_node_number)
            start_vector = create_vector_from_node_ids(random_nodes, ppi_graph)
            start_vector = np.array ([start_vector])
            start_vector_sparse = sparse.csr_matrix(start_vector)
            random_out_dictionarry[random_node_number][i] = np.dot(start_vector_sparse, kernel_sparse)

    # Comparision
    print("Calculating Z scores")
    for patient in patient_out_df_dic:
        number_ref = patient_out_df_dic[patient]["source_node_in_the_graph"]
        node_id = 0

        for node in ppi_graph.nodes:
            patient_out_df_dic[patient][node] = {}
            rv = patient_out_df_dic[patient]["propagation_vector"][0,node_id]
            value_list = []

            for i in random_out_dictionarry[number_ref]:
                value_list.append(random_out_dictionarry[number_ref][i][0,node_id])
            mean = float(np.mean(value_list))
            std = float(np.std(value_list))

            if std == 0:
                Z = 0.0

            else:
                Z = (float(rv) - mean) / std

            patient_out_df_dic[patient][node]["real_value"] = rv
            patient_out_df_dic[patient][node]["Z"] = Z
            patient_out_df_dic[patient][node]["mean"] = mean
            patient_out_df_dic[patient][node]["std"] = std
            node_id = node_id + 1

    patient_out_df_dic["running_parameters"] = {}
    patient_out_df_dic["running_parameters"]["kernel"] = kernel
    patient_out_df_dic["running_parameters"]["alpha"] = alpha

    print("Writing out output")
    # Output writing and formatting.
    final_result = DataFrame.from_dict(patient_out_df_dic, orient="index").T
    final_result.to_csv(outfile, sep = "\t")

    # This formatting creates the Z scores above all nodes and the real values above all the nodes
    # The former can be used to build up the IBD PPI network the latter can be used for constructing the
    # the TF-TG propagation which is basically one more step based on the alpha heat.
    z_dic = {}
    rv_dic = {}

    for patient in patient_out_df_dic:

        if patient != "running_parameters":
            z_dic[patient] = {}
            rv_dic[patient] = {}

            for node_id in patient_out_df_dic[patient]:

                if (node_id != "propagation_vector") & (node_id != "source_node_in_the_graph"):
                    z_dic[patient][node_id] = patient_out_df_dic[patient][node_id]["Z"]
                    rv_dic[patient][node_id] = patient_out_df_dic[patient][node_id]["real_value"]

    final_result_z = DataFrame.from_dict(z_dic, orient="columns")
    final_result_z.to_csv(outfile.replace(".txt", "_Z.txt"), sep = "\t")
    final_result_rv = DataFrame.from_dict(rv_dic, orient="columns")
    final_result_rv.to_csv(outfile.replace(".txt", "_rv.txt"), sep = "\t")

    if tf_tg_graph:
        print ("Calculating TF-TG part of the pipline")
        tf_tg_graph = get_edge_weight(tf_tg_graph)
        tf_tg_dictionarry = {}
        print("Calculating patient specific values")

        for patient in patient_out_df_dic:

            if patient != "running_parameters":
                tf_tg_dictionarry[patient] = {}

                for u, v in tf_tg_graph.edges:

                    if u in ppi_graph.nodes:
                        tf_real_weight_here = (1-alpha) * patient_out_df_dic[patient][u]["real_value"]*tf_tg_graph.edges[u,v]["weight"]

                        if v not in tf_tg_dictionarry[patient]:
                            tf_tg_dictionarry[patient][v] = {}
                            tf_tg_dictionarry[patient][v]["TF_TG_weight"] = tf_real_weight_here

                        else:
                            tf_tg_dictionarry[patient][v]["TF_TG_weight"] = tf_tg_dictionarry[patient][v]["TF_TG_weight"]\
                                                                             + tf_real_weight_here
        #Random runs
        random_tf_out_dictionarry = {}
        for random_node_number in number_of_affected_proteins_set:
            random_tf_out_dictionarry[random_node_number] = {}
            print("Running a case when", random_node_number, "nodes are in the TF-TG graph.")
            node_id = 0 # This is the node's index in the ppi_graph

            for u in ppi_graph.nodes:

                if u in tf_tg_graph.nodes:

                    if tf_tg_graph.nodes[u]["out_degree"] > 0:

                        for v in tf_tg_graph[u]:

                            if v not in random_tf_out_dictionarry[random_node_number]:
                                random_tf_out_dictionarry[random_node_number][v] = [0] * random_runs

                            for i in range(random_runs):
                                random_tf_weight = random_out_dictionarry[random_node_number][i][0, node_id]
                                weight_here = random_tf_weight * tf_tg_graph.edges[u,v]["weight"] * (1-alpha)
                                random_tf_out_dictionarry[random_node_number][v][i] = random_tf_out_dictionarry[random_node_number][v][i]+weight_here
                node_id=node_id+1

        # Comparison
        for patient in tf_tg_dictionarry:
            number_ref = patient_out_df_dic[patient]["source_node_in_the_graph"]

            for v in tf_tg_dictionarry[patient]:
                rv = tf_tg_dictionarry[patient][v]["TF_TG_weight"]
                value_list = random_tf_out_dictionarry[number_ref][v]
                mean = float(np.mean(value_list))
                std = float(np.std(value_list))

                if std == 0:
                    Z = 0.0

                else:
                    Z = (float(rv) - mean) / std

                tf_tg_dictionarry[patient][v]["Z"] = Z
                tf_tg_dictionarry[patient][v]["mean"] = mean
                tf_tg_dictionarry[patient][v]["std"] = std

        # to avoid memory issues:
        random_tf_out_dictionarry = 0

        #Output_wrting
        final_result = DataFrame.from_dict(tf_tg_dictionarry, orient="index").T
        final_result.to_csv(outfile.replace(".txt","TF_TG.txt"), sep = "\t")

        z_dic = {}
        rv_dic = {}

        for patient in tf_tg_dictionarry:

            if patient != "running_parameters":
                z_dic[patient] = {}
                rv_dic[patient] = {}

                for node_id in tf_tg_dictionarry[patient]:

                    if (node_id != "propagation_vector") & (node_id != "source_node_in_the_graph"):
                        z_dic[patient][node_id] = tf_tg_dictionarry[patient][node_id]["Z"]
                        rv_dic[patient][node_id] = tf_tg_dictionarry[patient][node_id]["TF_TG_weight"]

        final_result_z = DataFrame.from_dict(z_dic, orient="columns")
        final_result_z.to_csv(outfile.replace(".txt", "_Z_TF.txt"), sep = "\t")
        final_result_rv = DataFrame.from_dict(rv_dic, orient="columns")
        final_result_rv.to_csv(outfile.replace(".txt", "_rv_TF.txt"), sep = "\t")


def main(argv):
    args = parse_args(argv)
    input_df = pd.read_csv(args.patient_file, index_col = 0, sep = "\t", header = 0)
    graph = nx.read_edgelist(args.ncol_file, delimiter = " ", create_using = nx.DiGraph)

    # Creating giant component Every work is happening only on this graph.
    number_of_nodes = len(set(graph.nodes))
    print("The number of nodes in the graph are:", number_of_nodes)
    number_of_edges = len(graph.edges)
    print("The number of edges in the graph are:", number_of_edges)
    wg_nodes = max(nx.weakly_connected_components(graph), key = len)
    work_graph = graph.subgraph(wg_nodes)

    print("The overlap of the graph and the input data frame ids is:")
    print(len(set(input_df.index) & set(work_graph.nodes)))

    if len(set(input_df.index) & set(work_graph.nodes)) == 0:
        raise ValueError("The input graph and input dataframe has different IDs. Please check the IDs.")

    if not (len(re.findall(".tsv", args.output_file)) > 0):
        raise NameError("The outfile extension must be .tsv")

    random_runs = int(args.random_runs)

    if args.tf_target_network:
        tf_tg_graph = nx.read_edgelist(args.tf_target_network, delimiter = " ", create_using = nx.DiGraph)

        if len(set(tf_tg_graph.nodes) & set(work_graph.nodes)) == 0:
            raise(ValueError("The input ppi graph and TF-TG graph has different IDs. Please check the IDs."))

    else:
        tf_tg_graph = None

    calcaulate_network_propagation(work_graph, input_df, args.output_file, random_runs, tf_tg_graph)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
