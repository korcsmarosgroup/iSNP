
from collections import defaultdict
from pandas import DataFrame
from time import strftime
from scipy import sparse
import networkx as nx
import pandas as pd
import numpy as np
import argparse
import tempfile
import logging
import random
import sys
import re
import os


if os.path.isfile("Network_propagation.log"):
    os.remove("Network_propagation.log")

logging.basicConfig(filename = 'Network_propagation.log', level = logging.INFO)


def parse_args(argv):
    """ Command line interface for the module """
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--graph',
                        help="<path to the ncol network file [mandatory]",
                        dest="ncol_file",
                        action="store",
                        required=True)

    parser.add_argument("-pf", "--patient_file",
                        help="<path to a tab separated patient file where the first column is the same ID as in the "
                             "network nodes. The column headers are the patient Ids, the values are initial wheights "
                             "for propagation. In the test example we have them as binary values (0-1) [mandatory]",
                        dest="patient_file",
                        action="store",
                        required=True)

    parser.add_argument("-of", "--output_file",
                        help="<path to the outputfile [mandatory]",
                        dest="output_file",
                        action="store",
                        required=True)

    parser.add_argument("-tf", "--tf_target_network",
                        help="transcription factor target gene network as an ncol file. If given the program calultes "
                             "the TF-TG step.",
                        dest="tf_target_network",
                        action="store",
                        required=False)

    parser.add_argument("-rr", "--random_runs",
                        help="Random run number default 1000",
                        dest="random_runs",
                        default=1000,
                        type=int,
                        action="store",
                        required=False)

    parser.add_argument("-dg","--degree_propagation",
                        help= "Runs degree kept random proapgation, default = False",
                        dest="degree_propagation",
                        default=False,
                        type=bool,
                        required=False,
                        action="store")

    parser.add_argument("-dgsw", "--sortwide",
                        help="From the degree sorrting how much percentage up or down of the sorted nodes take the "
                             "random choice, default = 5",
                        dest="sortwide",
                        default = 5,
                        type=float,
                        required=False,
                        action="store")
    results = parser.parse_args(argv)
    return results


def calculate_alpha(edgenum, m = -0.02935302, b = 0.74842057):
    log_edge_count = np.log10(edgenum)
    alpha_val = round(m * log_edge_count + b, 3)

    if alpha_val <= 0:
        raise ValueError('Alpha <= 0 - Network Edge Count is too high')

    # There should never be a case where Alpha >= 1, as avg node degree will never be negative
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


def create_vector_from_node_ids(node_id_set, nx_graph, patient_series = False):
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


def calculate_network_propagation(
        ppi_graph,
        patient_df,
        outfile,
        random_runs,
        tf_tg_graph=None,
        degree_based_random_run = False,
        sortwide = 5
    ):
    """
    :param ppi_graph: networkx graph object
    :param patient_df: Dataframe of patients in the format of columns patient ids rows gene #IDs which are in the
    graph node IDs
    :param outfile: Output file where you store the data.
    :return:
    """
    number_of_nodes = len(set(ppi_graph.nodes))
    logging.info(f"### [{strftime('%H:%M:%S')}] The number of uniqe nodes in the giant component are: {number_of_nodes}")
    number_of_edges = len(ppi_graph.edges)
    logging.info(f"### [{strftime('%H:%M:%S')}] The number of uniqe edges in the giant component are: {number_of_edges}")
    ppi_names_set = set(ppi_graph.nodes)
    ppi_names_list = list(ppi_graph.nodes)

    # Create kernel
    ppi_graph = get_edge_weight(ppi_graph)
    A = np.array(nx.convert_matrix.to_numpy_array(ppi_graph, weight='weight'))
    alpha = calculate_alpha(number_of_edges)
    logging.info(f"### [{strftime('%H:%M:%S')}] Alpha: {alpha}") # calculated based on the preset of the Huang et al paper
    term2 = np.identity(A.shape[0]) - alpha * A
    term2_inv = np.linalg.inv(term2)
    term1 = (1.0 - alpha) * np.identity(A.shape[0])
    kernel = np.dot(term1, term2_inv)
    kernel_sparse = sparse.csr_matrix(kernel)
    logging.info(f"### [{strftime('%H:%M:%S')}] Kernel calculated")

    patinet_out_df_dic = {}

    patinet_out_df_dic["running_parameters"] = {}
    patinet_out_df_dic["running_parameters"]["kernel"] = kernel
    patinet_out_df_dic["running_parameters"]["alpha"] = alpha

    if degree_based_random_run == False:
        # Patient specific datasets
        number_of_affected_proteins_set = set()

        for patient_id in patient_df.columns.tolist():
            logging.info(f"### [{strftime('%H:%M:%S')}] This is patient: {patient_id}")
            SNP_affected_nodes = set(patient_df[patient_df[patient_id] > 0].index.tolist())
            logging.info(f"### [{strftime('%H:%M:%S')}] It has {len(SNP_affected_nodes)} SNPs.")
            source_nodes_in_graph = SNP_affected_nodes & ppi_names_set
            logging.info(f"### [{strftime('%H:%M:%S')}] From them there are {len(source_nodes_in_graph)} in the graph.")
            patinet_out_df_dic[patient_id] = {}
            number_of_affected_proteins_set.add(len(source_nodes_in_graph))
            start_vector = create_vector_from_node_ids(source_nodes_in_graph, ppi_graph,
                                                       patient_series=patient_df[patient_id])
            start_vector = np.array([start_vector])
            start_vector_sparse = sparse.csr_matrix(start_vector)
            patinet_out_df_dic[patient_id]["propagation_vector"] = np.dot(start_vector_sparse, kernel_sparse)
            patinet_out_df_dic[patient_id]["source_node_in_the_graph"] = len(source_nodes_in_graph)

        # Random datasets
        logging.info(f"### [{strftime('%H:%M:%S')}] Running random tests")
        random_out_dictionarry = {}
        for random_node_number in number_of_affected_proteins_set:
            random_out_dictionarry[random_node_number] = {}
            logging.info(f"### [{strftime('%H:%M:%S')}] Running a case when {random_node_number} nodes are in the graph.")

            for i in range(random_runs):
                # This sampling chooses from random nodes which are in the network.
                random_nodes = random.sample(ppi_names_list,
                                             random_node_number)
                start_vector = create_vector_from_node_ids(random_nodes, ppi_graph)
                start_vector = np.array([start_vector])
                start_vector_sparse = sparse.csr_matrix(start_vector)
                random_out_dictionarry[random_node_number][i] = np.dot(start_vector_sparse, kernel_sparse)

        # Comparision
        logging.info(f"### [{strftime('%H:%M:%S')}] Calculating Z scores")
        nodes = sorted(list(ppi_graph.nodes))

        for patient, patient_data in patinet_out_df_dic.items():

            if patient == "running_parameters":
                continue

            number_ref = patient_data["source_node_in_the_graph"]
            propagation_vector = patient_data["propagation_vector"]
            real_values = propagation_vector[0].toarray().flatten()

            random_outputs = random_out_dictionarry[number_ref]
            random_matrix = np.vstack([
                random_outputs[i][0].toarray().flatten() for i in random_outputs
            ])

            means = np.mean(random_matrix, axis=0)
            stds = np.std(random_matrix, axis=0)
            Z_scores = np.where(stds == 0, 0.0, (real_values - means) / stds)

            for node, rv, z, mean, std in zip(nodes, real_values, Z_scores, means, stds):
                patient_data[node] = {
                    "real_value": float(rv),
                    "Z": float(z),
                    "mean": float(mean),
                    "std": float(std)
                }

    if degree_based_random_run == True:
        # creating network specifc dataframe

        all_out_degree = nx.get_node_attributes(ppi_graph, "out_degree")
        all_out_degree_df = pd.DataFrame.from_dict(all_out_degree, orient='index', columns=["out_degree"])
        all_out_degree_df.sort_values(by=["out_degree"], inplace=True)
        all_out_degree_df["rank"] = all_out_degree_df["out_degree"].rank(method='first')

        # Cache for degree to node list mapping
        degree_to_nodes = all_out_degree_df.groupby("out_degree").apply(lambda df: df.index.tolist()).to_dict()
        ranked_nodes = all_out_degree_df.reset_index().rename(columns={"index": "node"})

        rank_difference = round(all_out_degree_df.shape[0] * (float(sortwide) / 100.0))
        random_out_dictionarry = {}
        patinet_out_df_dic = {}

        ppi_node_list = list(ppi_graph.nodes)

        for patient_id in patient_df.columns:
            logging.info(f"### [{strftime('%H:%M:%S')}] Patient: {patient_id}")
            SNP_nodes = set(patient_df.index[patient_df[patient_id] > 0])
            graph_nodes = SNP_nodes & ppi_names_set
            logging.info(f"### [{strftime('%H:%M:%S')}] SNPs in graph: {len(graph_nodes)}")

            patinet_out_df_dic[patient_id] = {}

            # Patient propagation
            start_vector = create_vector_from_node_ids(graph_nodes, ppi_graph, patient_series=patient_df[patient_id])
            start_vector_sparse = sparse.csr_matrix([start_vector])
            propagation_vector = start_vector_sparse @ kernel_sparse
            patinet_out_df_dic[patient_id]["propagation_vector"] = propagation_vector
            patinet_out_df_dic[patient_id]["source_node_in_the_graph"] = len(graph_nodes)

            # Precompute out-degrees
            patient_nodes_outdegree = dict(ppi_graph.out_degree(graph_nodes))
            selected_nodes_list = []

            for node, degree in patient_nodes_outdegree.items():

                if degree in degree_to_nodes and len(degree_to_nodes[degree]) >= 2 * rank_difference:
                    selected_nodes_list.append(degree_to_nodes[degree])

                else:
                    node_rank = all_out_degree_df.loc[node, "rank"]
                    min_rank = max(node_rank - rank_difference, 0)
                    max_rank = min(node_rank + rank_difference, all_out_degree_df.shape[0])
                    mask = (ranked_nodes["rank"] >= min_rank) & (ranked_nodes["rank"] <= max_rank)
                    selected_nodes = ranked_nodes[mask]["node"].tolist()
                    selected_nodes_list.append(selected_nodes)

            # Random propagation
            logging.info(f"### [{strftime('%H:%M:%S')}] Running random tests")
            random_matrix = np.zeros((random_runs, len(ppi_node_list)))

            for i in range(random_runs):
                random_nodes = [random.choice(nodes) for nodes in selected_nodes_list]
                rand_vector = create_vector_from_node_ids(random_nodes, ppi_graph)
                rand_sparse = sparse.csr_matrix([rand_vector])
                result = rand_sparse @ kernel_sparse
                random_matrix[i, :] = result.toarray()

            random_out_dictionarry[patient_id] = random_matrix

            # Calculate Z-scores
            real_values = propagation_vector.toarray().flatten()
            means = random_matrix.mean(axis=0)
            stds = random_matrix.std(axis=0)

            z_scores = np.where(stds == 0, 0.0, (real_values - means) / stds)

            for idx, node in enumerate(ppi_node_list):
                patinet_out_df_dic[patient_id][node] = {
                    "real_value": real_values[idx],
                    "Z": z_scores[idx],
                    "mean": means[idx],
                    "std": stds[idx],
                }

    logging.info(f"### [{strftime('%H:%M:%S')}] Writing out the output")
    # Output writing and formatting.
    final_result = pd.DataFrame.from_dict(patinet_out_df_dic, orient = "index").T
    final_result.to_csv(outfile, sep="\t")

    # This formatting creates the Z scores above all nodes and the real values above all the nodes
    # The former can be used to build up the IBD PPI network the latter can be used for constructing the
    # the TF-TG propagation which is basically one more step based on the alpha heat.
    z_dic = {}
    rv_dic = {}
    for patient in patinet_out_df_dic:

        if patient != "running_parameters":
            z_dic[patient] = {}
            rv_dic[patient] = {}

            for node_id in patinet_out_df_dic[patient]:

                if (node_id != "propagation_vector") & (node_id != "source_node_in_the_graph"):
                    z_dic[patient][node_id] = patinet_out_df_dic[patient][node_id]["Z"]
                    rv_dic[patient][node_id] = patinet_out_df_dic[patient][node_id]["real_value"]

    final_result_z = DataFrame.from_dict(z_dic, orient = "columns")
    final_result_z.to_csv(outfile.replace(".txt", "_Z.txt"), sep = "\t", index_label = "GeneID")

    final_result_rv = DataFrame.from_dict(rv_dic, orient = "columns")
    final_result_rv.to_csv(outfile.replace(".txt", "_rv.txt"), sep = "\t", index_label = "GeneID")

    if tf_tg_graph:
        logging.info(f"### [{strftime('%H:%M:%S')}] Calculating TF-TG part of the pipeline")
        tf_tg_graph = get_edge_weight(tf_tg_graph)
        tf_tg_dictionarry = {}
        logging.info(f"### [{strftime('%H:%M:%S')}] Calculating patinet specific values")

        for patient in patinet_out_df_dic:

            if patient != "running_parameters":
                tf_tg_dictionarry[patient] = {}

                for u, v in tf_tg_graph.edges:

                    if u in ppi_graph.nodes:
                        tf_real_weight_here = (1-alpha) * patinet_out_df_dic[patient][u]["real_value"]*\
                                               tf_tg_graph.edges[u,v]["weight"]

                        if v not in tf_tg_dictionarry[patient]:
                            tf_tg_dictionarry[patient][v] = {}
                            tf_tg_dictionarry[patient][v]["TF_TG_weight"] = tf_real_weight_here

                        else:
                            tf_tg_dictionarry[patient][v]["TF_TG_weight"] = tf_tg_dictionarry[patient][v]["TF_TG_weight"]\
                                                                             + tf_real_weight_here

        random_tf_out_dictionarry = {}

        if degree_based_random_run == False:
            print("FALSE")

            tf_nodes = [u for u in ppi_graph.nodes if u in tf_tg_graph and tf_tg_graph.nodes[u]["out_degree"] > 0]

            for random_node_number in number_of_affected_proteins_set:
                logging.info(f"### [{strftime('%H:%M:%S')}] Running a case when {random_node_number} nodes are in the TF-TG graph")
                tf_node_indices = {u: idx for idx, u in enumerate(tf_nodes)}
                
                # Preload random_out matrix: shape = (random_runs, num_tf_nodes)
                random_out_matrix = np.array([
                    [random_out_dictionarry[random_node_number][i][0, tf_node_indices[u]] for u in tf_nodes]
                    for i in range(random_runs)
                ])  # shape: (random_runs, len(tf_nodes))
                
                # Init output dict
                random_tf_out = defaultdict(lambda: np.zeros(random_runs))
                
                # Propagate signal through TF-TG graph
                for node_id, u in enumerate(tf_nodes):
                    tf_weights = random_out_matrix[:, node_id]  # shape: (random_runs,)
                    for v in tf_tg_graph[u]:
                        edge_weight = tf_tg_graph.edges[u, v]["weight"] * (1 - alpha)
                        random_tf_out[v] += tf_weights * edge_weight
                
                # Store as normal dict with lists
                random_tf_out_dictionarry[random_node_number] = {k: v.tolist() for k, v in random_tf_out.items()}

            # --- Z-score Computation ---
            for patient, tf_data in tf_tg_dictionarry.items():
                number_ref = patinet_out_df_dic[patient]["source_node_in_the_graph"]
                random_tf_out = random_tf_out_dictionarry[number_ref]
                
                for v, data in tf_data.items():
                    rv = data["TF_TG_weight"]
                    values = np.array(random_tf_out.get(v, [0.0] * random_runs))
                    mean = values.mean()
                    std = values.std()
                    Z = 0.0 if std == 0 else (rv - mean) / std

                    data["Z"] = Z
                    data["mean"] = mean
                    data["std"] = std

            # Free memory
            random_tf_out_dictionarry = 0

        if degree_based_random_run == True:
            print("TRUE")

            node_to_index = {node: idx for idx, node in enumerate(ppi_graph.nodes)}

            for patient, patient_data in patinet_out_df_dic.items():

                if patient == "running_parameters":
                    continue

                logging.info(f"### [{strftime('%H:%M:%S')}] Running degree-based random TF-TG test for patient {patient}")
                
                # Initialize random TF-TG results storage
                tf_output = defaultdict(lambda: np.zeros(random_runs))

                # Extract all patient random propagation vectors (shape: [random_runs, num_nodes])
                patient_random_matrix = np.vstack([
                    random_out_dictionarry[patient][i].flatten() for i in range(random_runs)
                ])  # Shape: (random_runs, num_nodes)

                # Only consider nodes in both graphs with non-zero out-degree in TF-TG
                for u in ppi_graph.nodes:

                    if u not in tf_tg_graph or tf_tg_graph.nodes[u]["out_degree"] == 0:
                        continue

                    u_idx = node_to_index[u]

                    # Get random weights for this node across all runs
                    random_weights_u = patient_random_matrix[:, u_idx]  # Shape: (random_runs,)

                    for v in tf_tg_graph[u]:
                        edge_weight = tf_tg_graph.edges[u, v]["weight"] * (1 - alpha)
                        tf_output[v] += random_weights_u * edge_weight  # Vectorized addition

                # Compute Z-scores for each target gene `v`
                for v, real_val in tf_tg_dictionarry[patient].items():
                    observed = real_val["TF_TG_weight"]
                    random_values = tf_output[v]
                    mean = random_values.mean()
                    std = random_values.std()
                    z = 0.0 if std == 0 else (observed - mean) / std

                    tf_tg_dictionarry[patient][v]["Z"] = z
                    tf_tg_dictionarry[patient][v]["mean"] = mean
                    tf_tg_dictionarry[patient][v]["std"] = std

            # Optional: clear large dict to free memory
            random_tf_out_dictionarry = None

        # Output wrting
        logging.info(f"### [{strftime('%H:%M:%S')}] Writing TF-TG output")
        final_result = DataFrame.from_dict(tf_tg_dictionarry, orient = "index").T
        final_result.to_csv(outfile.replace(".txt", "TF_TG.txt"), sep = "\t", index_label = "GeneID")

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

        final_result_z = DataFrame.from_dict(z_dic, orient = "columns")
        final_result_z.to_csv(outfile.replace(".txt", "_Z_TF.txt"), sep = "\t", index_label = "GeneID")

        final_result_rv = DataFrame.from_dict(rv_dic, orient = "columns")
        final_result_rv.to_csv(outfile.replace(".txt", "_rv_TF.txt"), sep = "\t", index_label = "GeneID")


def main(argv):
    logging.info(f"### [{strftime('%H:%M:%S')}] Starting the script")
    args = parse_args(argv)
    input_df = pd.read_csv(args.patient_file, index_col=0, sep="\t", header=0)

    with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_file:
        temp_file_path = temp_file.name

        with open(args.ncol_file, 'r') as original_file:
            for line in original_file:
                cleaned_line = line.rstrip()
                temp_file.write(cleaned_line + '\n')

    graph = nx.read_edgelist(temp_file_path, delimiter = " ", create_using = nx.DiGraph)

    # Creating giant component Every work is happening only on this graph.
    number_of_nodes = len(set(graph.nodes))
    logging.info(f"### [{strftime('%H:%M:%S')}] The number of nodes in the graph are: {number_of_nodes}")
    number_of_edges = len(graph.edges)
    logging.info(f"### [{strftime('%H:%M:%S')}] The number of edges in the graph are: {number_of_edges}")
    wg_nodes = max(nx.weakly_connected_components(graph), key=len)
    work_graph = graph.subgraph(wg_nodes)

    logging.info(f"### [{strftime('%H:%M:%S')}] The overlap of the graph and the input data frame ids is: "
                 f"{len(set(input_df.index) & set(work_graph.nodes))}")

    if len(set(input_df.index) & set(work_graph.nodes)) == 0:
        raise ValueError("The input graph and input dataframe has different IDs. Please check the IDs.")

    if not (len(re.findall(".txt", args.output_file)) > 0):
        raise NameError("The outfile extension must be .txt")

    random_runs = int(args.random_runs)

    if args.tf_target_network:

        with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as temp_file:
            temp_file_path = temp_file.name

            with open(args.tf_target_network, 'r') as tf_target_original:
                for line in tf_target_original:
                    cleaned_line = line.rstrip()
                    temp_file.write(cleaned_line + '\n')

        tf_tg_graph = nx.read_edgelist(temp_file_path, delimiter = " ", create_using = nx.DiGraph)

        if len(set(tf_tg_graph.nodes) & set(work_graph.nodes)) == 0:
            raise(ValueError("The input ppi graph and TF-TG graph has different IDs. Please check the IDs."))

    else:
        tf_tg_graph = None

    logging.info(f"### [{strftime('%H:%M:%S')}] Starting the Network Propagation")
    calculate_network_propagation(
        work_graph,
        input_df,
        args.output_file,
        random_runs,
        tf_tg_graph,
        degree_based_random_run = args.degree_propagation,
        sortwide = args.sortwide
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
