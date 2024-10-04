""" Network Enrichment Module """
import sys

import os
import argparse
import numpy as np
import pandas as pd

sys.path.append("/rds/general/user/bbohar/home/projects/iSNP/analytic-modules")
from common_libs.mitab_handler import mitab_handler as h


class InvalidNumberOfHops(Exception):
    """ Raised if the distance is < 0 or > 100 or not an integer """
    pass


def parse_args(argv=None):
    help_text = \
        """
        === Network Enrichment ===
        Name of the tool: Network Enrichment.
        \nDescription:
        Using a reference network, for a given number of hops this module will 
        enrich the patient network with interactions, of order n hops, found 
        within the reference network. 

        \nParameters:
        --mirna <path to an existing file> [mandatory]
        --genomic <path to the new output file> [mandatory]
        --output <path to the output file for the enriched network> [mandatory]
        """

    # New argument Parser
    parser = argparse.ArgumentParser(description=help_text)

    # Input file path to network file
    parser.add_argument("-i", "--input",
                        help="<path to MITAB Network file> [mandatory]",
                        dest="input",
                        action="store",
                        required=True)

    # Path distance for enrichment
    parser.add_argument("-d", "--distance",
                        help="<number of hops for the enrichment> [mandatory]",
                        dest="distance",
                        type=int,
                        action="store",
                        required=True)

    # Reference Network
    parser.add_argument("-r", "--reference-net",
                        help="<path to an new file> [mandatory]",
                        dest="reference_network",
                        action="store",
                        required=True)

    # Output file path
    parser.add_argument("-o", "--output",
                        help="<path to an new file> [mandatory]",
                        dest="output",
                        action="store",
                        required=True)

    results = parser.parse_args(argv)

    return results


def _check_args(network, distance, reference_network, output):
    """
    A function to check the arguments parsed from the CLI

    Parameters
    ----------
    network: str, path to the network file
    distance: int, number of hops to enrich the network
    reference_network: str, path to the reference network
    output: str, path to the output file

    Raises
    -------
    FileNotFoundError: rasied if any of the paths do not exisit
    InvalidNumberOfHops: raised if the distance is < 0 or > 100 or not an integer

    """
    if not os.path.exists(network):
        raise FileNotFoundError

    if not os.path.exists(reference_network):
        raise FileNotFoundError

    if not isinstance(distance, int):
        raise  InvalidNumberOfHops

    if distance < 0:
        raise InvalidNumberOfHops

    if distance > 100:
        raise InvalidNumberOfHops


def get_neighbours(ref_network, vertex_list, n_neighbours=1):
    """

    Parameters
    ----------
    ref_network: str, path to the reference network
    vertex_list: list, a list of vertices to be considered by the search
    n_neighbours: int, current level neighbours to consider

    Returns
    -------
    new_edges: pandas dataframe, a dataframe holding the additional nodes

    """

    inter_a = h.mitab_header[0]
    inter_b = h.mitab_header[1]

    vertex_frame = pd.DataFrame(vertex_list, columns=[inter_a])
    new_edges = pd.DataFrame.merge(vertex_frame, ref_network, how="left")
    new_edges = new_edges.dropna()
    n_neighbours -= 1

    if n_neighbours > 0:
        new_vertex = list(set(new_edges[inter_b]))
        edges_2 = get_neighbours(ref_network, new_vertex, n_neighbours)
        new_edges = pd.DataFrame.append(new_edges, edges_2)

    return new_edges


def enrich_network(network, reference_network, distance):
    """

    Parameters
    ----------
    network: str, path to the patient network file
    reference_network: str, path to the reference network
    distance: int, number of hops to enrich the network

    Returns
    -------
    full_net: pandas Dataframe, new merged/enriched network

    """

    inter_a = h.mitab_header[0]
    inter_b = h.mitab_header[1]

    aa_vertex = pd.unique(network[[inter_a, inter_b]].values.ravel('K'))
    aa_vertex = list(set(aa_vertex))

    out = get_neighbours(reference_network, aa_vertex, distance)

    full_net = pd.DataFrame.append(network, out)
    full_net['interaction_id'] = full_net.apply(lambda row: "".join(sorted([row[inter_a], row[inter_b]])), axis=1)
    full_net = full_net.drop_duplicates(subset='interaction_id')
    del(full_net['interaction_id'])

    return full_net


def load_network(network_file):
    """
    A helper method to parse the network(s) to a usable dataframe

    Parameters
    ----------
    network_file: str, path to the network file

    Returns
    -------
    network: pandas Dataframe, holding the parsed network file

    """
    handler = h.MiTabHandler()
    handler.parse(network_file)
    return handler.network


def run(network, distance, reference_network, output):
    """
    This function that controls the logic.
    args --> check_args --> run (load_network + enrich_network + serialise) --> exit

    Parameters
    ----------
    network: str, path to the patient network file
    distance: int, number of hops to enrich the network
    reference_network: str, path to the reference network
    output: str, path to the output file

    """

    if os.stat(network).st_size == 0:

        print(f'====== The network fasta file is empty! ======')
        open(output, "a").close()

    else:

        try:
            net = load_network(network)
            ref_net = load_network(reference_network)
            enriched_network = enrich_network(net, ref_net, distance)
            new_handler = h.MiTabHandler()
            new_handler.network = enriched_network
            new_handler.serialise_mitab(output)
        except RuntimeError:
            raise RuntimeError()


def main(argv):
    """ Main method - waits for exit code """
    args = parse_args(argv)
    _check_args(args.input, args.distance, args.reference_network, args.output)
    run(args.input, args.distance, args.reference_network, args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
