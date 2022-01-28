""" Network Enrichment Module Tests """
import os

import network_enrichment as ne
import pandas as pd
import pytest

from common_libs.mitab_handler import mitab_handler as h

# Expected results
expected_num_neighbours = [8, 9, 11]

inter_a = h.mitab_header[0]
inter_b = h.mitab_header[1]


@pytest.fixture()
def create_example_networks(tmpdir):
    """ Create a example test networks with redunent nodes and varying degrees """

    inter_a = h.mitab_header[0]
    inter_b = h.mitab_header[1]
    alt_a = h.mitab_header[2]
    alt_b = h.mitab_header[3]

    # Create base network
    network = pd.DataFrame({inter_a: ["A", "A", "A", "A", "A", "B", "X"],
                            inter_b: ["B", "C", "D", "E", "F", "A", "W"],
                            alt_a: ["A1", "A1", "A1", "A1", "A1", "B1", "X1"],
                            alt_b: ["B2", "C2", "D2", "E2", "F2", "A1", "W1"]})


    # Create reference network
    reference_network = pd.DataFrame({inter_a: ["F", "A", "1", "4", "4", "2", "5", "A", "123", "X"],
                                      inter_b: ["D", "1", "4", "5", "6", "3", "B", "C", "321", "W"],
                                      alt_a: ["alt:F", "alt:A1", "alt:1", "alt:4", "alt:4", "alt:2", "alt:5", "alt:A", "alt:123", "alt:X"],
                                      alt_b: ["alt:D", "alt:1", "alt:4", "alt:5", "alt:6", "alt:3", "alt:B", "alt:C", "alt:321", "alt:W"]})

    # Write the test networks to a local temp file
    network_path = os.path.join(str(tmpdir), 'network.tsv')
    reference_path = os.path.join(str(tmpdir), "reference_path.tsv")

    # Write the example patient network
    network_mitab = h.MiTabHandler()
    network_mitab.network[inter_a] = network[inter_a]
    network_mitab.network[inter_b] = network[inter_b]
    network_mitab.network[alt_a] = network[alt_a]
    network_mitab.network[alt_b] = network[alt_b]
    network_mitab.serialise_mitab(network_path)

    # Write the example reference network
    reference_mitab = h.MiTabHandler()
    reference_mitab.network[inter_a] = reference_network[inter_a]
    reference_mitab.network[inter_b] = reference_network[inter_b]
    reference_mitab.network[alt_a] = network[alt_a]
    reference_mitab.network[alt_b] = network[alt_b]
    reference_mitab.serialise_mitab(reference_path)

    return network_mitab.network, reference_mitab.network, network_path, reference_path


def load_network(network):
    """ Load the example networks into Navi-MiTab formats """
    handler = h.MiTabHandler()
    handler.network = network
    return handler.network


def test_failed_file():
    """ A test to ensure that incorrect file types can not be used in the module """
    wrong_file_path = "clearly_this_file_is_wrong.html"

    with pytest.raises(FileNotFoundError) as pytest_wrapped_e:
        ne._check_args(network=wrong_file_path,
                       distance=1,
                       reference_network=wrong_file_path,
                       output=wrong_file_path)

    assert pytest_wrapped_e.type == FileNotFoundError


@pytest.mark.parametrize('invalid_hops', [-12, "two", 10000])
def test_invalid_number_hops(tmpdir, create_example_networks, invalid_hops):
    """ A test to ensure that only valid number of hops can be inputted """

    output_file = os.path.join(str(tmpdir), "never_used_network_file.tsv")
    _, _, network_path, reference_path = create_example_networks

    with pytest.raises(ne.InvalidNumberOfHops) as pytest_wrapped_e:
        ne._check_args(network=network_path,
                       distance=invalid_hops,
                       reference_network=reference_path,
                       output=output_file)

    assert pytest_wrapped_e.type == ne.InvalidNumberOfHops


def test_duplicates_removal(tmpdir, create_example_networks):
    """ A component test to ensure the removal of duplicates from the network after enrichment """
    output_file = os.path.join(str(tmpdir), "example_enriched.tsv")
    network, reference_network, network_path, reference_path = create_example_networks

    # Run network enc
    ne.run(network=network_path,
           distance=3,
           reference_network=reference_path,
           output=output_file)

    # Find number of neighbours added
    output_frame = pd.read_csv(output_file, delimiter="\t", header=None)
    check_duplicates = output_frame.loc[(output_frame[0] == "b") & (output_frame[1] == "a")]
    print(output_file)
    assert len(check_duplicates) == 0


def test_duplicates_removal_keeps_original_order(tmpdir, create_example_networks):
    """ A component test to ensure that the removal of duplicates from the network after enrichment will not change the order of the interactors """
    output_file = os.path.join(str(tmpdir), "example_enriched.tsv")
    network, reference_network, network_path, reference_path = create_example_networks

    # Run network enc
    ne.run(network=network_path,
           distance=1,
           reference_network=reference_path,
           output=output_file)

    # Find number of neighbours added
    output_frame = pd.read_csv(output_file, delimiter="\t", header=None)
    check_duplicates = output_frame.loc[(output_frame[0] == "x") & (output_frame[1] == "w")]

    assert len(check_duplicates) == 1


@pytest.mark.parametrize('number_of_hops', [1, 2, 3])
def test_enrich_network(tmpdir, create_example_networks, number_of_hops):
    """ A component test to check the functionality of enrichment module """

    output_file = os.path.join(str(tmpdir), "example_enriched.tsv")
    network, reference_network, network_path, reference_path = create_example_networks

    # Run network enc
    ne.run(network=network_path,
           distance=number_of_hops,
           reference_network=reference_path,
           output=output_file)

    # Find number of neighbours added
    output_frame = pd.read_csv(output_file, delimiter="\t", header=None)
    enrichment_count = output_frame.count(axis=0)

    assert "123" not in output_frame.values
    assert "321" not in output_frame.values
    assert enrichment_count[0] == expected_num_neighbours[number_of_hops-1]


def test_enrich_network_alternatives(tmpdir, create_example_networks):
    """ A component test to check the functionality of enrichment module and alternative ids """

    output_file = os.path.join(str(tmpdir), "example_enriched.tsv")
    network, reference_network, network_path, reference_path = create_example_networks

    # Run network enc
    ne.run(network=network_path,
           distance=2,
           reference_network=reference_path,
           output=output_file)

    # Find number of neighbours added
    output_frame = pd.read_csv(output_file, delimiter="\t", header=None)

    assert "123" not in output_frame.values
    assert "321" not in output_frame.values
    assert len(output_frame[~output_frame.iloc[:, 2].str.contains('-')]) == expected_num_neighbours[1]

