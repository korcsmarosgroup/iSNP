import pandas as pd
import json
import pytest
from common_libs.mitab_handler import mitab_handler
from expected_structure import expected_dictionary


def json_reader(file):
    for line in open(file, mode="r"):
        yield json.loads(line)


@pytest.fixture(scope='session')
def new_parser():
    filepath = "example_mitab.tsv"
    parser = mitab_handler.MiTabHandler()
    parser.parse_mitab(filepath)
    return parser


def test_mitab_parser(new_parser):
    network = new_parser.network
    test_file = "example_mitab.tsv"

    with open(test_file, 'r') as e_file:
        expected = pd.read_csv(e_file, delimiter='\t', names=mitab_handler.mitab_header)
        expected = expected.fillna('-')

    pd.testing.assert_frame_equal(network, expected)
    assert isinstance(network, type(pd.DataFrame()))


@pytest.mark.parametrize("filename", ["tabsepfile.tsv", "mitabfile.mitab", "anotherfile.interactions"])
def test_any_format(tmpdir, filename):
    filename = tmpdir.join(filename)
    with open(filename.strpath, "w") as blank_file:
        blank_file.write("B\t"*41)

    handler = mitab_handler.MiTabHandler()
    handler.parse(filename, file_format='mitab')


def test_write_to_mitab(tmpdir, new_parser):
    expected_file = "example_mitab.tsv"
    tmp_path = tmpdir.join('test_mitab.tsv')
    new_parser.serialise_mitab(file_path=tmp_path.strpath, add_header=True)

    with open(expected_file, 'r') as e_file:
        expected = pd.read_csv(e_file, delimiter='\t', names=mitab_handler.mitab_header)
        expected = expected.fillna('-')

    with open(tmp_path.strpath, 'r') as a_file:
        actual = pd.read_csv(a_file, delimiter='\t')

    pd.testing.assert_frame_equal(actual, expected)


def test_sherlock_parsing(tmpdir):
    sherlock_json = "sherlock_format.json"
    handler = mitab_handler.MiTabHandler()
    handler.parse(sherlock_json, file_format='sherlock')

    assert expected_dictionary == handler.network.to_dict()


def test_sherlock_serialise(tmpdir):
    expected_path = "sherlock_format.json"
    tmp_path = tmpdir.join('network2json.json')
    handler = mitab_handler.MiTabHandler()
    handler.parse_sherlock(expected_path)
    handler.serialise_sherlock(tmp_path)

    ext = handler._parse_sherlock_structure(expected_path)
    act = handler._parse_sherlock_structure(tmp_path)

    edict = [ob for ob in ext]
    adict = [ob for ob in act]

    assert edict == adict
