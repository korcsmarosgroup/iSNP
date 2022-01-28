import pytest
import uniprot_id_formatter
import mock
import os


input_network_file = "example_files/test.tsv"
output_network_file = "example_files/results.tsv"
working_dir = f'{input_network_file.split("/")[0]}'


@mock.patch.object(uniprot_id_formatter, 'parse_args')
def test_input_network_file_exists(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = ['file', None, '-uc', '-ni', None]
        uniprot_id_formatter.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


@mock.patch.object(uniprot_id_formatter, 'parse_args')
def test_if_lower_case_and_upper_case_are_given(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [input_network_file, '-lc', '-uc', '-ni', None]
        uniprot_id_formatter.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 2


@mock.patch.object(uniprot_id_formatter, 'parse_args')
def test_case_one(args):
    args.return_value = [input_network_file, None, '-uc', '-ni', None]
    uniprot_id_formatter.main()

    num_lines_input = sum(1 for line in open(input_network_file))

    for filename in os.listdir(working_dir):
        key_word = "formatted"
        file = os.path.join(working_dir, filename)
        if key_word in file:
            num_lines_output = sum(1 for line in open(file))
            assert num_lines_input == num_lines_output


@mock.patch.object(uniprot_id_formatter, 'parse_args')
def test_case_two(args):
    args.return_value = [input_network_file, None, '-uc', '-ni', None]
    uniprot_id_formatter.main()

    int_a_output = ["mirbase:hsa-mir-99b-5p", "uniprotac:O95076", "uniprotac:P15172"]
    int_b_output = ["uniprotac:O95407", "uniprotac:O95150", "uniprotac:Q6S8E5"]

    for filename in os.listdir(working_dir):
        key_word = "formatted"
        file = os.path.join(working_dir, filename)
        if key_word in file:
            with open(file, 'r') as output:
                index = 0
                for line in output:
                    line = line.strip().split('\t')
                    assert line[0] == int_a_output[index]
                    assert line[1] == int_b_output[index]
                    index = index + 1


@mock.patch.object(uniprot_id_formatter, 'parse_args')
def test_case_three(args):
    args.return_value = [input_network_file, None, None, '-ni', None]
    uniprot_id_formatter.main()

    int_a_output = ["mirbase:hsa-mir-99b-5p", "uniprotac:o95076", "uniprotac:p15172"]
    int_b_output = ["uniprotac:o95407", "uniprotac:o95150", "uniprotac:q6s8e5"]

    for filename in os.listdir(working_dir):
        key_word = "formatted"
        file = os.path.join(working_dir, filename)
        if key_word in file:
            with open(file, 'r') as output:
                index = 0
                for line in output:
                    line = line.strip().split('\t')
                    assert line[0] == int_a_output[index]
                    assert line[1] == int_b_output[index]
                    index = index + 1


@mock.patch.object(uniprot_id_formatter, 'parse_args')
def test_case_four(args):
    args.return_value = [input_network_file, None, '-uc', None, None]
    uniprot_id_formatter.main()

    int_a_output = ["mirbase:hsa-mir-99b-5p", "uniprotac:O95076", "uniprotac:P15172-1"]
    int_b_output = ["uniprotac:O95407", "uniprotac:O95150-1", "uniprotac:Q6S8E5"]

    for filename in os.listdir(working_dir):
        key_word = "formatted"
        file = os.path.join(working_dir, filename)
        if key_word in file:
            with open(file, 'r') as output:
                index = 0
                for line in output:
                    line = line.strip().split('\t')
                    assert line[0] == int_a_output[index]
                    assert line[1] == int_b_output[index]
                    index = index + 1


@mock.patch.object(uniprot_id_formatter, 'parse_args')
def test_output_file_if_its_given(args):
    args.return_value = [input_network_file, None, '-uc', None, output_network_file]
    uniprot_id_formatter.main()

    assert len([name for name in os.listdir(working_dir)]) == 3
