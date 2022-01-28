import pytest
import network_id_mapper
import mock
import os


example_files_list = ["example_files/file1.tsv",
                      "example_files/mapping_file.json,example_files/mapping_file2.json",
                      "example_files/file2.tsv",
                      "example_files/mapping_file3.json"]


@mock.patch.object(network_id_mapper, 'parse_args')
def test_input_network_file_exists(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = ['file', '', '', '', example_files_list[1], '', False]
        network_id_mapper.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


@mock.patch.object(network_id_mapper, 'parse_args')
def test_remove_parameter_is_good(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [example_files_list[0], 'asdasd', '', '', example_files_list[1], '', False]
        network_id_mapper.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 2


@mock.patch.object(network_id_mapper, 'parse_args')
def test_mapping_files_exists(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [example_files_list[0], 'true', '', '', 'fake_file.txt', '', False]
        network_id_mapper.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 3


@mock.patch.object(network_id_mapper, 'parse_args')
def test_case_one(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    args.return_value = [example_files_list[0], 'true', '', 'uniprotac', example_files_list[1], '', output_file]
    network_id_mapper.main()

    num_lines = sum(1 for line in open(output_file))

    protein_number = 0
    ensg_number = 0
    with open(output_file) as output:
        for line in output:
            protein = line.find('p00007')
            ensg = line.find('ensg')
            if protein != -1:
                protein_number = protein_number + 1
            if ensg != -1:
                ensg_number = ensg_number + 1

    assert num_lines == 13
    assert protein_number == 2
    assert ensg_number == 11


@mock.patch.object(network_id_mapper, 'parse_args')
def test_case_two(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    args.return_value = [example_files_list[0], 'true', 'protein', 'uniprotac', example_files_list[1], '', output_file]
    network_id_mapper.main()

    ensg_number = 0
    with open(output_file) as output:
        for line in output:
            ensg = line.find('ensg')
            if ensg != -1:
                ensg_number = ensg_number + 1

    assert ensg_number == 14


@mock.patch.object(network_id_mapper, 'parse_args')
def test_case_three(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    args.return_value = [example_files_list[0], 'false', 'protein', 'uniprotac', example_files_list[1], '', output_file]
    network_id_mapper.main()

    ensg_number = 0
    with open(output_file) as output:
        for line in output:
            ensg = line.find('ensg')
            if ensg != -1:
                ensg_number = ensg_number + 1

    assert ensg_number == 16


@mock.patch.object(network_id_mapper, 'parse_args')
def test_case_four(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    args.return_value = [example_files_list[0], 'true', '', 'uniprotac', example_files_list[1], '', output_file]
    network_id_mapper.main()

    p11007_number = 0
    with open(output_file) as output:
        for line in output:
            p11007 = line.find('p11007')
            if p11007 != -1:
                p11007_number = p11007_number + 1

    assert p11007_number == 4


@mock.patch.object(network_id_mapper, 'parse_args')
def test_case_five(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    args.return_value = [example_files_list[0], 'false', '', 'test', example_files_list[1], '', output_file]
    network_id_mapper.main()

    num_lines_input_file = sum(1 for line in open(example_files_list[0]))
    num_lines_output_file = sum(1 for line in open(output_file))

    assert num_lines_input_file == num_lines_output_file


@mock.patch.object(network_id_mapper, 'parse_args')
def test_case_six(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    args.return_value = [example_files_list[0], 'true', '', 'test', example_files_list[1], '', output_file]
    network_id_mapper.main()

    assert os.stat(output_file).st_size == 0


@mock.patch.object(network_id_mapper, 'parse_args')
def test_case_if_uniquename_is_given(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    args.return_value = [example_files_list[2], 'false', '', 'uniprotac', example_files_list[3], 'true', output_file]
    network_id_mapper.main()

    outputfile_rows_int_a = ["uniprotac:p52879", "uniprotac:p29999", "uniprotac:g29999"]
    outputfile_rows_int_b = ["uniprotac:q29999", "uniprotac:p52879", "uniprotac:p52879"]

    result_a = []
    result_b = []

    with open(output_file) as output:
        for line in output:
            line = line.strip().split("\t")
            int_a = line[0]
            int_b = line[1]
            if int_a in outputfile_rows_int_a:
                result_a.append(int_a)
            if int_b in outputfile_rows_int_b:
                result_b.append(int_b)

    index = 0
    for member in outputfile_rows_int_a:
        assert member == result_a[index]
        index = index + 1

    index2 = 0
    for member in outputfile_rows_int_b:
        assert member == result_b[index2]
        index2 = index2 + 1
