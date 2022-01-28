import pytest
import network_combiner
import mock
import os


params_list = [
    "example_files/test1.tsv,example_files/test2.tsv,example_files/test3.tsv,example_files/empty.tsv",
    "example_files/output.tsv",
    "union"
]

not_empty_input_files = "example_files/test1.tsv,example_files/test2.tsv"


@mock.patch.object(network_combiner, 'parse_args')
def test_vcf_files_exists(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = ['file', False, 'asdasd']
        network_combiner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


@mock.patch.object(network_combiner, 'parse_args')
def test_method_exists(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [params_list[0], False, 'asdasd']
        network_combiner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 2


@mock.patch.object(network_combiner, 'parse_args')
def test_one_emtpy_input_file(args, tmpdir):
    fake_file = tmpdir.join("wrong_file.tsv")
    fake_file.write("")
    args.return_value = [str(fake_file), params_list[1], params_list[2]]
    network_combiner.main()

    assert os.stat(fake_file).st_size == 0


@mock.patch.object(network_combiner, 'parse_args')
def test_number_of_lines_in_output_file(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    output_file.write("")
    args.return_value = [not_empty_input_files, output_file, params_list[2]]
    network_combiner.main()

    num_lines = sum(1 for line in open(output_file))

    if params_list[2] == "union":
        assert num_lines == 5

    elif params_list[2] == "intersection":
        assert num_lines == 2

    elif params_list[2] == "difference":
        assert num_lines == 3


@mock.patch.object(network_combiner, 'parse_args')
def test_if_one_of_the_input_file_is_empty(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    output_file.write("")
    args.return_value = [params_list[0], output_file, params_list[2]]
    network_combiner.main()

    num_lines = sum(1 for line in open(output_file))

    if params_list[2] == "union":
        assert num_lines == 5

    elif params_list[2] == "intersection":
        assert num_lines == 0

    elif params_list[2] == "difference":
        assert num_lines == 5


@mock.patch.object(network_combiner, 'parse_args')
def test_snp_metadata(args, tmpdir):
    output_file = tmpdir.join("output.tsv")
    output_file.write("")
    args.return_value = [not_empty_input_files, output_file, params_list[2]]
    network_combiner.main()

    num_lines = sum(1 for line in open(output_file))

    assert num_lines == 5

    snp_metadata_array = ["origin:snp;dbsnp;rs00004", "origin:snp;dbsnp;rs00006", "origin:snp;dbsnp;rs00001",
                          "origin:snp;dbsnp;rs00007|origin:snp;dbsnp;rs00009", "origin:snp;dbsnp;rs00017"]

    index = 0
    with open(output_file, 'r') as f:

        for line in f:
            line = line.strip().split('\t')
            snp_metadata = line[27]

            assert snp_metadata == snp_metadata_array[index]

            index = index + 1

