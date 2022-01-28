import pytest
import vcf_cleaner
import mock
import os


example_files_list = ["example_files/test.vcf",
                      "example_files/",
                      "example_files/new.txt",
                      "example_files/references/reference_genome_file.bed",
                      "example_files/references/reference_vcf_file.vcf",
                      "example_files/output/"]


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_input_vcf_file_exists(args, tmpdir):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = ['file', False, False, example_files_list[3], False, tmpdir]
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_input_vcf_file_is_not_a_vcf_file(args, tmpdir):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        txt_file = "example_files/new.txt"
        args.return_value = [txt_file, False, False, example_files_list[3], False, tmpdir]
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 2


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_input_folder_exists(args, tmpdir):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [False, 'input_folder', False, example_files_list[3], False, tmpdir]
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 3


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_input_snp_list_exists(args, tmpdir):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [False, False, 'snp_list.txt', False, example_files_list[4], tmpdir]
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 4


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_reference_genome_file_exists(args, tmpdir):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [example_files_list[0], False, False, 'reference.bed', False, tmpdir]
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 5


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_reference_vcf_file_exists(args, tmpdir):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [False, False, example_files_list[2], False, 'reference.vcf', tmpdir]
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 6


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_output_folder_exists(args):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [example_files_list[0], False, False, example_files_list[3], False, 'output']
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 7


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_if_input_vcf_file_and_input_folder_at_the_same_time(args, tmpdir):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [example_files_list[0], example_files_list[1], False, example_files_list[3], False, tmpdir]
        vcf_cleaner.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 8


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_case_one(args, tmpdir):
    output_folder = tmpdir.mkdir('output')
    args.return_value = [example_files_list[0], False, False, example_files_list[3], False, output_folder]
    vcf_cleaner.main()

    print(example_files_list[0])
    print(example_files_list[3])

    assert len([name for name in os.listdir(output_folder)]) == 1

    for filename in os.listdir(output_folder):
        file = os.path.join(output_folder, filename)
        print(file)
        num_lines = sum(1 for line in open(file))
        assert num_lines == 31


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_case_two(args, tmpdir):
    output_folder = tmpdir.mkdir('output')
    args.return_value = [False, example_files_list[1], False, example_files_list[3], False, output_folder]
    vcf_cleaner.main()

    assert len([name for name in os.listdir(output_folder)]) == 2

    for filename in os.listdir(output_folder):
        file = os.path.join(output_folder, filename)
        num_lines = sum(1 for line in open(file))
        assert num_lines == 31


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_case_three(args, tmpdir):
    output_folder = tmpdir.mkdir('output')
    args.return_value = [example_files_list[0], False, False, example_files_list[3], False, output_folder]
    vcf_cleaner.main()

    for filename in os.listdir(output_folder):
        assert len([name for name in os.listdir(output_folder)]) == 1
        file = os.path.join(output_folder, filename)
        locations = []
        with open(file, "r") as output_file:
            for line in output_file:
                line = line.strip()
                if line[:1] != "#":
                    info = line.split("\t")
                    location = info[1]
                    locations.append(location)
        assert locations == ['40972192', '41019095', '41132429']


@mock.patch.object(vcf_cleaner, 'parse_args')
def test_case_four(args, tmpdir):
    output_folder = tmpdir.mkdir('output')
    args.return_value = [False, False, example_files_list[2], False, example_files_list[4], output_folder]
    vcf_cleaner.main()

    assert len([name for name in os.listdir(output_folder)]) == 1

    for filename in os.listdir(output_folder):
        file = os.path.join(output_folder, filename)
        num_lines = sum(1 for line in open(file))
        assert num_lines == 3

