import os
import sys
import mock
import pytest
import vcf_iterator

sys.path.append("./analytic-modules/")


params_list = [
    "example_files/vcf_data/margaret.vcf,example_files/vcf_data/norman.vcf,example_files/vcf_data/bucky.vcf",
    "example_files/patients.tar.gz",
    "example_files/example_bim_fam/norwich.ibd.ichip.bed,example_files/example_bim_fam/norwich.ibd.ichip.bim,example_files/example_bim_fam/norwich.ibd.ichip.fam",
    "example_files/output/"
]


@mock.patch.object(vcf_iterator, 'parse_args')
def test_vcf_files_exists(args):
    """ Exception handling test: no vcf file found """
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = ['file', False, False, params_list[3]]
        vcf_iterator.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


@mock.patch.object(vcf_iterator, 'parse_args')
def test_tar_file_exists(args):
    """ Exception handling test: no tar file found """
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [params_list[0], 'no_tar', False, params_list[3]]
        vcf_iterator.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 2


@mock.patch.object(vcf_iterator, 'parse_args')
def test_plink_file_exists(args):
    """ Exception handling test: no plink file found """
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = [params_list[0], False, 'no_plink', params_list[3]]
        vcf_iterator.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 3


@mock.patch.object(vcf_iterator, 'parse_args')
def test_no_output_folder(args):
    """ Exception handling test: no output folder found"""
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        args.return_value = ['file', False, False, 'output']
        vcf_iterator.main()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 4


@mock.patch.object(vcf_iterator, 'parse_args')
def test_component_vcf_iterator(args, tmpdir):
    """ Exception handling test: component test for just vcf files """
    output_folder = tmpdir.mkdir('output_test')
    args.return_value = [params_list[0], False, False, output_folder]
    vcf_iterator.main()

    assert len([name for name in os.listdir(output_folder)]) == 3


def test_plink_extraction(tmpdir):
    output_folder = tmpdir.mkdir('output_test')
    fam_file = params_list[2].split(",")[2]
    num_lines = sum(1 for line in open(fam_file))

    vcf_iterator.plink_extraction(params_list[2], output_folder, index_plink=1)
    num_files = len([name for name in os.listdir(output_folder)])

    assert num_lines == num_files


def test_plink_output_vcf_files_exists(tmpdir):
    output_folder = tmpdir.mkdir('output_test')
    vcf_iterator.plink_extraction(params_list[2], output_folder, index_plink=1)
    for file in os.listdir(output_folder):

        assert file.endswith(".vcf")

