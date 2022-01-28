""" Tests for vcf filtering """
import pytest
from vcf_main_filter import vcf_filter


INPUT_FILE = 'sample.vcf'
FILE_NAME = 'test.vcf'
ANNOTATION = 'sample_nochr.bed'
SNP = 'snp_list.txt'
EXPECTED_FILE_FORMAT = 'example_filter.vcf'


def test_vcf_filtering_component(tmpdir):
    """ Component test for vcf filtering """
    actual_output_file = tmpdir.join(str(FILE_NAME))
    vcf_filter(input=INPUT_FILE,
               output=actual_output_file,
               annotation=ANNOTATION,
               snp=SNP)
    assert open(actual_output_file).read() == open(EXPECTED_FILE_FORMAT).read()


@pytest.mark.parametrize('next_file', ['filenotfound.vcf'])
def test_vcf_filtering_wrong_file(tmpdir, next_file):
    actual_output_file = tmpdir.join(str(FILE_NAME))
    with pytest.raises(FileNotFoundError) as excinfo:
        vcf_filter(input=next_file,
                   output=actual_output_file,
                   annotation=ANNOTATION,
                   snp=SNP)
    assert excinfo.type == FileNotFoundError
