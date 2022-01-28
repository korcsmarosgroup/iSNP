import os
import pytest
import mutated_sequence_generator as mut

input_vcf = "./test/sample_shorter.vcf"
genome = "./test/human-annotation_renamed.fasta"
region_length = 5
output_wild_type = "./test/out/wild.fasta"
output_mutated = "./test/out/mutated.fasta"


def test_no_fasta():
    genome = "./test/wrong.fasta"
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 202


def test_no_vcf():
    input_vcf = "./test/wrong.vcf"
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 201


def test_wrong_type_vcf():
    input_vcf = 2
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 101


def test_wrong_type_genome():
    genome = 2
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 102


def test_wrong_type_length():
    region_length = "10"
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 105


def test_wrong_type_output_wild():
    output_wild_type = 2

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 103


def test_wrong_type_output_mutated():
    output_mutated = 4
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 104


def test_negative_length():
    region_length = -1

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 106


def test_mutate_sequence():
    mut.mutate(input_vcf, genome, output_wild_type, output_mutated, region_length)
    with open(output_wild_type) as fin:
        first_line = fin.readline()
        second_line = fin.readline()
        assert first_line == ">origin:snp1;AC=2;AN=2;DB;DP=182;H2;NS=65 | mutated:False\n"
        assert second_line == "AAGACTGAAAT\n"

    with open(output_mutated) as fin:
        first_line = fin.readline()
        second_line = fin.readline()
        assert first_line == ">origin:snp1;AC=2;AN=2;DB;DP=182;H2;NS=65 | mutated:True\n"
        assert second_line == "AAGACCGAAAT\n"







