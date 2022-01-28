""" Automated tests for mirna interaciton prediction """
import time
import pytest
import mirna_interaction_predictor as mirna
from common_libs.mitab_handler import mitab_handler


def test_fasta_parser(tmpdir):
    """ Check that invalid sequence alphabets are not read """
    example_fasta = "example/0001-patient-msg-protein-coding-mutant.fasta"
    clean_file, seq_info = mirna.parse_sequences(example_fasta)
    acutal_file = open(clean_file).readlines()

    assert len(acutal_file) == 6
    assert "X" not in acutal_file


@pytest.mark.parametrize('params', [['a', 'b'], [100, -100]])
def test_paramaters(params):
    """ A test to ensure only valid parameters can be passed """

    energy, score = params

    with pytest.raises(mirna.InvalidMirandaParameter) as pytest_wrapped_e:
        mirna._check_threshold(score, energy)

    assert pytest_wrapped_e.type == mirna.InvalidMirandaParameter


@pytest.mark.parametrize('test_file', ['example/mature-example.fasta'])
def test_mirbase_parser(test_file):
    """ Test the reading of the database file but a standard fasta and a database embl .dat file """
    db, details = mirna.parse_database(test_file)

    with open(db) as mirbase:
        first_line = mirbase.readline()

    assert db.endswith('.fasta')
    assert first_line.__contains__('>')


@pytest.mark.parametrize('test_file', ['random_file.fasta', 'not_file.txt'])
def test_incorrect_file(test_file):
    """ Test invalid file format and none existant files """
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        mirna.check_fasta(test_file)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


def test_create_network_file(tmpdir):
    """ Testing the creation of the network file """

    sequences = "example/0001-patient-msg-protein-coding-mutant.fasta"
    database = "example/mature-example.fasta"
    output = tmpdir.mkdir("sub").join("output_file.inter")

    score = 1
    energy = -1

    database, database_info = mirna.parse_database(database)
    genomic, sequences_info = mirna.parse_sequences(sequences)

    predictions = mirna._predictor(sequences=genomic,
                                   database=database,
                                   score=score,
                                   energy=energy,
                                   strict=False)

    start = time.time()
    mirna_preds = mirna.extract_results(predictions)
    end = time.time()
    print(f"mirand time: {end - start}")

    start = time.time()
    mirna.create_network_file(mirna_preds, sequences_info, output)
    end = time.time()
    print(f"combining time: {end-start}")

    new_h = mitab_handler.MiTabHandler()
    new_h.parse(output, file_format='mitab')

    with open(output, 'r') as test_file:
        check_line = test_file.readline()

    assert check_line.islower()
    assert check_line[62:88] == "taxid:9606('homo sapiens')"


def test_miranda_predictor():
    """ A test to check the miranda sub-process call """

    # Miranda default parameters for testing
    sequences = "example/0001-patient-msg-protein-coding-mutant.fasta"
    database = "example/mature-example.fasta"
    score = 1
    energy = -1

    database, database_info = mirna.parse_database(database)
    genomic, sequences_info = mirna.parse_sequences(sequences)

    predictions = mirna._predictor(sequences=genomic,
                                   database=database,
                                   score=score,
                                   energy=energy,
                                   strict=False)

    mirna_preds = mirna.extract_results(predictions)

    # Get scores and energies
    scores = []
    energies = []

    for preds in mirna_preds:
        scores.append(float(preds.Max_Score))
        energies.append(float(preds.Max_Energy))

    assert len(mirna_preds) == 12
    assert all(s > score for s in scores)
    assert any(e < energy for e in energies)
