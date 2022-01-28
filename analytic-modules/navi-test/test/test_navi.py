""" A simple Test for navi_test, Author: Madgwick """

import time

import pytest
from navi_test import NaviTest


@pytest.fixture(scope='session')
def navi():
    """ Create a testing instance for Navi """
    # Testing parameters
    input_file = "input_test.csv"
    exit_code = 42
    new_navi = NaviTest(input_file, exit_code)

    return new_navi


@pytest.fixture(scope='session')
def navi_not_existing_input_file():
    """ Create a testing instance for Navi """
    # Testing parameters
    input_file = "this_file_does_not_exists.csv"
    exit_code = 42
    new_navi = NaviTest(input_file, exit_code)

    return new_navi


@pytest.fixture()
def expected_file_format():
    """ Expected output csv file as a string """
    expected = "structureId,chainId,residueCount,macromoleculeType\n" \
               "100D,A,20,DNA/RNA Hybrid\n" \
               "100D,B,20,DNA/RNA Hybrid"

    return expected


@pytest.mark.parametrize('wait_sec', [1, 2, 3], ids=['wait 1s', 'wait 2s', 'wait 3s'])
def test_wait(navi, wait_sec):
    """ Test the wait time within the wait time list with ids """
    start = time.time()
    navi.wait(wait_sec)
    seconds = int(round(time.time() - start))
    assert seconds == wait_sec


def test_input_file_check(navi_not_existing_input_file):
    """ Test if check_input_file exits for not existing files """
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        navi_not_existing_input_file.check_input_file()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


def test_write_output(tmpdir, expected_file_format, navi):
    """ Test the writing of the file by creating a new temp directory """
    output_file = 'out_test.csv'
    file = tmpdir.join(output_file)
    navi.copy_to_output(str(file))
    print("outfile: " + str(file))
    assert file.read() == expected_file_format


def test_exit_code(navi):
    """ Test the exit code """
    exit_code = 42
    assert exit_code == navi.exit_code


def test_exit(navi):
    """ Testing system exit codes """
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        navi.navi_exit()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 42
