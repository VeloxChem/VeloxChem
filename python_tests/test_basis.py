import pathlib
import sys
import pytest
import os
sys.path.insert(0, 'sup_files')
import make_file, md5_p

cases = [f.name for f in pathlib.Path("../basis").iterdir() if f.is_file]
cases.remove('guess')
cases.remove('MIN-CC-PVDZ')
@pytest.mark.parametrize('case', cases)
def test_acceptance(case):
        this_basis = "../basis/" + case
        make_file.make_file(this_basis,"mm")
        # calculae md5 for mm
        calculated_md5 = md5_p.mdd5('mm')
        # get reference md5 string
        expected_md5 = open(this_basis).readlines()[-1]
        os.unlink("mm")
        assert expected_md5 == calculated_md5
        

