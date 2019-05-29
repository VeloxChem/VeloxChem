import pathlib
import sys
import pytest
import os
import hashlib, sys
from veloxchem.inputparser import InputParser

cases = [f.name for f in pathlib.Path("../basis").iterdir() if f.is_file]
cases.remove('guess')
cases.remove('MIN-CC-PVDZ')
@pytest.mark.parametrize('case', cases)
def test_acceptance(case):
        this_basis = "../basis/" + case
        make_file(this_basis,"mm")
        # calculae md5 for mm
        calculated_md5 = mdd5('mm')
        # get reference md5 string
        expected_md5 = open(this_basis).readlines()[-1]
        if "\n" in expected_md5:
            expected_md5 = expected_md5.split("\n")[0]
        os.unlink("mm")
        assert expected_md5 == calculated_md5

def make_file(filename,fout):
            basis_parser = InputParser(filename)
            basis_dict = basis_parser.get_dict()
            L = list(basis_dict.keys())
            L.pop(-1)
            fo = open(fout, "w")
            fo.write("@BASIS_SET " + str(filename.split("/")[2].upper()) + "\n\n")
            count = int(0)
            for l in L:
                fo.write(
                    "@"
                    + str(l.split("_")[0].upper())
                    + " "
                    + str(l.split("_")[1].upper())
                    + "\n"
                )
                for i in range((len(basis_dict[str(l)]))):
                    if basis_dict[l][i].split()[0] in "SPDFGHIKL":
                            sentence = str(basis_dict[l][i]).split()[0] + ' '+ str(basis_dict[l][i]).split()[1] + '  ' + str(basis_dict[l][i]).split()[2]
                            fo.write(str(sentence))
                    for z in range(len(basis_dict[l][i].split())):
                            if basis_dict[l][i].split()[0] not in "SPDFGHIKL":
                                    if z == 0:
                                            fo.write("%18s" % basis_dict[l][i].split()[z])
                                    else:
                                            fo.write(" %19s" % basis_dict[l][i].split()[z])
                    fo.write("\n")
                count = count + 1
                if count < len(L):
                    fo.write("@END\n\n")
                else:
                    fo.write("@END\n")        

def mdd5(filename):
    hasher = hashlib.md5()
    f = open(filename, "rb")
    content = f.read()
    hasher.update(content)
    #print(hasher.hexdigest())
    f.close()
    return hasher.hexdigest()
