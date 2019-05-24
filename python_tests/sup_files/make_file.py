from veloxchem.inputparser import InputParser

def make_file(filename,fout):
            basis_parser = InputParser(filename)
            basis_dict = basis_parser.get_dict()
            L = list(basis_dict.keys())
            L.pop(-1)
            fo = open(fout, "w")
            fo.write("@BASIS_SET " + str(filename.split("/")[2].upper()) + "\n\n")
            #print(str(filename.split("/")[2].upper()))
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
                    #fo.write(str(basis_dict[l][i])+"\n")
                    #print(str(basis_dict[l][i]).split())
                    if basis_dict[l][i].split()[0] in "SPDFGHIKL":
                            #print(basis_dict[l][i])
                            fo.write(str(basis_dict[l][i]))
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

if __name__ == "__main__":
        make_file(filename,fout)
