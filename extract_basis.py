# Script from extracting and rewriting basis sets, either from BSSE or using a .json-file
import matplotlib.pyplot as plt
import numpy as np
import basis_set_exchange
import hashlib
import json
import sys

# Element names, numbers, and angular momenta
Z_charge = np.arange(117)+1
Z_label  = ['H','HE','LI','BE','B','C','N','O','F','NE','NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE','I','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','ES','FM','MD','NO','LR','RF','DB','SG','BH','HS','MT','DS','RG','CN','NH','FL','MC','LV','TS','OG']
ang_mom_vect = ['S','P','D','F','G','H','I','J','K']

# Control: print contraction of basis set
def print_gtos(pgtos,cgtos,ang_mom_vect):
    pgto_info,cgto_info = '',''
    for n in np.arange(len(pgtos)):
        if pgtos[n] > 0:
            pgto_info = pgto_info+str(int(pgtos[n]))+ang_mom_vect[n]
            if pgtos[n+1] > 0: pgto_info = pgto_info+','
    for n in np.arange(len(cgtos)):
        if cgtos[n] > 0:
            cgto_info = cgto_info+str(int(cgtos[n]))+ang_mom_vect[n]
            if cgtos[n+1] > 0: cgto_info = cgto_info+','
    print('('+pgto_info.lower()+') -> ['+cgto_info.lower()+']')
    return False

# Basis set name, in bsse or as .json
basis     = '6-31G'
json_file = False

# Load data
if json_file:
    f = open(basis,); data = json.load(f); f.close(); bs_dict = data
else:
    bs_dict = basis_set_exchange.get_basis(basis, header=False)
elements = bs_dict['elements']; element_avail = []; element_label = []
for Z in Z_charge:
    try:
        element = elements['{}'.format(int(Z))]; element_label.append(Z_label[Z-1]); element_avail.append(element)
    except:
        continue
print('Found basis for\n',element_label)

# Rewrite to vlx format
original_stdout = sys.stdout     
with open('basis/{}'.format(basis), 'w') as f:
    sys.stdout = f 
    print('@BASIS_SET',basis)
    Z_indx = 0
    for element in element_avail:
        pgtos,cgtos = np.zeros(len(ang_mom_vect)), np.zeros(len(ang_mom_vect))
        sys.stdout = f 
        basis_info = element['electron_shells']
        print('\n@ATOMBASIS',element_label[Z_indx])
        for i in np.arange(len(basis_info)):
            block = basis_info[i]
            l_mom = block['angular_momentum']
            if len(l_mom) == 1:
                pgtos[l_mom[0]] += len(block['exponents'])
                ang_mom = ang_mom_vect[l_mom[0]]
                for l in np.arange(len(block['coefficients'])):
                    cgtos[l_mom[0]] += 1
                    tmp_print = ''
                    n_prim = 0
                    for k in np.arange(len(block['exponents'])):
                        if np.float(block['coefficients'][l][k]) != 0.0:
                            expo = '{:.12e}'.format(np.float(block['exponents'][k].replace('E','e')))
                            coef = '{:.12e}'.format(np.float(block['coefficients'][l][k].replace('E','e')))
                            if np.abs(np.float(coef)) > 0.0:
                                if n_prim != 0:
                                    tmp_print += '\n'
                                if np.float(coef) >= 0.0: tmp_print += expo+'  '+coef
                                else: tmp_print += expo+' '+coef
                                n_prim += 1
                    print(ang_mom,n_prim,'  1')
                    print(tmp_print)
            else:
                for j in np.arange(len(l_mom)):
                    cgtos[l_mom[j]] += 1
                    pgtos[l_mom[j]] += len(block['coefficients'])
                    ang_mom = ang_mom_vect[l_mom[j]]
                    tmp_print = ''
                    n_prim = 0
                    for k in np.arange(len(block['exponents'])):
                        if np.float(block['coefficients'][l][k]) != 0.0:
                            expo = '{:.12e}'.format(np.float(block['exponents'][k].replace('E','e')))
                            coef = '{:.12e}'.format(np.float(block['coefficients'][l+j][k].replace('E','e')))
                            if np.abs(np.float(coef)) > 0.0:
                                if n_prim != 0:
                                    tmp_print += '\n'
                                if np.float(coef) >= 0.0: tmp_print += expo+'  '+coef
                                else: tmp_print += expo+' '+coef
                                n_prim += 1
                    print(ang_mom,n_prim,'  1')
                    print(tmp_print)
        print('@END')
        sys.stdout = original_stdout 
        print(element_label[Z_indx],'contraction:')
        print_gtos(pgtos,cgtos,ang_mom_vect)
        Z_indx += 1
    sys.stdout = original_stdout 
    f.close()

# Add md5-sum
f = open('basis/{}'.format(basis), 'r')
data = ''
for x in f:
    data += x
f.close() 
md5 = hashlib.md5(data.encode('utf-8')).hexdigest()
f = open('basis/{}'.format(basis), 'a')
f.write(md5)
f.close()
