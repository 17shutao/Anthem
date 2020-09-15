###read predict peptide input 
import sys, os, json
from collections import OrderedDict

def readpredictfile(length, predictfile_peptide, allmatrix = None):
    if os.path.exists(predictfile_peptide) == False:
        print('Error: "' + predictfile_peptide + '" does not exist.')
        sys.exit(1)    
    upscaleAA = ['A', 'R', 'N', 'D', 'V', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'C']
    unusualAA = ['B', 'J', 'O', 'U', 'X', 'Z', 'b', 'j', 'o', 'u', 'x', 'z']
    peptides = []

    with open(predictfile_peptide, 'rt') as f:
        content = f.readlines()
        f.seek(0, 0)
        content_txt = f.read()
        
    ###check unusual amino acid
    dict_nameseqs = OrderedDict()
    for eachAA in unusualAA:
        if eachAA in content_txt:
            print('Peptides contain unnatural amino acid, eg: B, J, O, U, X, Z')
            sys.exit(1)

    for eachline in content:
        eachline = (eachline.strip()).upper()
        for eachAA in eachline:
            if eachAA not in upscaleAA:
                print('Unrecognised character exists in your peptide file')###check characters
                sys.exit(1) 
        if len(eachline) == int(length):
            peptides.append(eachline)
        elif len(eachline) == 0:
            continue
        else:
            print('The peptide length is not same with the length you choose')
            sys.exit(1)
    dict_nameseqs['input-sequences'] = peptides
                    
    if allmatrix == None:
        return(dict_nameseqs)
        
    elif allmatrix != None:       
        for HLA, matrixs in allmatrix.items():
            HLAname = HLA
            for eachmatrixname, eachmatrix in matrixs.items():
                if 'AAF' in eachmatrixname:
                    for values in eachmatrix.values():
                        matrixlen = len(values)
                        break                        
        for peptides in dict_nameseqs.values():
            for eachpep in peptides:
                if len(eachpep) != matrixlen:
                    print('Peptide in file should be the same with the length that model can predict')
                    sys.exit(1)    
    
        return(dict_nameseqs, HLAname)