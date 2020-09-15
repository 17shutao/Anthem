###read training peptide input 
import sys, os, re
from collections import OrderedDict

def readtrainfile(length, trainfile_peptide):
    if os.path.exists(trainfile_peptide) == False:
        print('Error: "' + trainfile_peptide + '" does not exist.')
        sys.exit(1)    
    
    upscaleAA = ['A', 'R', 'N', 'D', 'V', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'C']
    unusualAA = ['B', 'J', 'O', 'U', 'X', 'Z', 'b', 'j', 'o', 'u', 'x', 'z']
    peptides, labels = [], []

    with open(trainfile_peptide, 'rt') as f:
        content = f.readlines()
        f.seek(0, 0)
        content_txt = f.read()

    ###check unusual amino acid
    for eachAA in unusualAA:
        if eachAA in content_txt:
            print('Peptides contain unnatural amino acid, eg: B, J, O, U, X, Z')
            sys.exit(1)
    
    rownum, posstart, numpos, negstart = 1, 0, 0, 0
    for eachline in content:
        ###check pos label at first row
        if rownum == 1:
            if '1' in eachline:
                posstart += 1
            else:
                print("Please put the label '1' at the first row.")
                sys.exit(1)
            rownum -= 1
        elif '-1' in eachline:
            negstart += 1
        elif posstart == 1 and negstart == 0:
            peptide = (eachline.split('\n')[0]).upper()
            ###check unrecognised character
            for eachAA in peptide:
                if eachAA not in upscaleAA:
                    print('Unrecognized character exists in your training file')###check characters
                    sys.exit(1)
            ###check peptide length
            if len(peptide) == int(length):
                peptides.append(peptide)
                labels.append('1')
                numpos += 1
            elif len(peptide) == 0:
                continue
            else:
                print('The length of training peptide is not same with the length you choose')
                sys.exit(1)
        else:
            peptide = (eachline.split('\n')[0]).upper()
            ###check unrecognised character
            for eachAA in peptide:
                if eachAA not in upscaleAA:
                    print('Unrecognised character exists in your training file')###check characters
                    sys.exit(1)         
            ###check peptide length    
            if len(peptide) == int(length):
                peptides.append(peptide)
                labels.append('-1')
            elif len(peptide) == 0:
                continue
            else:
                print('The length of training peptide is not same with the length you choose')
                sys.exit(1)                       
                          
    ###check number of training data
    if len(set(peptides)) < 5:
        print('The number of peptide in training set should be >= 5 unique instances')
        sys.exit(1)

    dict_labelpep = OrderedDict()
    for items in enumerate(peptides):
        dict_labelpep[items[1]] = labels[items[0]]

    if '-1' not in labels:
        return(dict_labelpep, 'NO')
    else:
        dict_possetneg = OrderedDict()
        dict_possetneg['neg'] = OrderedDict()
        for labelitems in enumerate(labels):
            if labelitems[1] == '1':
                dict_possetneg.setdefault('pos', []).append(peptides[labelitems[0]])
            elif labelitems[1] == '-1':
                dict_possetneg['neg'].setdefault('0', []).append(peptides[labelitems[0]])                
        return(dict_possetneg, 'YES')