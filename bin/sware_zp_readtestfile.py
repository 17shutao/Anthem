###read test peptide input 
import sys, os
from collections import OrderedDict

def readtestfile(length, test_peptide):
    if test_peptide != None:
        if os.path.exists(test_peptide) == False:
            print('Error: "' + test_peptide + '" does not exist.')
            sys.exit(1)    
        
        upscaleAA = ['A', 'R', 'N', 'D', 'V', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'C']
        unusualAA = ['B', 'J', 'O', 'U', 'X', 'Z', 'b', 'j', 'o', 'u', 'x', 'z']
        peptides, labels = [], []
    
        with open(test_peptide, 'rt') as f:
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
                    print("Please put the lable '1' at the first row.")
                    sys.exit(1)
                rownum -= 1     
            elif '-1' in eachline:
                negstart += 1                
            elif posstart == 1 and negstart == 0:
                peptide = eachline.split('\n')[0]
                ###check unrecognised character
                for eachAA in peptide:
                    if eachAA not in upscaleAA:
                        print('Unrecognised character exists in your training file')###check characters
                        sys.exit(1)         
                ###check peptide length    
                if len(peptide) == int(length):
                    peptides.append(peptide)
                    labels.append('1')
                    numpos += 1
                elif len(peptide) == 0:
                    continue
                else:
                    print('The length of test peptide is not same with the length you choose')
                    sys.exit(1)                
            else:
                peptide = eachline.split('\n')[0]
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
                    print('The length of test peptide is not same with the length you choose')
                    sys.exit(1)                
    
        postestpeptide, negtestpeptide = [], []
        for items in enumerate(peptides):
            if labels[items[0]] == '1':
                postestpeptide.append(items[1])
            elif labels[items[0]] == '-1':
                negtestpeptide.append(items[1])
                            
        dict_seqlabel = OrderedDict()
        for eachpos in postestpeptide:
            dict_seqlabel[eachpos] = '1'        
        for eachneg in negtestpeptide:
            dict_seqlabel[eachneg] = '-1'                
    
        if '-1' not in labels or '1' not in labels:
            print('The test file must contain both positive and negative labeled peptide')
            sys.exit(1)
        else:
            return(dict_seqlabel)
    else:
        return(None)