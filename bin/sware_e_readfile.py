###read the HLA and peptide input
import sys, os, re
from collections import OrderedDict
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_i_checkHLA, sware_t_checkHLAformat

def readfile(HLAtype, peptide_file):
    if os.path.exists(peptide_file) == False:
        print('Error: "' + peptide_file + '" does not exist.')
        sys.exit(1)    
    upscaleAA = ['A', 'R', 'N', 'D', 'V', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'C']        
    unusualAA = ['B', 'J', 'O', 'U', 'X', 'Z', 'b', 'j', 'o', 'u', 'x', 'z']
    peptides = []

    with open(peptide_file, 'rt') as f:
        content = f.readlines()
        f.seek(0, 0)
        content_txt = f.read()   
        
    ###check unusual amino acid
    dict_nameseqs = OrderedDict()
    for eachAA in unusualAA:
        if eachAA in content_txt:
            print('Peptides contain unnatural amino acid, eg: B, J, O, U, X, Z')
            sys.exit(1)

    linerow = 0
    for eachline in content:
        eachline = (eachline.strip()).upper()
        if linerow == 0:
            peptidelength = len(eachline)
        linerow += 1
        for eachAA in eachline:
            if eachAA not in upscaleAA:
                print('Unrecognised character exists in your peptide file')###check characters
                sys.exit(1) 
        if len(eachline) == peptidelength:
            peptides.append(eachline)
        elif len(eachline) == 0:
            continue
        else:
            print('The length of peptides in your file is not same')
            sys.exit(1)
    dict_nameseqs['input-sequences'] = peptides
                                
    ###check input HLA format                                
    HLAlist = sware_t_checkHLAformat.checkHLAinput(HLAtype, None)      
    
    ###check available HLA in the length
    sware_i_checkHLA.checkHLA(str(peptidelength), HLAlist)
        
    return(HLAlist, dict_nameseqs, peptidelength)