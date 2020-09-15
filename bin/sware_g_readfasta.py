###read the HLA and sequence with fasta format
import os, sys, re
from collections import OrderedDict
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_i_checkHLA, sware_t_checkHLAformat

def readfasta(length, HLAtype, fasta_file):
    if os.path.exists(fasta_file) == False:
        print('Error: "' + fasta_file + '" does not exist.')
        sys.exit(1)    
    upscaleAA = ['A', 'R', 'N', 'D', 'V', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'C']                
    unusualAA = ['B', 'J', 'O', 'U', 'X', 'Z', 'b', 'j', 'o', 'u', 'x', 'z']
        
    with open(fasta_file) as f:
        content = f.read()
    
    ###check fasta format
    if '>' not in content:
        print('The input file seems not in fasta format.')
        sys.exit(1)
            
    records = content.split('>')[1:]
    dict_HLAfasta = OrderedDict()
    seqname, sequence = [], []
    for fasta in records:
        array = fasta.split('\n')
        while '|' in array[0]:
            array[0] = array[0].replace('|', '-')    
        while '_' in array[0]:
            array[0] = array[0].replace('_', '-')                     
        thissequence = ''.join(array[1:])
        try:
            seqname.append(array[0].split()[0])
        except:
            print('Something wrong in fasta file, please check. Probably check the name of each sequence.')
            sys.exit(1)
        thissequence = (thissequence.strip()).upper()
        for eachAA in thissequence:
            if eachAA not in upscaleAA:
                print('Unrecognised character exists in your peptide file')###check characters
                sys.exit(1)
        ###check unusuall amino acid
        for eachAA in unusualAA:
            if eachAA in thissequence:
                print('Sequence contain unnatural amino acid, eg: B, J, O, U, X, Z')
                sys.exit(1)                          
        if len(thissequence) < int(length):
            print('Sequence is shorter than the predicted length chosen')
            sys.exit(1)                         
        sequence.append(thissequence)            
            
    dict_nameseqs = OrderedDict()
    for items in enumerate(sequence):
        startpos, stoppos = 0, int(length)
        seqlist = []
        while stoppos < len(items[1]) +1 :
            seqlist.append(items[1][startpos: stoppos])
            startpos += 1
            stoppos += 1
        dict_nameseqs[seqname[items[0]]] = seqlist
    
    ###check input HLA format
    HLAlist = sware_t_checkHLAformat.checkHLAinput(HLAtype, None) 
    for eachHLA in HLAlist:
        dict_HLAfasta[eachHLA] = OrderedDict()
        for items in enumerate(seqname):
            dict_HLAfasta[eachHLA][items[1]] = sequence[items[0]]

    ###check available HLA in the length
    sware_i_checkHLA.checkHLA(length, HLAlist)  
            
    return(HLAlist, dict_nameseqs, dict_HLAfasta)