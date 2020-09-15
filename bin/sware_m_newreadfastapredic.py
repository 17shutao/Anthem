###read predict peptide with fasta format
import os, sys, json
from collections import OrderedDict

def readpredictfasta(length, trainfile_fasta, allmatrix = None):
    if os.path.exists(trainfile_fasta) == False:
        print('Error: "' + trainfile_fasta + '" does not exist.')
        sys.exit(1)    
    upscaleAA = ['A', 'R', 'N', 'D', 'V', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'C']                
    unusualAA = ['B', 'J', 'O', 'U', 'X', 'Z', 'b', 'j', 'o', 'u', 'x', 'z']      
        
    with open(trainfile_fasta) as f:
        content = f.read()

    ###check fasta format
    if '>' not in content:
        print('The input file seems not in fasta format.')
        sys.exit(1)
        
    records = content.split('>')[1:] 
    dict_fasta = OrderedDict()    
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
        ###check unusual amino acid
        for eachAA in unusualAA:
            if eachAA in thissequence:
                print('Sequence contain unnatural amino acid, eg: B, J, O, U, X, Z')
                sys.exit(1)                          
        for eachAA in thissequence:
            if eachAA not in upscaleAA:
                print('Unrecognised character exists in your peptide file')###check characters
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
    
    for items in enumerate(seqname):
        dict_fasta[items[1]] = sequence[items[0]]

    if allmatrix == None:
        return(dict_nameseqs, dict_fasta)
    
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
                    print('The length you choose is not the same with the length that model can predict')
                    sys.exit(1)           
        
        return(dict_nameseqs, HLAname, dict_fasta)