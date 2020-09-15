###generate BSI62 for train peptide
from collections import OrderedDict
import json

def BSI62(predictpeptide, length, trainpeptide, testpeptide, matrix):
    blosum62 = {
        'A': [4,  -1, -2, -2, 0,  -1, -1, 0, -2,  -1, -1, -1, -1, -2, -1, 1,  0,  -3, -2, 0],  # A
        'R': [-1, 5,  0,  -2, -3, 1,  0,  -2, 0,  -3, -2, 2,  -1, -3, -2, -1, -1, -3, -2, -3], # R
        'N': [-2, 0,  6,  1,  -3, 0,  0,  0,  1,  -3, -3, 0,  -2, -3, -2, 1,  0,  -4, -2, -3], # N
        'D': [-2, -2, 1,  6,  -3, 0,  2,  -1, -1, -3, -4, -1, -3, -3, -1, 0,  -1, -4, -3, -3], # D
        'C': [0,  -3, -3, -3, 9,  -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], # C
        'Q': [-1, 1,  0,  0,  -3, 5,  2,  -2, 0,  -3, -2, 1,  0,  -3, -1, 0,  -1, -2, -1, -2], # Q
        'E': [-1, 0,  0,  2,  -4, 2,  5,  -2, 0,  -3, -3, 1,  -2, -3, -1, 0,  -1, -3, -2, -2], # E
        'G': [0,  -2, 0,  -1, -3, -2, -2, 6,  -2, -4, -4, -2, -3, -3, -2, 0,  -2, -2, -3, -3], # G
        'H': [-2, 0,  1,  -1, -3, 0,  0,  -2, 8,  -3, -3, -1, -2, -1, -2, -1, -2, -2, 2,  -3], # H
        'I': [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4,  2,  -3, 1,  0,  -3, -2, -1, -3, -1, 3],  # I
        'L': [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2,  4,  -2, 2,  0,  -3, -2, -1, -2, -1, 1],  # L
        'K': [-1, 2,  0,  -1, -3, 1,  1,  -2, -1, -3, -2, 5,  -1, -3, -1, 0,  -1, -3, -2, -2], # K
        'M': [-1, -1, -2, -3, -1, 0,  -2, -3, -2, 1,  2,  -1, 5,  0,  -2, -1, -1, -1, -1, 1],  # M
        'F': [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0,  0,  -3, 0,  6,  -4, -2, -2, 1,  3,  -1], # F
        'P': [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7,  -1, -1, -4, -3, -2], # P
        'S': [1,  -1, 1,  0,  -1, 0,  0,  0,  -1, -2, -2, 0,  -1, -2, -1, 4,  1,  -3, -2, -2], # S
        'T': [0,  -1, 0,  -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1,  5,  -2, -2, 0],  # T
        'W': [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1,  -4, -3, -2, 11, 2,  -3], # W
        'Y': [-2, -2, -2, -3, -2, -1, -2, -3, 2,  -1, -1, -2, -1, 3,  -3, -2, -2, 2,  7,  -1], # Y
        'V': [0,  -3, -3, -3, -1, -2, -2, -3, -3, 3,  1,  -2, 1,  -1, -2, -2, 0,  -3, -1, 4],  # V
        '-': [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],  # -
    }    
    AA = 'ARNDCQEGHILKMFPSTWYV-'         
    
    if matrix == None:
        postrainpep, dict_setnegtrainpep = [], OrderedDict()
        for eachpos in trainpeptide['pos']:
            postrainpep.append(eachpos)
        for eachset, negpeps in trainpeptide['neg'].items():
            for eachneg in negpeps:
                dict_setnegtrainpep.setdefault(eachset, []).append(eachneg)
               
        dict_AAlenfre = {
                'A': [], 'R': [], 'N': [], 'D': [], 'C': [],
                'Q': [], 'E': [], 'G': [], 'H': [], 'I': [],
                'L': [], 'K': [], 'M': [], 'F': [], 'P': [],
                'S': [], 'T': [], 'W': [], 'Y': [], 'V': [],
                }
        for i in range(int(length)):
            dict_AAfre = {
                'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0,
                'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
                'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0,
                'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0,
                }
            for eachseq in postrainpep:
                if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                    continue
                else:
                    dict_AAfre[eachseq[i]] += blosum62[eachseq[i]][AA.index(eachseq[i])]
            for eachAA, eachvallist in dict_AAlenfre.items():
                eachvallist.append(dict_AAfre[eachAA]/len(postrainpep))
                               
        dict_trainBSI62 = OrderedDict()
        for eachseq in postrainpep:
            scoreseq = 0
            for i in range(int(length)):
                if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                    scoreseq += 0
                else:
                    scoreseq += dict_AAlenfre[eachseq[i]][i]
            dict_trainBSI62.setdefault('pos', []).append((eachseq, scoreseq))                        
            
        dict_trainBSI62['neg'] = OrderedDict()
        for eachset, negpeps in dict_setnegtrainpep.items():
            for eachseq in negpeps:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAlenfre[eachseq[i]][i]
                dict_trainBSI62['neg'].setdefault(eachset, []).append((eachseq, scoreseq))
            
        BSI62matrix = {}
        BSI62matrix['BSI62'] = dict_AAlenfre

        if predictpeptide == None:
            dict_predictnameseqBSI62 = None
        else:
            dict_predictnameseqBSI62 = OrderedDict()
            for eachname, seqs in predictpeptide.items():
                dict_predictnameseqBSI62[eachname] = OrderedDict()
                for eachseq in seqs:
                    scoreseq = 0
                    for i in range(int(length)):
                        if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                            scoreseq += 0
                        else:
                            scoreseq += dict_AAlenfre[eachseq[i]][i]
                    dict_predictnameseqBSI62[eachname][eachseq] = scoreseq
                
        if testpeptide == None:
            dict_testBSI62 = None
        else:
            dict_testBSI62 = OrderedDict()
            postestpep, negtestpep = [], []
            for key, value in testpeptide.items():
                if value == '1':
                    postestpep.append(key)
                else:
                    negtestpep.append(key)
                    
            for eachseq in postestpep:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAlenfre[eachseq[i]][i]
                dict_testBSI62[eachseq] = (scoreseq, 1)            
            for eachseq in negtestpep:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAlenfre[eachseq[i]][i]
                dict_testBSI62[eachseq] = (scoreseq, -1)                                         
                
        return(dict_trainBSI62, dict_predictnameseqBSI62, dict_testBSI62, BSI62matrix)
    
    else:
        for HLA, matrixs in matrix.items():
            dict_AAlenfre = matrixs['BSI62']
            
        lengthmatrix = len(dict_AAlenfre['A']) 
        
        dict_predictnameseqBSI62 = OrderedDict()
        for eachname, seqs in predictpeptide.items():
            dict_predictnameseqBSI62[eachname] = OrderedDict()
            for eachseq in seqs:
                scoreseq = 0
                for i in range(lengthmatrix):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAlenfre[eachseq[i]][i]
                dict_predictnameseqBSI62[eachname][eachseq] = scoreseq 
                
        return(dict_predictnameseqBSI62)                        