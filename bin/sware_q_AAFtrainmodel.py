###generate AAF for train peptide
from collections import OrderedDict
import json

def AAF(predictpeptide, length, trainpeptide, testpeptide, matrix):
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
                    dict_AAfre[eachseq[i]] += 1
            for eachAA, eachvallist in dict_AAlenfre.items():
                eachvallist.append((dict_AAfre[eachAA]/len(postrainpep))#each AA percentage
                                   /
                                   (max(dict_AAfre.values())/len(postrainpep)))#AA percentage which is largest    
                
        dict_trainAAF = OrderedDict()
        for eachseq in postrainpep:
            scoreseq = 0
            for i in range(int(length)):
                if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                    scoreseq += 0
                else:
                    scoreseq += dict_AAlenfre[eachseq[i]][i]
            dict_trainAAF.setdefault('pos', []).append((eachseq, scoreseq))
            
        dict_trainAAF['neg'] = OrderedDict()
        for eachset, negpeps in dict_setnegtrainpep.items():
            for eachseq in negpeps:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAlenfre[eachseq[i]][i]
                dict_trainAAF['neg'].setdefault(eachset, []).append((eachseq, scoreseq))            
            
        AAFmatrix = {}
        AAFmatrix['AAF'] = dict_AAlenfre                        
        
        if predictpeptide == None:
            dict_predictnameseqAAF = None
        else:            
            dict_predictnameseqAAF = OrderedDict()
            for eachname, seqs in predictpeptide.items():
                dict_predictnameseqAAF[eachname] = OrderedDict()
                for eachseq in seqs:
                    scoreseq = 0
                    for i in range(int(length)):
                        if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                            scoreseq += 0
                        else:
                            scoreseq += dict_AAlenfre[eachseq[i]][i]
                    dict_predictnameseqAAF[eachname][eachseq] = scoreseq
        
        if testpeptide == None:
            dict_testAAF = None
        else:
            dict_testAAF = OrderedDict()
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
                dict_testAAF[eachseq] = (scoreseq, 1)            
            for eachseq in negtestpep:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAlenfre[eachseq[i]][i]
                dict_testAAF[eachseq] = (scoreseq, -1)            
                                            
        return(dict_trainAAF, dict_predictnameseqAAF, dict_testAAF, AAFmatrix)
    
    else:        
        for HLA, matrixs in matrix.items():
            dict_AAlenfre = matrixs['AAF']
        
        lengthmatrix = len(dict_AAlenfre['A']) 
        
        dict_predictnameseqAAF = OrderedDict()
        for eachname, seqs in predictpeptide.items():
            dict_predictnameseqAAF[eachname] = OrderedDict()
            for eachseq in seqs:
                scoreseq = 0
                for i in range(lengthmatrix):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAlenfre[eachseq[i]][i]
                dict_predictnameseqAAF[eachname][eachseq] = scoreseq 
                
        return(dict_predictnameseqAAF)    