###generate IC50 for train peptide
from collections import OrderedDict
import json, math, os

def IC50(predictpeptide, length, trainpeptide, testpeptide, matrix):
    currentpath = os.getcwd()
    if matrix == None:
        with open('%s/source/AAUniProt.json' % currentpath, 'r') as f:
            dict_AAUniProtAF = json.load(f)
        
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
                if dict_AAfre[eachAA] == 0:
                    eachvallist.append(1/len(postrainpep))#set the minimum value
                else:
                    eachvallist.append(dict_AAfre[eachAA]/len(postrainpep))
                    
        dict_AAPSSM = {
                'A': [], 'R': [], 'N': [], 'D': [], 'C': [],
                'Q': [], 'E': [], 'G': [], 'H': [], 'I': [],
                'L': [], 'K': [], 'M': [], 'F': [], 'P': [],
                'S': [], 'T': [], 'W': [], 'Y': [], 'V': [],
                }            
        for eachAA, eachAF in dict_AAlenfre.items():
            for eachposAF in eachAF:
                powerval = (eachposAF)/dict_AAUniProtAF[eachAA]
                posPSSM = math.log(powerval, 10)
                dict_AAPSSM[eachAA].append(posPSSM)
                
        dict_trainIC50 = OrderedDict()
        for eachseq in postrainpep:    
            PSSMscore = 0
            for i in range(int(length)):
                PSSMscore += dict_AAPSSM[eachseq[i]][i]
            bindingscore = PSSMscore/int(length)
            IC50score = math.pow(50000, (0.8 - bindingscore)/1.6)
            dict_trainIC50.setdefault('pos', []).append((eachseq, IC50score))      
            
        dict_trainIC50['neg'] = OrderedDict()
        for eachset, negpeps in dict_setnegtrainpep.items():
            for eachseq in negpeps:            
                PSSMscore = 0
                for i in range(int(length)):
                    PSSMscore += dict_AAPSSM[eachseq[i]][i]
                bindingscore = PSSMscore/int(length)
                IC50score = math.pow(50000, (0.8 - bindingscore)/1.6)
                dict_trainIC50['neg'].setdefault(eachset, []).append((eachseq, IC50score))
            
        IC50matrix = {}
        IC50matrix['IC50'] = dict_AAPSSM
        
        if predictpeptide == None:
            dict_predictnameseqIC50 = None
        else:
            dict_predictnameseqIC50 = OrderedDict()
            for eachname, seqs in predictpeptide.items():
                dict_predictnameseqIC50[eachname] = OrderedDict()
                for eachseq in seqs:
                    PSSMscore = 0
                    for i in range(int(length)):
                        PSSMscore += dict_AAPSSM[eachseq[i]][i]
                    bindingscore = PSSMscore/int(length)
                    IC50score = math.pow(50000, (0.8 - bindingscore)/1.6)
                    dict_predictnameseqIC50[eachname][eachseq] = IC50score
                
        if testpeptide == None:
            dict_testIC50 = None
        else:
            dict_testIC50 = OrderedDict()
            postestpep, negtestpep = [], []
            for key, value in testpeptide.items():
                if value == '1':
                    postestpep.append(key)
                else:
                    negtestpep.append(key)            
            
            for eachseq in postestpep:    
                PSSMscore = 0
                for i in range(int(length)):
                    PSSMscore += dict_AAPSSM[eachseq[i]][i]
                bindingscore = PSSMscore/int(length)
                IC50score = math.pow(50000, (0.8 - bindingscore)/1.6)
                dict_testIC50[eachseq] = (IC50score, 1)               
            for eachseq in negtestpep:
                PSSMscore = 0
                for i in range(int(length)):
                    PSSMscore += dict_AAPSSM[eachseq[i]][i]
                bindingscore = PSSMscore/int(length)
                IC50score = math.pow(50000, (0.8 - bindingscore)/1.6)
                dict_testIC50[eachseq] = (IC50score, -1)
                            
        return(dict_trainIC50, dict_predictnameseqIC50, dict_testIC50, IC50matrix)

    else:            
        for HLA, matrixs in matrix.items():
            dict_AAPSSM = matrixs['IC50']
        
        lengthmatrix = len(dict_AAPSSM['A'])
        
        dict_predictnameseqIC50 = OrderedDict()
        for eachname, seqs in predictpeptide.items():
            dict_predictnameseqIC50[eachname] = OrderedDict()
            for eachseq in seqs:
                PSSMscore = 0
                for i in range(lengthmatrix):
                    PSSMscore += dict_AAPSSM[eachseq[i]][i]
                bindingscore = PSSMscore/lengthmatrix
                IC50score = math.pow(50000, (0.8 - bindingscore)/1.6)
                dict_predictnameseqIC50[eachname][eachseq] = IC50score
        
        return(dict_predictnameseqIC50)                           