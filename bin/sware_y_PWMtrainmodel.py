###generate PWM for train peptide
from collections import OrderedDict
import json, math

def PWM(predictpeptide, length, trainpeptide, testpeptide, matrix):
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
                if dict_AAfre[eachAA] == 0:
                    eachvallist.append(1/len(postrainpep))#set the minimum value
                else:
                    eachvallist.append(dict_AAfre[eachAA]/len(postrainpep))       
    
        dict_AAPWM = {
                'A': [], 'R': [], 'N': [], 'D': [], 'C': [],
                'Q': [], 'E': [], 'G': [], 'H': [], 'I': [],
                'L': [], 'K': [], 'M': [], 'F': [], 'P': [],
                'S': [], 'T': [], 'W': [], 'Y': [], 'V': [],
                }    
        for eachAA, eachAF in dict_AAlenfre.items():
            for eachposAF in eachAF:
                powerval = (eachposAF)/0.05
                posPWM = math.log(powerval, 2)
                dict_AAPWM[eachAA].append(posPWM)
        
        dict_trainPWM = OrderedDict()
        for eachseq in postrainpep:
            scoreseq = 0
            for i in range(int(length)):
                if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                    scoreseq += 0
                else:
                    scoreseq += dict_AAPWM[eachseq[i]][i]
            PWMscore = scoreseq/int(length)
            dict_trainPWM.setdefault('pos', []).append((eachseq, PWMscore))
            
        dict_trainPWM['neg'] = OrderedDict()
        for eachset, negpeps in dict_setnegtrainpep.items():
            for eachseq in negpeps:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAPWM[eachseq[i]][i]
                PWMscore = scoreseq/int(length)
                dict_trainPWM['neg'].setdefault(eachset, []).append((eachseq, PWMscore))
            
        PWMmatrix = {}
        PWMmatrix['PWM'] = dict_AAPWM

        if predictpeptide == None:
            dict_predictnameseqPWM = None
        else:
            dict_predictnameseqPWM = OrderedDict()
            for eachname, seqs in predictpeptide.items():
                dict_predictnameseqPWM[eachname] = OrderedDict()
                for eachseq in seqs:
                    scoreseq = 0
                    for i in range(int(length)):
                        if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                            scoreseq += 0
                        else:
                            scoreseq += dict_AAPWM[eachseq[i]][i]
                    PWMscore = scoreseq/int(length) 
                    dict_predictnameseqPWM[eachname][eachseq] = PWMscore
        
        if testpeptide == None:
            dict_testPWM = None
        else:
            dict_testPWM = OrderedDict()
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
                        scoreseq += dict_AAPWM[eachseq[i]][i]
                PWMscore = scoreseq/int(length)
                dict_testPWM[eachseq] = (PWMscore, 1)
            for eachseq in negtestpep:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAPWM[eachseq[i]][i]
                PWMscore = scoreseq/int(length)
                dict_testPWM[eachseq] = (PWMscore, -1)            
                        
        return(dict_trainPWM, dict_predictnameseqPWM, dict_testPWM, PWMmatrix)
    
    else:        
        for HLA, matrixs in matrix.items():
            dict_AAPWM = matrixs['PWM']        
            
        lengthmatrix = len(dict_AAPWM['A']) 
            
        dict_predictnameseqPWM = OrderedDict()
        for eachname, seqs in predictpeptide.items():
            dict_predictnameseqPWM[eachname] = OrderedDict()
            for eachseq in seqs:
                scoreseq = 0
                for i in range(lengthmatrix):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += dict_AAPWM[eachseq[i]][i]
                PWMscore = scoreseq/lengthmatrix
                dict_predictnameseqPWM[eachname][eachseq] = PWMscore
                
        return(dict_predictnameseqPWM)    