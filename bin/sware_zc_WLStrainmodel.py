###generate WLS for train peptide
from collections import OrderedDict
import json, os, subprocess
from random import randint

def WLS(predictpeptide, length, trainpeptide, testpeptide, matrix):
    currentpath = os.getcwd()
    if matrix == None:
        while True:
            templatefolder = str(randint(0, 10000))
            if os.path.exists('%s/%s' % (currentpath, templatefolder)):
                continue
            else:
                break
        os.system('mkdir %s/%s' % (currentpath, templatefolder))

        postrainpep, dict_setnegtrainpep = [], OrderedDict()
        for eachpos in trainpeptide['pos']:
            postrainpep.append(eachpos)
        for eachset, negpeps in trainpeptide['neg'].items():
            for eachneg in negpeps:
                dict_setnegtrainpep.setdefault(eachset, []).append(eachneg)
        
        with open('%s/%s/trainpeptide.fasta' % (currentpath, templatefolder), 'w') as f:
            rownum = 1
            for eachtrainpeptide in postrainpep:
                f.write('>' + str(rownum) + '\n')
                f.write(eachtrainpeptide + '\n')
                rownum += 1
        
        subprocess.call(['%s/source/weblogo/seqlogo' % currentpath, ('-f %s/%s/trainpeptide.fasta' % (currentpath, templatefolder)), ('-o %s/%s/trainpeptide' % (currentpath, templatefolder))])        

        dict_WLS = {
            'A': [], 'R': [], 'N': [], 'D': [], 'C': [],
            'Q': [], 'E': [], 'G': [], 'H': [], 'I': [],
            'L': [], 'K': [], 'M': [], 'F': [], 'P': [],
            'S': [], 'T': [], 'W': [], 'Y': [], 'V': [],
            }
        with open('%s/%s/trainpeptide.eps' % (currentpath, templatefolder), 'r') as f:
            valuelist = []
            startpoint = 0
            addempty = 1
            for eachline in f:
                if ') startstack' in eachline:            
                    startpoint += 1
                elif startpoint == 1 and 'endstack' not in eachline:
                    valuelist.append(eachline.split(' ')[1:3])
                elif startpoint == 1 and 'endstack' in eachline:
                    for eachvalue in valuelist:
                        dict_WLS[eachvalue[1][1]].append(eachvalue[0])
                    for eachAAtemp, eachwistemp in dict_WLS.items():
                        if len(eachwistemp) != addempty:
                            eachwistemp.append('0')                        
                    startpoint -= 1
                    valuelist = []
                    addempty += 1
                else:
                    continue
                   
        dict_trainWLS = OrderedDict()
        for eachseq in postrainpep:
            scoreseq = 0
            for i in range(int(length)):
                if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                    scoreseq += 0
                else:
                    scoreseq += float(dict_WLS[eachseq[i]][i])
            dict_trainWLS.setdefault('pos', []).append((eachseq, scoreseq))
               
        dict_trainWLS['neg'] = OrderedDict()
        for eachset, negpeps in dict_setnegtrainpep.items():
            for eachseq in negpeps:        
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += float(dict_WLS[eachseq[i]][i])
                dict_trainWLS['neg'].setdefault(eachset, []).append((eachseq, scoreseq))
               
        WLSmatrix = {}
        WLSmatrix['WLS'] = dict_WLS
        os.system('rm -r %s/%s' % (currentpath, templatefolder))
           
        if predictpeptide == None:
            dict_predictnameseqWLS = None
        else:
            dict_predictnameseqWLS = OrderedDict()
            for eachname, seqs in predictpeptide.items():
                dict_predictnameseqWLS[eachname] = OrderedDict()
                for eachseq in seqs:
                    scoreseq = 0
                    for i in range(int(length)):
                        if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                            scoreseq += 0
                        else:
                            scoreseq += float(dict_WLS[eachseq[i]][i])
                    dict_predictnameseqWLS[eachname][eachseq] = scoreseq            
               
        if testpeptide == None:
            dict_testWLS = None
        else:
            dict_testWLS = OrderedDict()
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
                        scoreseq += float(dict_WLS[eachseq[i]][i])
                dict_testWLS[eachseq] = (scoreseq, 1)
            for eachseq in negtestpep:
                scoreseq = 0
                for i in range(int(length)):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += float(dict_WLS[eachseq[i]][i])
                dict_testWLS[eachseq] = (scoreseq, -1)               
                              
        return(dict_trainWLS, dict_predictnameseqWLS, dict_testWLS, WLSmatrix)
    
    else:        
        for HLA, matrixs in matrix.items():
            dict_WLS = matrixs['WLS']
            
        lengthmatrix = len(dict_WLS['A']) 
            
        dict_predictnameseqWLS = OrderedDict()
        for eachname, seqs in predictpeptide.items():
            dict_predictnameseqWLS[eachname] = OrderedDict()
            for eachseq in seqs:
                scoreseq = 0
                for i in range(lengthmatrix):
                    if eachseq[i] in ['B', 'J', 'O', 'U', 'X', 'Z']:
                        scoreseq += 0
                    else:
                        scoreseq += float(dict_WLS[eachseq[i]][i])
                dict_predictnameseqWLS[eachname][eachseq] = scoreseq                 
                    
        return(dict_predictnameseqWLS)