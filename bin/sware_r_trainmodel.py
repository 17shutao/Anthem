###use train peptide to train model
import sys, os, re, math, time, json, tarfile
from collections import OrderedDict
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_w_writeresult, sware_zn_seqlogo, sware_zq_writearff, sware_zr_calperformance

def trainmodel(predic_feas, length, resultfolder, CVtime, train_feaslabel, newHLA, modelpath, fastafile, test_feaslabel, thresholdvalueadjust, allmatrix, negprovided):
    featurelist = ['AAF', 'WLS', 'PSSM', 'PWM', 'BSI62']
    currentpath = os.getcwd()
    ###create folder
    while True:
        templatefolder = time.strftime("%Y%m%d%H%M%S", time.localtime(time.time()))
        if os.path.exists('%s/%s' % (currentpath, templatefolder)):
            continue
        else:
            break
    os.system('mkdir %s/%s' % (currentpath, templatefolder))        
   
   
    if train_feaslabel != None:
        os.system('mkdir /%s/%s/%s' % (currentpath, templatefolder, newHLA))
        os.chdir('source/weka-3-9-3')  
        ###train arff
        for eachset, seqfeaslist in train_feaslabel['neg'].items():
            with open('%s/%s/train_%s_beforeselect.arff' % (currentpath, templatefolder, str(eachset)), 'w') as f:
                f.write('@RELATION peptide\n')
                f.write('\n')
                for featurename in featurelist:   
                    f.write('@ATTRIBUTE' + ' ' + featurename + ' ' + 'REAL\n')  
                f.write('@ATTRIBUTE' + ' ' + 'class' + ' ' + '{1, -1}\n')
                f.write('\n')
                f.write('@DATA\n')           
                for eachtupleseqfeas in train_feaslabel['pos']:
                    for eachfea in eachtupleseqfeas[1]:
                        f.write(str(eachfea) + ' ')
                    f.write('1')
                    f.write('\n')
                for eachtupleseqfeas in seqfeaslist:
                    for eachfea in eachtupleseqfeas[1]:
                        f.write(str(eachfea) + ' ')
                    f.write('-1')
                    f.write('\n')         
        filelist = os.listdir('%s/%s' % (currentpath, templatefolder))
          
        ###feature selection
        for eachfile in filelist:
            if 'beforeselect' in eachfile:
                filename = 'train' + '_' + eachfile.split('_')[1] + '_afterselect.txt'
                os.system('java weka.attributeSelection.WrapperSubsetEval -i %s/%s/%s -s "weka.attributeSelection.BestFirst -D 2" -E AUC -B weka.classifiers.meta.FilteredClassifier -- -F weka.filters.supervised.attribute.Discretize -W weka.classifiers.bayes.AveragedNDependenceEstimators.A1DE > %s/%s/%s' % (currentpath, templatefolder, eachfile, currentpath, templatefolder, filename))
                
        filelist = os.listdir('%s/%s' % (currentpath, templatefolder))        
          
        ###get feature selection information
        dict_setfeaselected = OrderedDict()
        for eachfile in filelist:
            if 'afterselect' in eachfile:
                modelnum = eachfile.split('_')[1]
                with open('%s/%s/%s' % (currentpath, templatefolder, eachfile), 'r') as f:
                    for eachline in f:
                        if 'Selected' in eachline:
                            templist = eachline.split(' ')
                            while '' in templist:
                                templist.remove('')
                            featureselected = ''.join(templist[2].split(','))
                            dict_setfeaselected[modelnum] = featureselected
        allmatrix[newHLA]['modelfeature'] = dict_setfeaselected          
          
        ###write arff after feature selection
        for eachset, seqfeaslist in train_feaslabel['neg'].items():
            featuretoselect = list(dict_setfeaselected[eachset])
            with open('%s/%s/train_%s_selected.arff' % (currentpath, templatefolder, str(eachset)), 'w') as f:
                f.write('@RELATION peptide\n')
                f.write('\n')
                for eachfeaturetoselect in featuretoselect:
                    f.write('@ATTRIBUTE' + ' ' + featurelist[int(eachfeaturetoselect) - 1] + ' ' + 'REAL\n')  
                f.write('@ATTRIBUTE' + ' ' + 'class' + ' ' + '{1, -1}\n')
                f.write('\n')
                f.write('@DATA\n')           
                for eachtupleseqfeas in train_feaslabel['pos']:
                    for eachfeaturetoselect in featuretoselect:
                        f.write(str(eachtupleseqfeas[1][int(eachfeaturetoselect) - 1]) + ' ')
                    f.write('1')
                    f.write('\n')
                for eachtupleseqfeas in train_feaslabel['neg'][eachset]:
                    for eachfeaturetoselect in featuretoselect:
                        f.write(str(eachtupleseqfeas[1][int(eachfeaturetoselect) - 1]) + ' ')
                    f.write('-1')
                    f.write('\n')             
        filelist = os.listdir('%s/%s' % (currentpath, templatefolder))
        
        ###train model
        for eachfile in filelist:
            if 'selected' in eachfile:
                attributenum = 0
                with open('%s/%s/%s' % (currentpath, templatefolder, eachfile), 'r') as f:
                    rownum = 0
                    for eachline in f:
                        if rownum < 10:
                            if '@ATTRIBUTE' in eachline:
                                attributenum += 1
                            rownum += 1 
                modelname = '%s_%s_%s_model' % (length, newHLA, eachfile.split('_')[1])
                resultname = modelname + '.txt'
                if CVtime != 'LOO':
                    os.system('java weka.classifiers.meta.FilteredClassifier -x %s -F weka.filters.supervised.attribute.Discretize -W weka.classifiers.bayes.AveragedNDependenceEstimators.A1DE -d %s/%s/%s/%s -p %s -t %s/%s/%s > %s/%s/%s' % (CVtime, currentpath, templatefolder, newHLA, modelname, str(attributenum), currentpath, templatefolder, eachfile, currentpath, templatefolder, resultname))
                else:
                    os.system('java weka.classifiers.meta.FilteredClassifier -x %s -F weka.filters.supervised.attribute.Discretize -W weka.classifiers.bayes.AveragedNDependenceEstimators.A1DE -d %s/%s/%s/%s -p %s -t %s/%s/%s > %s/%s/%s' % (str(len(train_feaslabel)), currentpath, templatefolder, newHLA, modelname, str(attributenum), currentpath, templatefolder, eachfile, currentpath, templatefolder, resultname))
        os.chdir(os.path.pardir)
        os.chdir(os.path.pardir)
        dict_performanceCV = sware_zr_calperformance.calper('Yes', templatefolder, resultfolder)
        modelspecificity = dict_performanceCV['Specificity']
        allmatrix[newHLA]['Specificity'] = modelspecificity
        with open('%s/%s/%s/allmatrix.json' % (currentpath, templatefolder, newHLA), 'w') as f:
            json.dump(allmatrix, f)  
            
        os.chdir('/%s/%s/%s' % (currentpath, templatefolder, newHLA))
        with tarfile.open('/%s/%s/%s' % (currentpath, resultfolder, newHLA), 'w') as tar:
            listdir = os.listdir()
            for eachfile in listdir:
                if 'allmatrix.json' in eachfile or '_model' in eachfile:
                    tar.add(eachfile)
        os.chdir(os.path.pardir)
        os.chdir(os.path.pardir)
        
        if predic_feas != None:
            ###write predict arff
            sware_zq_writearff.writearff(predic_feas, templatefolder, length, newHLA, featurelist, dict_setfeaselected, negprovided)
            ###use model
            os.chdir('source/weka-3-9-3')
            filelist = os.listdir('%s/%s' % (currentpath, templatefolder))
            resultfile = []
            for eachfile in filelist:
                if 'predict' in eachfile:
                    modelindex = (eachfile.split('_')[-1]).split('.')[0]
                    eachname = eachfile.split('_')[3]
                    resultname = '%s_%s_%s_%s_result.txt' % (length, newHLA, eachname, modelindex)                    
                    modelname = '%s_%s_%s_model' % (length, newHLA, modelindex)
                    attributenum = len(dict_setfeaselected[modelindex]) + 1
                    resultfile.append(resultname)
                    os.system('java weka.classifiers.meta.FilteredClassifier -p %s -l %s/%s/%s/%s -T %s/%s/%s > %s/%s/%s' % (str(attributenum), currentpath, templatefolder, newHLA, modelname, currentpath, templatefolder, eachfile, currentpath, templatefolder, resultname))
            os.chdir(os.path.pardir)
            os.chdir(os.path.pardir)

        if test_feaslabel != None:
            os.chdir('source/weka-3-9-3')
            if negprovided == None:
                modelnum = 5
            else:
                modelnum = 1   
            for eachmodelnum in range(modelnum):
                featuretoselect = list(dict_setfeaselected[str(eachmodelnum)])                         
                with open('%s/%s/test_%s.arff' % (currentpath, templatefolder, str(eachmodelnum)), 'w') as f:
                    f.write('@RELATION peptide\n')
                    f.write('\n')
                    for eachfeaturetoselect in featuretoselect:   
                        f.write('@ATTRIBUTE' + ' ' + featurelist[int(eachfeaturetoselect) -1] + ' ' + 'REAL\n')  
                    f.write('@ATTRIBUTE' + ' ' + 'class' + ' ' + '{1, -1}\n')
                    f.write('\n')
                    f.write('@DATA\n')           
                    for eachseq, feaslabel in test_feaslabel.items():
                        for eachfeaturetoselect in featuretoselect:
                            f.write(str(feaslabel[0][int(eachfeaturetoselect) - 1]) + ' ')
                        f.write(str(feaslabel[1]))
                        f.write('\n')               
            ###use model
            filelist = os.listdir('%s/%s' % (currentpath, templatefolder))
            for eachfile in filelist:
                if 'test' in eachfile:            
                    modelindex = (eachfile.split('_')[-1]).split('.')[0]
                    resultname = '%s_testresult.txt' % modelindex
                    modelname = '%s_%s_%s_model' % (length, newHLA, modelindex)
                    attributenum = len(dict_setfeaselected[modelindex]) + 1
                    os.system('java weka.classifiers.meta.FilteredClassifier -p %s -l %s/%s/%s/%s -T %s/%s/%s > %s/%s/%s' % (str(attributenum), currentpath, templatefolder, newHLA, modelname, currentpath, templatefolder, eachfile, currentpath, templatefolder, resultname))
            os.chdir(os.path.pardir)
            os.chdir(os.path.pardir)
            dict_testperformanceCV = sware_zr_calperformance.calper(None, templatefolder, resultfolder)
        else:
            dict_testperformanceCV = None
    else:
        dict_performanceCV = None
        dict_testperformanceCV = None

    if modelpath != None:
        os.system('mv %s/model %s/%s' % (currentpath, currentpath, templatefolder))
        dict_setfeaselected = allmatrix[newHLA]['modelfeature']
        modelnum = len(dict_setfeaselected.keys())
        modelspecificity = allmatrix[newHLA]['Specificity']
        if modelnum == 1:
            sware_zq_writearff.writearff(predic_feas, templatefolder, length, newHLA, featurelist, dict_setfeaselected, 'negprovided')
        else:
            sware_zq_writearff.writearff(predic_feas, templatefolder, length, newHLA, featurelist, dict_setfeaselected, None)
        ###use model
        os.chdir('source/weka-3-9-3')
        filelist = os.listdir('%s/%s' % (currentpath, templatefolder))                
        resultfile = []
        for eachfile in filelist:
            if 'predict' in eachfile:
                modelindex = (eachfile.split('_')[-1]).split('.')[0]
                eachname = eachfile.split('_')[3]
                resultname = '%s_%s_%s_%s_result.txt' % (length, newHLA, eachname, modelindex)
                modelname = '%s_%s_%s_model' % (length, newHLA, modelindex)
                attributenum = len(dict_setfeaselected[modelindex]) + 1
                resultfile.append(resultname)
                os.system('java weka.classifiers.meta.FilteredClassifier -p %s -l %s/%s/model/%s -T %s/%s/%s > %s/%s/%s' % (str(attributenum), currentpath, templatefolder, modelname, currentpath, templatefolder, eachfile, currentpath, templatefolder, resultname))        
        os.chdir(os.path.pardir)
        os.chdir(os.path.pardir)

    if predic_feas != None:
        dict_nameseqscore = OrderedDict()
        dict_HLAfastanamescorelist = OrderedDict()
        for eachname, seqfea in predic_feas.items():
            dict_nameseqscore[eachname] = OrderedDict()
            dict_HLAfastanamescorelist[eachname] = OrderedDict()
            for eachresultfile in resultfile:
                if eachname == eachresultfile.split('_')[2]:
                    with open('%s/%s/%s' % (currentpath, templatefolder, eachresultfile), 'r') as f:
                        checkline = 0
                        binderlist = []
                        for eachline in f:
                            eachlenwithoutblank = eachline.split(' ')
                            while '' in eachlenwithoutblank:
                                eachlenwithoutblank.remove('')                    
                            if 'inst#' in eachline:
                                checkline += 1
                            elif checkline == 1 and len(eachline) > 1 and '-' not in eachline:
                                dict_HLAfastanamescorelist[eachname].setdefault(eachlenwithoutblank[0], []).append(float(eachline.split(' ')[-2]))
                            elif checkline == 1 and len(eachline) > 1 and '-' in eachline:
                                dict_HLAfastanamescorelist[eachname].setdefault(eachlenwithoutblank[0], []).append((1 - float(eachline.split(' ')[-2])))
        if thresholdvalueadjust == None:
            HLAspecificity = modelspecificity
        else:
            HLAspecificity = modelspecificity*thresholdvalueadjust        
        for eachname, peptidescorelist in dict_HLAfastanamescorelist.items():
            pepscorecal, binderlist = [], []
            for eachpep, scores in peptidescorelist.items():
                scoresum = round(sum(scores)/5, 3)
                pepscorecal.append(str(scoresum))
                if scoresum >= HLAspecificity:
                    binderlist.append('yes')
                else:
                    binderlist.append(' ')    
            peptidelist = list(predic_feas[eachname].keys())   
            dict_nameseqscore[eachname] = list(zip(peptidelist, pepscorecal, binderlist))  
    
        binderforlogo = OrderedDict()
        for eachname, pepscoreklevel in dict_nameseqscore.items():
            if eachname == 'input-sequences':
                eachnameforseqlogo = 'motif-from-predicted-binders'
            else:
                eachnameforseqlogo = 'motif-from-predicted-binders-of-%s' % eachname            
            binderforlogo[eachnameforseqlogo] = []
            for eachtuple in pepscoreklevel:
                if eachtuple[2] != ' ':
                    binderforlogo[eachnameforseqlogo].append(eachtuple[0])
        
        # sware_zn_seqlogo.seqlogo(templatefolder, resultfolder, None, binderforlogo)
        # os.chdir(os.path.pardir)
        
        ###create binder visualization for fasta
        if fastafile != None:
            dict_nameposlevel = OrderedDict()
            for eachname, pepscorelevel in dict_nameseqscore.items():
                dict_nameposlevel[eachname] = []
                thisfastaseq = fastafile[eachname]
                samebinderinseq = 0
                for eachtuple in pepscorelevel:
                    if eachtuple[2] != ' ':
                        binderseq = eachtuple[0]
                        matchedseq = re.findall(binderseq, thisfastaseq)
                        if len(matchedseq) == 1:
                            pos = thisfastaseq.find(binderseq) + math.ceil(int(length)/2)
                            dict_nameposlevel[eachname].append((pos, eachtuple[2], eachtuple[1]))
                        else:
                            samebinderinseq += 1
                            poslist = list(re.finditer(binderseq, thisfastaseq))
                            for eachposlist in poslist:
                                postuple = eachposlist.span()
                                startpos = postuple[0]
                                pos = startpos + math.ceil(int(length)/2)
                                dict_nameposlevel[eachname].append((pos, eachtuple[2], eachtuple[1]))
                if samebinderinseq != 0:
                    poslevelsort = sorted(dict_nameposlevel[eachname], key = lambda t : t[0])
                    dict_nameposlevel[eachname] = []
                    for eachposlevel in poslevelsort:
                        dict_nameposlevel[eachname].append(eachposlevel)
        else:
            fastafile, dict_nameposlevel = None, None
    else:
        dict_nameseqscore, dict_nameposlevel = None, None

    os.system('rm -r %s/%s' % (currentpath, templatefolder))

    if train_feaslabel != None:       
        sware_w_writeresult.writeresult(resultfolder, length, newHLA, dict_nameseqscore, None, dict_performanceCV, dict_testperformanceCV, 'TrainYourModel', CVtime, dict_nameposlevel, thresholdvalueadjust)
    else:        
        sware_w_writeresult.writeresult(resultfolder, length, newHLA, dict_nameseqscore, None, dict_performanceCV, dict_testperformanceCV, 'useYourOwnModel', CVtime, dict_nameposlevel, thresholdvalueadjust) 