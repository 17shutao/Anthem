###make prediction
import sys, os, re, math, time, json
from collections import OrderedDict
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_w_writeresult, sware_zn_seqlogo

def pred(features, length, fastafile = None, thresholdvalueadjust = None):
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

    while True:
        resultfolder = time.strftime("%Y%m%d%H%M%S", time.localtime(time.time()))
        if os.path.exists('%s/%s' % (currentpath, resultfolder)):
            continue
        else:
            break
    os.system('mkdir %s/%s' % (currentpath, resultfolder))
    
    ###get model's best features
    with open('%s/source/20200511_dict_lenHLAnumselect.json' % currentpath, 'r') as f:
        fstr = f.read()
        dict_lenHLAnumselect = json.loads(fstr, object_pairs_hook = OrderedDict)    

    ###get model's specificity
    with open('%s/source/trainA1DEperformance.json' % currentpath, 'r') as f:
        fstr = f.read()
        dict_lenHLAperformance = json.loads(fstr, object_pairs_hook = OrderedDict)

    dict_HLAspeci = OrderedDict()
    for eachHLA, names in features.items():
        dict_HLAspeci[eachHLA] = sum(dict_lenHLAperformance[str(length)][eachHLA]['Specificity'])/5    
    
    ###create arff
    dict_HLAnameseqscore = OrderedDict()
    for eachHLA, names in features.items():
        dict_HLAnameseqscore[eachHLA] = OrderedDict()
        for eachname, seqfea in names.items():
            dict_HLAnameseqscore[eachHLA][eachname] = OrderedDict()
            for eachmodelnum in range(5):
                featuretoselect = list(dict_lenHLAnumselect[length][eachHLA][str(eachmodelnum)])            
                with open('%s/%s/%s_%s_%s.arff' % (currentpath, templatefolder, eachHLA, eachname, str(eachmodelnum)), 'w') as f:
                    f.write('@RELATION peptide\n')
                    f.write('\n')
                    for eachfeaturetoselect in featuretoselect:   
                        f.write('@ATTRIBUTE' + ' ' + featurelist[int(eachfeaturetoselect) - 1] + ' ' + 'REAL\n')  
                    f.write('@ATTRIBUTE' + ' ' + 'class' + ' ' + '{1, -1}\n')
                    f.write('\n')
                    f.write('@DATA\n')
                    for eachseq, feas in seqfea.items():
                        for eachfeaturetoselect in featuretoselect:
                            f.write(str(feas[int(eachfeaturetoselect) - 1]) + ' ')
                        f.write('?')
                        f.write('\n') 

    ###use model
    filelist = os.listdir('%s/%s' % (currentpath, templatefolder))
    os.chdir('source/weka-3-9-3')
    dict_HLAfastanamescorelist = OrderedDict()
    for eachfile in filelist:
        if '.arff' in eachfile:
            HLAname = eachfile.split('_')[0]
            fastaname = eachfile.split('_')[1]
            if HLAname not in dict_HLAfastanamescorelist.keys():
                dict_HLAfastanamescorelist[HLAname] = OrderedDict()
            if fastaname not in dict_HLAfastanamescorelist[HLAname].keys():
                dict_HLAfastanamescorelist[HLAname][fastaname] = OrderedDict()            
            modelname = length + '_' + eachfile[4] + eachfile[6:8] + eachfile[9:11] + '_' + (eachfile.split('_')[-1]).split('.')[0]
            models = os.listdir('%s/models/%s' % (currentpath, length))
            if modelname in models:
                resultend = modelname.split('_')[1] + '_' + modelname.split('_')[-1]                            
                try:
                    attributenum = 0
                    with open('%s/%s/%s' % (currentpath, templatefolder, eachfile), 'r') as f:
                        rownum = 0
                        for eachline in f:
                            if rownum < 10:
                                if '@ATTRIBUTE' in eachfile:
                                    attributenum += 1
                                rownum += 1
                    os.system('java weka.classifiers.meta.FilteredClassifier -p %s -l %s/models/%s/%s -T %s/%s/%s > %s/%s/%s_result_%s.txt' % (str(attributenum), currentpath, length, modelname, currentpath, templatefolder, eachfile, currentpath, templatefolder, eachfile.split('.')[0], resultend))                    
                    with open('%s/%s/%s_result_%s.txt' % (currentpath, templatefolder, eachfile.split('.')[0], resultend), 'r') as f:                
                        checkline = 0
                        for eachline in f:
                            eachlenwithoutblank = eachline.split(' ')
                            while '' in eachlenwithoutblank:
                                eachlenwithoutblank.remove('')                            
                            if 'inst#' in eachline:
                                checkline += 1
                            elif checkline == 1 and len(eachline) > 1 and '-' not in eachline:
                                dict_HLAfastanamescorelist[HLAname][fastaname].setdefault(eachlenwithoutblank[0], []).append(float(eachline.split(' ')[-2]))
                            elif checkline == 1 and len(eachline) > 1 and '-' in eachline:
                                dict_HLAfastanamescorelist[HLAname][fastaname].setdefault(eachlenwithoutblank[0], []).append((1 - float(eachline.split(' ')[-2])))
                except:
                    os.system('rm -r %s/%s' % (currentpath, templatefolder))
                    os.system('rm %s/*_result.txt' % currentpath)
                    sys.exit(1)
                    print('Prediction failed, please check the input file')
    for eachHLA, names in dict_HLAfastanamescorelist.items():
        if thresholdvalueadjust == None:
            HLAspecificity = dict_HLAspeci[eachHLA]
        else:
            HLAspecificity = dict_HLAspeci[eachHLA]*thresholdvalueadjust
        for eachname, peptidescorelist in names.items():
            pepscorecal, binderlist = [], []
            for eachpep, scores in peptidescorelist.items():
                scoresum = round(sum(scores)/5, 3)
                pepscorecal.append(str(scoresum))
                if scoresum >= HLAspecificity:
                    binderlist.append('yes')
                else:
                    binderlist.append(' ')    
            peptidelist = list(features[eachHLA][eachname].keys())
            dict_HLAnameseqscore[eachHLA][eachname] = list(zip(peptidelist, pepscorecal, binderlist))            
    os.chdir(os.path.pardir)
    os.chdir(os.path.pardir)

    ###for fasta format
    binderforlogo = OrderedDict()
    for eachHLA, names in dict_HLAnameseqscore.items():
        binderforlogo[eachHLA] = OrderedDict()
        for eachname, pepscorelevel in names.items():
            if eachname == 'input-sequences':
                eachnameforseqlogo = 'motif-from-predicted-binders'
            else:
                eachnameforseqlogo = 'motif-from-predicted-binders-of-%s' % eachname            
            binderforlogo[eachHLA][eachnameforseqlogo] = []
            for eachtuple in pepscorelevel:
                if eachtuple[2] != ' ':
                    binderforlogo[eachHLA][eachnameforseqlogo].append(eachtuple[0])
    
    # sware_zn_seqlogo.seqlogo(templatefolder, resultfolder, binderforlogo, None)
    # os.chdir(os.path.pardir)
    
    ###create binder visualizatino for fasta
    if fastafile != None:
        dict_HLAnameposlevel = OrderedDict()
        for eachHLA, names in dict_HLAnameseqscore.items():
            dict_HLAnameposlevel[eachHLA] = OrderedDict()
            for eachname, pepscorelevel in names.items():
                dict_HLAnameposlevel[eachHLA][eachname] = []
                thisfastaseq = fastafile[eachHLA][eachname]
                samebinderinseq = 0
                for eachtuple in pepscorelevel:
                    if eachtuple[2] != ' ':
                        binderseq = eachtuple[0]
                        matchedseq = re.findall(binderseq, thisfastaseq)
                        if len(matchedseq) == 1:
                            pos = thisfastaseq.find(binderseq) + math.ceil(int(length)/2)
                            dict_HLAnameposlevel[eachHLA][eachname].append((pos, eachtuple[2], eachtuple[1]))
                        else:
                            samebinderinseq += 1
                            poslist = list(re.finditer(binderseq, thisfastaseq))
                            for eachposlist in poslist:
                                postuple = eachposlist.span()
                                startpos = postuple[0]
                                pos = startpos + math.ceil(int(length)/2)
                                dict_HLAnameposlevel[eachHLA][eachname].append((pos, eachtuple[2], eachtuple[1]))
                if samebinderinseq != 0:
                    poslevelsort = sorted(dict_HLAnameposlevel[eachHLA][eachname], key = lambda t : t[0])
                    dict_HLAnameposlevel[eachHLA][eachname] = []
                    for eachposlevel in poslevelsort:
                        dict_HLAnameposlevel[eachHLA][eachname].append(eachposlevel)

    os.system('rm -r %s/%s' % (currentpath, templatefolder))

    if fastafile != None:
        sware_w_writeresult.writeresult(resultfolder, length, None, None, dict_HLAnameseqscore, None, None, 'prediction', None, dict_HLAnameposlevel, thresholdvalueadjust)
    else:
        sware_w_writeresult.writeresult(resultfolder, length, None, None, dict_HLAnameseqscore, None, None, 'prediction', None, None, thresholdvalueadjust)
 