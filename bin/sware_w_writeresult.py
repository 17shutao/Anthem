###write the prediction result
import os, re, sys
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_zo_writebinderlocation

def writeresult(resultfolder, length, newHLA, dict_nameseqscore, dict_HLAnameseqscore, dict_performanceCV, dict_testperformanceCV, mode, CVtime, dict_HLAnameposlevel, thresholdvalueadjust):
    currentpath = os.getcwd()
    if mode == 'prediction':
        with open('%s/%s/length_%s_prediction_result.txt' % (currentpath, resultfolder, length), 'w') as f:
            f.write('HLA allele' + '\t' + 'Set' + '\t' + '\t' + 'Peptide' + '\t' + '\t' + 'Binder?' + '\t' + 'Prediction_score' + '\n')
            f.write('---------------------------------------------------------------------' + '\n')
            for eachHLA, names in dict_HLAnameseqscore.items():
                f.write(eachHLA + '\n')
                f.write('          ' + '\t' + '------------------------------------------------' + '\n')
                for eachname, pepscore in names.items():
                    f.write('          ' + '\t' + eachname + '\n')
                    numbinder = 0
                    for eachtuple in pepscore:
                        if eachtuple[2] == 'yes':
                            numbinder += 1
                        f.write('          ' + '\t' + '   ' + '\t' + '\t' + eachtuple[0] + '\t' + '\t' + eachtuple[2] + '\t' + '\t' + eachtuple[1] + '\n')
                    f.write('          ' + '\t' + '------------------------------------------------' + '\n')        
                    if eachname == 'input-sequences':
                        if thresholdvalueadjust != None:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total input sequences from allele %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachHLA, str(int(thresholdvalueadjust*100))) + '\n')
                        else:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total input sequences from allele %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachHLA, '100') + '\n')
                    else:
                        if thresholdvalueadjust != None:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total segments of protein %s from allele %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachname, eachHLA, str(int(thresholdvalueadjust*100))) + '\n')
                        else:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total segments of protein %s from allele %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachname, eachHLA, '100') + '\n')        
                    f.write('\n')
        if dict_HLAnameposlevel != None:
            sware_zo_writebinderlocation.writebinder(resultfolder, length, None, dict_HLAnameposlevel, None)

    elif mode == 'TrainYourModel':
        with open('%s/%s/trainmodel_result.txt' % (currentpath, resultfolder), 'w') as f:
            f.write('Model Performance' + '\t' + '\t' + 'Cross validation: %s' % CVtime + '\n')
            for key, value in dict_performanceCV.items():
                f.write(key + '\t' + '\t' + str(value) + '\n') 
            if dict_testperformanceCV != None:
                f.write('\n')
                f.write('Test performance' + '\t' + '\t' + 'Test Data' + '\n')
                for key, value in dict_testperformanceCV.items():
                    f.write(key + '\t' + '\t' + str(value) + '\n')
            if dict_nameseqscore != None:                
                f.write('\n')
                f.write('Prediction result of predict file')
                f.write('\n')                
                f.write('Model name' + '\t' + 'Set' + '\t' + '\t' + 'Peptide' + '\t' + '\t' + 'Binder?' + '\t' + 'Prediction_score' + '\n')
                f.write('---------------------------------------------------------------------' + '\n')
                f.write(newHLA + '\n')
                f.write('          ' + '\t' + '------------------------------------------------' + '\n')        
                for eachname, pepscore in dict_nameseqscore.items():
                    f.write('          ' + '\t' + eachname + '\n')
                    numbinder = 0
                    for eachtuple in pepscore:
                        if eachtuple[2] == 'yes':
                            numbinder += 1
                        f.write('          ' + '\t' + '   ' + '\t' + '\t' + eachtuple[0] + '\t' + '\t' + eachtuple[2] + '\t' + '\t' + eachtuple[1] + '\n')
                    f.write('          ' + '\t' + '------------------------------------------------' + '\n')
                    if eachname == 'input-sequences':
                        if thresholdvalueadjust != None:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total input sequences from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), newHLA, str(int(thresholdvalueadjust*100))) + '\n')
                        else:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total input sequences from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), newHLA, '100') + '\n')
                    else:
                        if thresholdvalueadjust != None:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total segments of protein %s from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachname, newHLA, str(int(thresholdvalueadjust*100))) + '\n')
                        else:
                            f.write('          ' + '\t' + '%s predicted binders out of %s total segments of protein %s from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachname, newHLA, '100') + '\n')                    
                    f.write('\n')                
                if dict_HLAnameposlevel != None:
                    sware_zo_writebinderlocation.writebinder(resultfolder, length, newHLA, None, dict_HLAnameposlevel)

    elif mode == 'useYourOwnModel':
        with open('%s/%s/length_%s_%s_prediction_result.txt' % (currentpath, resultfolder, length, newHLA), 'w') as f:
            f.write('Model name' + '\t' + 'Set' + '\t' + '\t' + 'Peptide' + '\t' + '\t' + 'Binder?' + '\t' + 'Prediction_score' + '\n')
            f.write('---------------------------------------------------------------------' + '\n')
            f.write(newHLA + '\n')
            f.write('          ' + '\t' + '------------------------------------------------' + '\n')        
            for eachname, pepscore in dict_nameseqscore.items():
                f.write('          ' + '\t' + eachname + '\n')
                numbinder = 0
                for eachtuple in pepscore:
                    if eachtuple[2] == 'yes':
                        numbinder += 1
                    f.write('          ' + '\t' + '   ' + '\t' + '\t' + eachtuple[0] + '\t' + '\t' + eachtuple[2] + '\t' + '\t' + eachtuple[1] + '\n')
                f.write('          ' + '\t' + '------------------------------------------------' + '\n')        
                if eachname == 'input-sequences':
                    if thresholdvalueadjust != None:
                        f.write('          ' + '\t' + '%s predicted binders out of %s total input sequences from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), newHLA, str(int(thresholdvalueadjust*100))) + '\n')
                    else:
                        f.write('          ' + '\t' + '%s predicted binders out of %s total input sequences from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), newHLA, '100') + '\n')
                else:
                    if thresholdvalueadjust != None:
                        f.write('          ' + '\t' + '%s predicted binders out of %s total segments of protein %s from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachname, newHLA, str(int(thresholdvalueadjust*100))) + '\n')
                    else:
                        f.write('          ' + '\t' + '%s predicted binders out of %s total segments of protein %s from model %s (at threshold of %s%% specificity)' % (str(numbinder), str(len(pepscore)), eachname, newHLA, '100') + '\n')                
                f.write('\n')        
        if dict_HLAnameposlevel != None:
            sware_zo_writebinderlocation.writebinder(resultfolder, length, newHLA, None, dict_HLAnameposlevel)                