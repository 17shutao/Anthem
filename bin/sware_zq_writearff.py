###write the arff file
import os

def writearff(predic_feas, templatefolder, length, newHLA, featurelist, dict_setfeaselected, negprovided):
    if negprovided == None:
        modelnum = 5
    else:
        modelnum = 1    
    currentpath = os.getcwd()
    
    for eachname, pepfea in predic_feas.items():
        for eachmodelnum in range(modelnum):
            featuretoselect = list(dict_setfeaselected[str(eachmodelnum)])        
            with open('%s/%s/predict_%s_%s_%s_%s.arff' % (currentpath, templatefolder, length, newHLA, eachname, str(eachmodelnum)), 'w') as f:
                f.write('@RELATION peptide\n')
                f.write('\n')
                for eachfeaturetoselect in featuretoselect:  
                    f.write('@ATTRIBUTE' + ' ' + featurelist[int(eachfeaturetoselect) - 1] + ' ' + 'REAL\n')  
                f.write('@ATTRIBUTE' + ' ' + 'class' + ' ' + '{1, -1}\n')
                f.write('\n')
                f.write('@DATA\n')
                for eachseq, feas in pepfea.items():
                    for eachfeaturetoselect in featuretoselect:
                        f.write(str(feas[int(eachfeaturetoselect) - 1]) + ' ')
                    f.write('?')
                    f.write('\n')        
