###calculate performance
from collections import OrderedDict
import numpy as np
import math, os
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, auc

def calper(trainperformancefile, templatefolder, resultfolder):
    currentpath = os.getcwd()
    if trainperformancefile != None:
        filename = 'model.txt' 
        piclabel = 'Training ROC curve'
        picname = 'Model_ROC'
    else:
        filename = 'testresult.txt'
        piclabel = 'Test ROC curve'
        picname = 'Test_ROC'
        
    dict_performanceCV = OrderedDict()
    labelscorelist = []
    filelist = os.listdir('%s/%s' % (currentpath, templatefolder))
    modelnum = 0
    for eachfile in filelist:
        if filename in eachfile:
            modelnum += 1
            with open('%s/%s/%s' % (currentpath, templatefolder, eachfile), 'r') as f:
                checkline = 0
                tp, tn, fn, fp = 0, 0, 0, 0
                for eachline in f:
                    if 'inst#' in eachline:
                        checkline += 1
                    elif checkline == 1 and len(eachline) > 1 and '+' not in eachline:                    
                        labelcontent = eachline.split(':')[1]
                        if '-' not in labelcontent:
                            tp += 1
                            ###如果是正例，预测为正例，直接用，总之是要得到预测为正例的p来画AUC
                            labelscorelist.append(('1', eachline.split(' ')[-2]))
                        else:
                            tn += 1
                            ###如果是负例，预测为负例，用1-p
                            scorethispep = 1 - float(eachline.split(' ')[-2])
                            labelscorelist.append(('-1', str(scorethispep)))
                    elif checkline == 1 and len(eachline) > 1 and '+' in eachline:
                        labelcontent = eachline.split(':')[1]
                        if '-' not in labelcontent:
                            fn += 1
                            ###如果是正例，预测为负列，用1-p
                            scorethispep = 1 - float(eachline.split(' ')[-2])
                            labelscorelist.append(('1', str(scorethispep)))
                        else:
                            fp += 1
                            ###如果是负例，预测为正例，直接用
                            labelscorelist.append(('-1', eachline.split(' ')[-2]))
                sensitivity = round(tp/(tp + fn), 3)
                specificity = round(tn/(tn + fp), 3)      
                accuracy = round((tp + tn)/(tp + fn + fp + tn), 3)                    
                mcc = round((tp * tn - fp * fn)/math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)), 3)
    ###只把正例当作pos计算tpr，fpr    
                dict_performanceCV.setdefault('Sensitivity', []).append(sensitivity)
                dict_performanceCV.setdefault('Specificity', []).append(specificity)
                dict_performanceCV.setdefault('Accuracy', []).append(accuracy)
                dict_performanceCV.setdefault('MCC', []).append(mcc)
    dict_performanceCVcal = OrderedDict()
    for eachmetric, valuelist in dict_performanceCV.items():
        dict_performanceCVcal[eachmetric] = round(sum(valuelist)/modelnum, 3)
    labels, probs = [], []
    for eachtuple in labelscorelist:
        labels.append(int(eachtuple[0]))
        probs.append(float(eachtuple[1]))
         
    y_labels = np.array(labels)
    y_scores = np.array(probs)
    dict_performanceCVcal['Area Under the Curve (AUC)'] = round(roc_auc_score(y_labels, y_scores), 3)
    if dict_performanceCVcal['Area Under the Curve (AUC)'] < 0.5:
        labels_reverse = []
        for eachlabel in labels:
            if eachlabel == 1:
                labels_reverse.append(-1)
            else:
                labels_reverse.append(1)
        labels = labels_reverse
        y_labels = np.array(labels)
        dict_performanceCVcal['Area Under the Curve (AUC)'] = round(1 - dict_performanceCV['Area Under the Curve (AUC)'], 3)

    ###draw ROC curve
    fpr, tpr, threshold = roc_curve(y_labels, y_scores)
    figure, ax = plt.subplots(figsize = (10, 10))
    plt.tick_params(which = 'major', width = 1, length = 10, pad = 10, labelsize = 16, zorder = 100)
    ax.set(xlim = (-0.05, 1.05), ylim = (-0.05, 1.05))
    x_major_locator = plt.MultipleLocator(0.1)
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(x_major_locator)    
    ax.plot(fpr, tpr, color = 'red', lw = 3, label = '%s (AUC) = %s' % (piclabel, str(dict_performanceCVcal['Area Under the Curve (AUC)'])))
    ax.plot(ax.get_xlim(), ax.get_ylim(), color = 'darkgray', lw = 1, linestyle = 'dotted')    
    plt.grid(True, linestyle = 'dotted', color = 'darkgray')
    plt.xlabel('False Positive Rate', family = 'Tahoma', size = 20, weight = 'light', labelpad = 16)
    plt.ylabel('True Positive Rate', family = 'Tahoma', size = 20, weight = 'light', labelpad = 16)
    plt.legend(loc = 'lower right', fontsize = 'xx-large', shadow = True, edgecolor = 'black', borderaxespad = 2)
    figure.savefig('%s/%s/%s.png' % (currentpath, resultfolder, picname))
    
    return(dict_performanceCVcal)