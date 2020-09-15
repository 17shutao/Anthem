###combine feature for training new model
from collections import OrderedDict
import sys, os, re, time, json
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_q_AAFtrainmodel, sware_v_IC50trainmodel, sware_y_PWMtrainmodel, sware_zc_WLStrainmodel, sware_zk_BSI62trainmodel

def combfeanew(length,  predictpeptide, trainpeptide, newHLA, matrix, testpeptide):
    currentpath = os.getcwd()
    ###create folder to save model,matrix and prediction result
    while True:
        resultfolder = time.strftime("%Y%m%d%H%M%S", time.localtime(time.time()))
        if os.path.exists('%s/%s' % (currentpath, resultfolder)):
            continue
        else:
            break
    os.system('mkdir %s/%s' % (currentpath, resultfolder))
    
    if matrix == None:
        allmatrix = OrderedDict()
        allmatrix[newHLA] = OrderedDict()
    
        [train_AAF, dict_predictnameseqAAF, test_AAF, AAFmatrix] = sware_q_AAFtrainmodel.AAF(predictpeptide, length, trainpeptide, testpeptide, None)
        [train_IC50, dict_predictnameseqIC50, test_IC50, IC50matrix] = sware_v_IC50trainmodel.IC50(predictpeptide, length, trainpeptide, testpeptide, None)
        [train_PWM, dict_predictnameseqPWM, test_PWM, PWMmatrix] = sware_y_PWMtrainmodel.PWM(predictpeptide, length, trainpeptide, testpeptide, None)        
        [train_WLS, dict_predictnameseqWLS, test_WLS, WLSmatrix] = sware_zc_WLStrainmodel.WLS(predictpeptide, length, trainpeptide, testpeptide, None)
        [train_BSI62, dict_predictnameseqBSI62, test_BSI62, BSI62matrix] = sware_zk_BSI62trainmodel.BSI62(predictpeptide, length, trainpeptide, testpeptide, None)

        allmatrix[newHLA]['AAF'] = AAFmatrix['AAF']
        allmatrix[newHLA]['IC50'] = IC50matrix['IC50']
        allmatrix[newHLA]['PWM'] = PWMmatrix['PWM']
        allmatrix[newHLA]['WLS'] = WLSmatrix['WLS']
        allmatrix[newHLA]['BSI62'] = BSI62matrix['BSI62']
    
        dict_trainfeaslabel = OrderedDict()
        for eachtrainfeaitems in enumerate(train_AAF['pos']):
            allfeatures = [eachtrainfeaitems[1][1], train_WLS['pos'][eachtrainfeaitems[0]][1], train_IC50['pos'][eachtrainfeaitems[0]][1], train_PWM['pos'][eachtrainfeaitems[0]][1], train_BSI62['pos'][eachtrainfeaitems[0]][1]]
            dict_trainfeaslabel.setdefault('pos', []).append((eachtrainfeaitems[1][0], allfeatures))

        dict_trainfeaslabel['neg'] = OrderedDict()
        for eachset, seqfealist in train_AAF['neg'].items():
            for eachtrainfeaitems in enumerate(seqfealist):
                allfeatures = [eachtrainfeaitems[1][1], train_WLS['neg'][eachset][eachtrainfeaitems[0]][1], train_IC50['neg'][eachset][eachtrainfeaitems[0]][1], train_PWM['neg'][eachset][eachtrainfeaitems[0]][1], train_BSI62['neg'][eachset][eachtrainfeaitems[0]][1]]
                dict_trainfeaslabel['neg'].setdefault(eachset, []).append((eachtrainfeaitems[1][0], allfeatures))

        dict_predictfeas, dict_testfeaslabel = OrderedDict(), OrderedDict()
        if predictpeptide == None:
            dict_predictfeas = None
        else:
            for eachname, prefea in dict_predictnameseqAAF.items():
                dict_predictfeas[eachname] = OrderedDict()  
                for eachpredictseq, fea in prefea.items():
                    allfeatures = [fea, dict_predictnameseqWLS[eachname][eachpredictseq], dict_predictnameseqIC50[eachname][eachpredictseq],
                                   dict_predictnameseqPWM[eachname][eachpredictseq], dict_predictnameseqBSI62[eachname][eachpredictseq]]
                    dict_predictfeas[eachname][eachpredictseq] = allfeatures
                
        if testpeptide == None:
            dict_testfeaslabel = None
        else:
            for eachtest, fealabel in test_AAF.items():
                allfeatures = [fealabel[0], test_WLS[eachtest][0], test_IC50[eachtest][0], test_PWM[eachtest][0], test_BSI62[eachtest][0]]
                dict_testfeaslabel[eachtest] = (allfeatures, fealabel[1])

        return(dict_trainfeaslabel, dict_predictfeas, resultfolder, dict_testfeaslabel, allmatrix)
    
    else:             
        dict_predictnameseqAAF = sware_q_AAFtrainmodel.AAF(predictpeptide, None, None, None, matrix)
        dict_predictnameseqIC50 = sware_v_IC50trainmodel.IC50(predictpeptide, None, None, None, matrix)
        dict_predictnameseqPWM = sware_y_PWMtrainmodel.PWM(predictpeptide, None, None, None, matrix)
        dict_predictnameseqWLS = sware_zc_WLStrainmodel.WLS(predictpeptide, None, None, None, matrix)
        dict_predictnameseqBSI62 = sware_zk_BSI62trainmodel.BSI62(predictpeptide, None, None, None, matrix)
        
        dict_predictfeas = OrderedDict()
        for eachname, prefea in dict_predictnameseqAAF.items():
            dict_predictfeas[eachname] = OrderedDict()  
            for eachpredictseq, fea in prefea.items():
                allfeatures = [fea, dict_predictnameseqWLS[eachname][eachpredictseq], dict_predictnameseqIC50[eachname][eachpredictseq],
                               dict_predictnameseqPWM[eachname][eachpredictseq], dict_predictnameseqBSI62[eachname][eachpredictseq]]
                dict_predictfeas[eachname][eachpredictseq] = allfeatures        
    
        return(dict_predictfeas, resultfolder)            