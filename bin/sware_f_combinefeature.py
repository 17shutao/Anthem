###combine features
from collections import OrderedDict
import sys, os, re
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_c_AAFpredict, sware_u_IC50predict, sware_x_PWMpredict, sware_zb_WLSfpredict, sware_zd_BSI62predict

def combfea(length, HLAlist, peptides):
    currentpath = os.getcwd()
    AAFfeature = sware_c_AAFpredict.AAF(length, HLAlist, peptides, currentpath)
    IC50feature = sware_u_IC50predict.IC50(length, HLAlist, peptides, currentpath)
    PWMfeature = sware_x_PWMpredict.PWM(length, HLAlist, peptides, currentpath)
    WLSfeature = sware_zb_WLSfpredict.WLS(length, HLAlist, peptides, currentpath)
    BSI62feature = sware_zd_BSI62predict.BSI62(length, HLAlist, peptides, currentpath)
    
    allfeaturecomb = OrderedDict()
    for eachHLA, names in AAFfeature.items():
        allfeaturecomb[eachHLA] = OrderedDict()
        for eachname, seqfea in names.items():
            allfeaturecomb[eachHLA][eachname] = OrderedDict()
            for eachseq, fea in seqfea.items():
                allfeatures = [fea, WLSfeature[eachHLA][eachname][eachseq],
                               IC50feature[eachHLA][eachname][eachseq],
                               PWMfeature[eachHLA][eachname][eachseq],
                               BSI62feature[eachHLA][eachname][eachseq]]
                allfeaturecomb[eachHLA][eachname][eachseq] = allfeatures
            
    return allfeaturecomb