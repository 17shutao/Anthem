###calculate PWM feature
from collections import OrderedDict
import json

def PWM(length, HLAlist, namepeptides, currentpath): 
    dict_HLAmatrix = OrderedDict()
    for eachHLA in HLAlist:
        with open('%s/feature_matrics/PWM/%s/%s.json' % (currentpath, str(length), eachHLA), 'r') as f:
            fstr = f.read()
            dict_PWMmatrix = json.loads(fstr, object_pairs_hook = OrderedDict)
        dict_HLAmatrix[eachHLA] = dict_PWMmatrix
                                    
    dict_HLAnameseqPWM = OrderedDict()
    for eachHLA in HLAlist:
        dict_HLAnameseqPWM[eachHLA] = OrderedDict()
        for eachname, peptides in namepeptides.items():
            dict_HLAnameseqPWM[eachHLA][eachname] = OrderedDict()
            for eachpeptide in peptides:
                PWMscore = 0
                for i in range(int(length)):
                    PWMscore += dict_HLAmatrix[eachHLA][eachpeptide[i]][i]
                dict_HLAnameseqPWM[eachHLA][eachname][eachpeptide] = PWMscore
                 
    return dict_HLAnameseqPWM