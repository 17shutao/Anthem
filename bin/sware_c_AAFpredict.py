###calculate AAF feature
from collections import OrderedDict
import json, os

def AAF(length, HLAlist, namepeptides, currentpath):
    dict_HLAmatrix = OrderedDict()
    for eachHLA in HLAlist:
        with open('%s/feature_matrics/AAF/%s/%s.json' % (currentpath, str(length), eachHLA), 'r') as f:
            fstr = f.read()
            dict_AAFmatrix = json.loads(fstr, object_pairs_hook = OrderedDict)        
        dict_HLAmatrix[eachHLA] = dict_AAFmatrix
    
    dict_HLAnameseqAAF = OrderedDict()
    for eachHLA in HLAlist:
        dict_HLAnameseqAAF[eachHLA] = OrderedDict()
        for eachname, peptides in namepeptides.items():
            dict_HLAnameseqAAF[eachHLA][eachname] = OrderedDict()
            for eachpeptide in peptides:
                AAFscore = 0
                for i in range(int(length)):
                    AAFscore += dict_HLAmatrix[eachHLA][eachpeptide[i]][i]
                dict_HLAnameseqAAF[eachHLA][eachname][eachpeptide] = AAFscore
      
    return dict_HLAnameseqAAF
    