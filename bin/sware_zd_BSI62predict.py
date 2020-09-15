###calculate BSI62 feature
from collections import OrderedDict
import json, os

def BSI62(length, HLAlist, namepeptides, currentpath):
    dict_HLAmatrix = OrderedDict()
    for eachHLA in HLAlist:
        with open('%s/feature_matrics/BSI62/%s/%s.json' % (currentpath, str(length), eachHLA), 'r') as f:
            fstr = f.read()
            dict_AAFmatrix = json.loads(fstr, object_pairs_hook = OrderedDict)        
        dict_HLAmatrix[eachHLA] = dict_AAFmatrix
    
    dict_HLAnameseqBSI62 = OrderedDict()
    for eachHLA in HLAlist:
        dict_HLAnameseqBSI62[eachHLA] = OrderedDict()
        for eachname, peptides in namepeptides.items():
            dict_HLAnameseqBSI62[eachHLA][eachname] = OrderedDict()
            for eachpeptide in peptides:
                AAFscore = 0
                for i in range(int(length)):
                    AAFscore += dict_HLAmatrix[eachHLA][eachpeptide[i]][i]
                dict_HLAnameseqBSI62[eachHLA][eachname][eachpeptide] = AAFscore
      
    return dict_HLAnameseqBSI62
    