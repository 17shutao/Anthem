###calculate WLS feature
from collections import OrderedDict
import json

def WLS(length, HLAlist, namepeptides, currentpath):
    dict_HLAmatrix = OrderedDict()
    for eachHLA in HLAlist:
        with open('%s/feature_matrics/WLS/%s/%s.json' % (currentpath, str(length), eachHLA), 'r') as f:
            fstr = f.read()
            dict_WLSmatrix = json.loads(fstr, object_pairs_hook = OrderedDict)
        dict_HLAmatrix[eachHLA] = dict_WLSmatrix
                
    dict_HLAnameseqWLS = OrderedDict()
    for eachHLA in HLAlist:
        dict_HLAnameseqWLS[eachHLA] = OrderedDict()     
        for eachname, peptides in namepeptides.items():
            dict_HLAnameseqWLS[eachHLA][eachname] = OrderedDict()
            for eachpeptide in peptides:
                WLSscore = 0
                for i in range(int(length)):
                    WLSscore += float(dict_HLAmatrix[eachHLA][eachpeptide[i]][i])
                dict_HLAnameseqWLS[eachHLA][eachname][eachpeptide] = WLSscore        
    
    return dict_HLAnameseqWLS