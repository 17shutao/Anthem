###calculate IC50 feature
from collections import OrderedDict
import json, math

def IC50(length, HLAlist, namepeptides, currentpath):    
    dict_HLAmatrix = OrderedDict()
    for eachHLA in HLAlist:
        with open('%s/feature_matrics/IC50/%s/%s.json' % (currentpath, str(length), eachHLA), 'r') as f:
            fstr = f.read()
            dict_IC50matrix = json.loads(fstr, object_pairs_hook = OrderedDict)
        dict_HLAmatrix[eachHLA] = dict_IC50matrix
        
    dict_HLAnameseqIC50 = OrderedDict()
    for eachHLA in HLAlist:
        dict_HLAnameseqIC50[eachHLA] = OrderedDict()
        for eachname, peptides in namepeptides.items():
            dict_HLAnameseqIC50[eachHLA][eachname] = OrderedDict()
            for eachpeptide in peptides:
                PSSMscore = 0
                for i in range(int(length)):
                    PSSMscore += dict_HLAmatrix[eachHLA][eachpeptide[i]][i]
                bindingscore = PSSMscore/int(length)
                IC50score = math.pow(50000, (0.8 - bindingscore)/1.6) 
                dict_HLAnameseqIC50[eachHLA][eachname][eachpeptide] = IC50score
                                   
    return dict_HLAnameseqIC50 