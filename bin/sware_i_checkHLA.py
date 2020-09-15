###check whether HLA is available in the selected length
import json, sys, os

def checkHLA(length, HLAlist):
    currentpath = os.getcwd()
    with open('%s/source/lenghHLA.json' % currentpath) as f:
        dict_lenHLA = json.load(f)
        lenghHLAavailable = dict_lenHLA[length]
        
    for eachHLA in HLAlist:
        if eachHLA not in lenghHLAavailable:
            print('Currently, we can not predict peptide bind to %s in the length you choose, please select another' % eachHLA +
                  ' allotype or train a new model of that HLA of the length you choose.\nFor choosing HLA that can be predict, please refer to the HLA list file.' + '\n' + 
                  'For training new model, please see README file.')
            sys.exit(1)
    
    return()
                