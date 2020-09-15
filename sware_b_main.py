#main function
#!/usr/bin/env python
#_*_coding:utf-8_*_
import argparse
from bin import *

parser = argparse.ArgumentParser(usage="it's usage tip.")
parser.add_argument('--length', type=int, help='the length of the peptide, from 8 - 14')
parser.add_argument('--threshold', type=int, help='the threshold to define predicted binder, integer from 50 - 100')
parser.add_argument('--peptide_file', type=str, help='the path of the file contains peptides')
parser.add_argument('--fasta_file', type=str, help='the path of the fasta file contains sequence')
parser.add_argument('--HLA', type=str, help='the HLA allele you want to predict peptides bind to, eg: HLA-A*01:01 or HLA-A*01:01,HLA-A*02:01' + '\n'
                    'Or if you want to build a new model, please name it. Note: you can only choose letters from A-Z or a-z or numbers 0-9. There should no space in the model name')
parser.add_argument('--trainingfile', type=str, help='the path of training file to train model')
parser.add_argument('--CV', type=str, help='set the number of folds for cross-validation, either choose 5, 10 or LOO(leave one out), the default is 5')
parser.add_argument('--testfile', type=str, help='the path of test file to test model')
parser.add_argument('--predictfile_peptide', type=str, help='the path of predict file in peptide format')
parser.add_argument('--predictfile_fasta', type=str, help='the path of predict file in fasta format')
parser.add_argument('--mode', type=str, help='choose the mode, three modes can be used: prediction, TrainYourModel, useYourOwnModel, only select one mode each time')
parser.add_argument('--modelfile', type=str, help='the file named "modelfile" generated from the TrainYourModel mode')
args = parser.parse_args()

if args.mode == 'prediction':
    ###check input
    [HLAlist, namepeptides, peptidelength, fastafile, thresholdvalueadjust] = sware_h_checkargs.checkarg(args)

    if peptidelength == None:
        ###generate features
        featurecombine = sware_f_combinefeature.combfea(str(args.length), HLAlist, namepeptides)
        ###use model
        sware_d_prediction.pred(featurecombine, str(args.length), fastafile, thresholdvalueadjust)
    else:
        ###generate features
        featurecombine = sware_f_combinefeature.combfea(str(peptidelength), HLAlist, namepeptides)        
        ###wirte in arff and predict
        sware_d_prediction.pred(featurecombine, str(peptidelength), fastafile, thresholdvalueadjust)

elif args.mode == 'TrainYourModel':
    ###check training and predict peptide
    [trainpeptide, neg, CVtime, predictpeptide, fastafile, testpeptide, thresholdvalueadjust] = sware_n_checknewargs.checknewargs(args)
    
    if neg == 'NO':
        ###generate neg
        trainpeptide = sware_o_negpeptide.generateneg(args.length, trainpeptide)
        ###generate features
        [train_feaslabel, predic_feas, resultfolder, test_feaslabel, allmatrix] = sware_p_combinefeaturefornew.combfeanew(str(args.length), predictpeptide, trainpeptide, args.HLA, None, testpeptide)
        ##train model
        sware_r_trainmodel.trainmodel(predic_feas, args.length, resultfolder, CVtime, train_feaslabel, args.HLA, None, fastafile, test_feaslabel, thresholdvalueadjust, allmatrix, None)
    else:
        ###generate features
        [train_feaslabel, predic_feas, resultfolder, test_feaslabel, allmatrix] = sware_p_combinefeaturefornew.combfeanew(str(args.length), predictpeptide, trainpeptide, args.HLA, None, testpeptide)
        ###train model
        sware_r_trainmodel.trainmodel(predic_feas, args.length, resultfolder, CVtime, train_feaslabel, args.HLA, None, fastafile, test_feaslabel, thresholdvalueadjust, allmatrix, 'negprovided')

elif args.mode == 'useYourOwnModel':
    ###check predict peptide
    [predictpeptide, HLA, lengthmatrix, fastafile, thresholdvalueadjust, allmatrix] = sware_s_checkusemodelargs.chaeckusemodelargs(args)

    ###generate features for predict peptide
    [predic_feas, resultfolder] = sware_p_combinefeaturefornew.combfeanew(lengthmatrix, predictpeptide, None, None, allmatrix, None)    
    
    ###use model to predict
    sware_r_trainmodel.trainmodel(predic_feas, lengthmatrix, resultfolder, None, None, HLA, args.modelfile, fastafile, None, thresholdvalueadjust, allmatrix, None)         
elif args.mode != 'prediction' or args.mode != 'TrainYourModel' or args.mode != 'useYourOwnModel':
    print('Please at least select one mode to use the tool!' + '\n'
          + 'You can choose either "prediction", "TrainYourModel" or "useYourOwnModel"')
    