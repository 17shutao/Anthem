###check the arguments for training new model
import sys, os, re
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_j_newreadfiletrain, sware_l_newreadfilepredic, sware_m_newreadfastapredic, sware_t_checkHLAformat, sware_zp_readtestfile

def checknewargs(args):
    if not args.length:
        print('Please set the peptide length to build modle and predict, eg: --length 8')
        sys.exit(1)
    elif int(args.length) < 8 or int(args.length) > 14 or not int(args.length):
        print('The length can only be set to 8, 9, 10, 11, 12, 13, 14')
        sys.exit(1)              
    elif not args.trainingfile:
        print('Training file is missing, please specify the path of training file\neg: --trainingfile your_path/trainpeptide.txt')
        sys.exit(1)          
    elif args.predictfile_fasta and args.predictfile_peptide:
        print('Please only give prediction file with either peptide or fasta format')
        sys.exit(1)            
    elif not args.HLA:
        print('Please name your model')
        sys.exit(1)
    elif args.CV:
        if args.CV not in ['5', '10', 'LOO']:
            print('The CV can only be 5, 10, LOO')
            sys.exit(1)
    if not args.CV:
        CVtime = '5'
    else:
        CVtime = args.CV

    sware_t_checkHLAformat.checkHLAinput(args.HLA, 'trainmodel')   
    
    if args.testfile:
        testpeptide = sware_zp_readtestfile.readtestfile(args.length, args.testfile)
    else:
        testpeptide = sware_zp_readtestfile.readtestfile(args.length, None)

    [trainpeptide, neg] = sware_j_newreadfiletrain.readtrainfile(args.length, args.trainingfile)        

    if args.threshold:
        thresholdlist = list(range(50, 101))
        try:
            thresholdvalue = int(args.threshold)
        except:
            print('Please choose a integer threshold between 50 - 100')
            sys.exit(1)
        else:
            if thresholdvalue not in thresholdlist:
                print('The range of the threshold should between 50 - 100')
                sys.exit(1)

    if args.predictfile_peptide:
        predictpeptide = sware_l_newreadfilepredic.readpredictfile(args.length, args.predictfile_peptide)
        if args.threshold:
            thresholdvalueadjust = int(args.threshold)/100        
            return(trainpeptide, neg, CVtime, predictpeptide, None, testpeptide, thresholdvalueadjust)
        else:
            return(trainpeptide, neg, CVtime, predictpeptide, None, testpeptide, None)
    elif args.predictfile_fasta:
        [predictpeptide, fastafile] = sware_m_newreadfastapredic.readpredictfasta(args.length, args.predictfile_fasta)
        if args.threshold:
            thresholdvalueadjust = int(args.threshold)/100        
            return(trainpeptide, neg, CVtime, predictpeptide, fastafile, testpeptide, thresholdvalueadjust)
        else:
            return(trainpeptide, neg, CVtime, predictpeptide, fastafile, testpeptide, None)
    else:
        predictpeptide = None
        if args.threshold:
            thresholdvalueadjust = int(args.threshold)/100        
            return(trainpeptide, neg, CVtime, predictpeptide, None, testpeptide, thresholdvalueadjust)
        else:
            return(trainpeptide, neg, CVtime, predictpeptide, None, testpeptide, None)