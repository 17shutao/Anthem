###check the arguments for using the trained model
import sys, os, re, json, tarfile
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_l_newreadfilepredic, sware_m_newreadfastapredic

def chaeckusemodelargs(args):
    if not args.predictfile_fasta and not args.predictfile_peptide:
        print('prediction file is missing, please specify the path of prediction file\neg: --predictfile_peptide path/predictpeptide.txt')
        sys.exit(1)        
    elif args.predictfile_fasta and args.predictfile_peptide:
        print('Please only give prediction file with either peptide or fasta format')
        sys.exit(1)
    if not args.modelfile:
        print('The modelfile is missing, please provide it\neg: --modelfile your_path/modelfile')
        sys.exit(1)

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

    ###get length of matrix
    with tarfile.open(args.modelfile, 'r') as tar:
        tar.extractall('model')
    with open('model/allmatrix.json', 'r') as f:
        allmatrix = json.load(f)
        
    for HLA, matrixs in allmatrix.items():
        dict_AAlenfre = matrixs['AAF']        
        
    lengthmatrix = len(dict_AAlenfre['A']) 

    if args.predictfile_peptide:
        [predictpeptide, HLA] = sware_l_newreadfilepredic.readpredictfile(lengthmatrix, args.predictfile_peptide, allmatrix)
        if args.threshold:
            thresholdvalueadjust = int(args.threshold)/100        
            return(predictpeptide, HLA, lengthmatrix, None, thresholdvalueadjust, allmatrix)
        else:
            return(predictpeptide, HLA, lengthmatrix, None, None, allmatrix)
    elif args.predictfile_fasta:
        [predictpeptide, HLA, fastafile] = sware_m_newreadfastapredic.readpredictfasta(lengthmatrix, args.predictfile_fasta, allmatrix)
        if args.threshold:
            thresholdvalueadjust = int(args.threshold)/100        
            return(predictpeptide, HLA, lengthmatrix, fastafile, thresholdvalueadjust, allmatrix)
        else:
            return(predictpeptide, HLA, lengthmatrix, fastafile, None, allmatrix)