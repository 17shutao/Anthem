###check the arguments for predict
import sys, os, re
pPath = re.sub(r'bin$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)
from bin import sware_e_readfile, sware_g_readfasta

def checkarg(args):
    if not args.peptide_file and not args.fasta_file:
        print('Please give the peptide or sequence file to predict,\neg: --peptide_file your_path/peptide.txt or --fasta_file your_path/sequence.txt')
        sys.exit(1)        
    elif args.peptide_file and args.fasta_file:
        print('Please only give one file either in peptide sequence or in fasta format')
        sys.exit(1)            
    elif not args.HLA:
        print('Please set the HLA allele to predict, eg: --HLA HLA-A*01:01.\nIf multiple alleles, use comma to separate without space, eg: HLA-A*01:01,HLA-A*02:01')
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
                                    
    if args.peptide_file:
        [HLAlist, namepeptides, peptidelength] = sware_e_readfile.readfile(args.HLA, args.peptide_file)
        if args.threshold:
            thresholdvalueadjust = int(args.threshold)/100        
            return([HLAlist, namepeptides, peptidelength, None, thresholdvalueadjust])
        else:
            return([HLAlist, namepeptides, peptidelength, None, None])
    elif args.fasta_file:
        if not args.length:
            print('Please set the peptide length to predict, eg: --length 8')
            sys.exit(1)
        elif int(args.length) < 8 or int(args.length) > 14 or not int(args.length):
            print('The length can only be set to 8, 9, 10, 11, 12, 13, 14')
            sys.exit(1)              
        [HLAlist, namepeptides, fastafile] = sware_g_readfasta.readfasta(str(args.length), args.HLA, args.fasta_file)
        if args.threshold:
            thresholdvalueadjust = int(args.threshold)/100
            return([HLAlist, namepeptides, None, fastafile, thresholdvalueadjust])
        else:
            return([HLAlist, namepeptides, None, fastafile, None])