###check HLA input format
import sys

def checkHLAinput(HLAtype, trainmodel):
    if trainmodel == None:
        HLAlist = []
        if ',' not in HLAtype:
            if 'HLA-' not in HLAtype or '*' not in HLAtype or ':' not in HLAtype:
                print('Please provide correct HLA allele name, eg: HLA-A*01:01')
                sys.exit(1)
            elif HLAtype.count('HLA') > 1:
                print('If multiple HLA alleles, use comma to separate and without space, eg: HLA-A*01:01,HLA-A*01:01')
                sys.exit(1)
            elif len(HLAtype) != 11:
                print('Please provide correct HLA allele name, eg: HLA-A*01:01')
                sys.exit(1)            
            HLAlist.append(HLAtype.strip())            
        elif ',' in HLAtype:                        
            if HLAtype.count('HLA') < 2:
                print('If multiple HLA alleles, use comma to separate and without space, eg: HLA-A*01:01,HLA-A*01:01')
                sys.exit(1)
            HLAtype = HLAtype.split(',')
            for eachHLA in HLAtype:
                if 'HLA-' not in eachHLA or '*' not in eachHLA or ':' not in eachHLA:
                    print('Please provide correct HLA allele name, eg: HLA-A*01:01')
                    sys.exit(1)
                elif len(eachHLA) != 11:
                    print('Please provide correct HLA allele name, eg: HLA-A*01:01')
                    sys.exit(1)                   
                HLAlist.append(eachHLA.strip())       
        
        return(HLAlist)
    else:
        namellist = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                     'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
                     '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
        if ' ' in HLAtype:
            print("There should be no space in the model name")
            sys.exit(1)
        elif len(HLAtype) > 20:
            print("The length of the model name should be less than 20 characters")
            sys.exit(1)
        else:
            for eachone in HLAtype:
                if eachone not in namellist:
                    print('You can only choose letters from Aa-Zz or 0-9')
                    sys.exit(1)    