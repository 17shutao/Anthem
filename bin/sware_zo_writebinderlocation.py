###write the binder location when fasta format
import os

def writebinder(resultfolder, length, newHLA, dict_HLAnameposlevel = None, dict_nameposlevel = None):
    currentpath = os.getcwd()
    if dict_HLAnameposlevel != None:
        with open('%s/%s/length_%s_binderlocinfasta.txt' % (currentpath, resultfolder, length), 'w') as f:
            numHLA = 1                
            for eachHLA, names in  dict_HLAnameposlevel.items():
                numname = 1
                f.write(str(len(names.keys())) + '\t' + str(length) + '\t' + eachHLA + '\n')
                for eachname, poslevellist in names.items():
                    if len(poslevellist) != 0:
                        for eachtuple in poslevellist:
                            f.write(str(eachtuple[0]) + '\t' + eachtuple[1] + '\t' + eachtuple[2] + '\t')
                        if numHLA == len(dict_HLAnameposlevel.keys()) and numname == len(names.keys()):
                            continue
                        else:
                            f.write('\n')
                    else:
                        if numHLA == len(dict_HLAnameposlevel.keys()) and numname == len(names.keys()):
                            continue
                        else:
                            f.write('\n')
                    numname += 1
                numHLA += 1
    
    elif dict_nameposlevel != None:
        with open('%s/%s/usermodel_length_%s_binderlocinfasta.txt' % (currentpath, resultfolder, length), 'w') as f:
            numname = 1
            f.write(str(len(dict_nameposlevel.keys())) + '\t' + str(length) + '\t' + newHLA + '\n')
            for eachname, poslevellist in dict_nameposlevel.items():
                if len(poslevellist) != 0:
                    for eachtuple in poslevellist:
                        f.write(str(eachtuple[0]) + '\t' + eachtuple[1] + '\t' + eachtuple[2] + '\t')
                    if numname == len(dict_nameposlevel.keys()):
                        continue
                    else:
                        f.write('\n')
                else:
                    if numname == len(dict_nameposlevel.keys()):
                        continue
                    else:
                        f.write('\n')
                numname += 1          