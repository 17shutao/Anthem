###generate seqlogo
import os

def seqlogo(templatefolder, resultfolder, predictfunctionlogo, trainfunctionlogo):
    currentpath = os.getcwd()
    if trainfunctionlogo == None:
        for eachHLA, names in predictfunctionlogo.items():
            for eachname, peptidelist in names.items():
                if len(peptidelist) != 0:
                    rownum = 1
                    fastaname1 = eachHLA.split(':')[0]
                    fastaname2 = eachHLA.split(':')[1]
                    os.system('cp %s/source/logopng/%s-001.png %s/%s' % (currentpath, eachHLA, currentpath, resultfolder))
                    with open('%s/%s/%s_%s_%s.fasta' % (currentpath, templatefolder, fastaname1, fastaname2, eachname), 'w') as f:
                        for eachpep in peptidelist:
                            f.write('>' + str(rownum) + '\n')
                            f.write(eachpep + '\n')
                            rownum += 1
        
        fastalist = os.listdir('%s/%s' % (currentpath, templatefolder))
        os.chdir('%s/%s' % (currentpath, templatefolder))
        for eachfasta in fastalist:
            if '.fasta' in eachfasta:
                filename = eachfasta.split('.')[0]
                os.system('python %s/source/seq2logo/Seq2Logo.py -S 1 -Z on --format PNG -f %s -o %s/%s/%s' % (currentpath, eachfasta, currentpath, resultfolder, filename))
        
    elif predictfunctionlogo == None:

        for eachname, peptidelist in trainfunctionlogo.items():
            if len(peptidelist) != 0:
                rownum = 1
                with open('%s/%s/%s.fasta' % (currentpath, templatefolder, eachname), 'w') as f:
                    for eachpep in peptidelist:
                        f.write('>' + str(rownum) + '\n')
                        f.write(eachpep + '\n')
                        rownum += 1
         
        fastalist = os.listdir('%s/%s' % (currentpath, templatefolder))
        os.chdir('%s/%s' % (currentpath, templatefolder))
        for eachfasta in fastalist:
            if '.fasta' in eachfasta:
                filename = eachfasta.split('.')[0]
                os.system('python %s/source/seq2logo/Seq2Logo.py -S 1 -Z on --format PNG -f %s -o %s/%s/%s' % (currentpath, eachfasta, currentpath, resultfolder, filename))