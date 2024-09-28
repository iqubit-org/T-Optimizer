import sys

#print(len(sys.argv))
#print(sys.argv[0])
#print(sys.argv[1])
if len(sys.argv)!=2:
    print('python3 tfc2qc.py filename.tfc')
else:
    with open (sys.argv[1], "r") as fin, open(sys.argv[1].replace('.tfc','.qc'),'w') as fout:
        for line in fin:
            oline=line.replace(',',' ')
            osplit=oline.split(' ')
            if osplit[0]=='t3':
                fout.write('H '+osplit[3])
                fout.write(oline.replace('t3 ','Z '))
                fout.write('H '+osplit[3])
            elif osplit[0]=='t2':
                fout.write(oline.replace('t2 ', 'CX '))
            else:
                fout.write(oline)
