import os
installdir=os.environ['CHESS_ROOT']
CHESS_TOOLBOX=os.path.join(installdir,'chess_toolbox')

def densify(file_in,file_out,mpirun=None):
    command=CHESS_TOOLBOX+' convert-matrix-format bigdft_to_dense '+\
             file_in+' '+file_out
    if mpirun: command=mpirun+' '+command
    os.system(command)

def reduce_dense_matrix(lut,f):
    """Create a new density matrix out of the dense matrix indicated by file f"""
    ff=open(f,'r')
    ff.readline() #read first line
    fw=[] #open(f+'-reduced','w')
    #fw.write('# 1 '+str(len(lut))+'\n')
    for l in ff.readlines():
        entry=l.split()
        i,j,val=entry[0:3]
        i=int(i)-1
        j=int(j)-1
        if i in lut and j in lut:
            inew=lut.index(i)+1
            jnew=lut.index(j)+1
            #fw.write(str(inew)+' '+str(jnew)+' '+str(val)+'\n')
            fw.append((inew,jnew,val))
            #print 'dumpof',i,j,'into',str(i/4+1),str(j/4+1)
    ff.close()
    return fw #fw.close()

def get_atomic_lookup(metadata_file):
    alut=[]
    f=open(metadata_file,'r')
    for l in f.readlines():
        if 'on_which_atom' in l:
            entry=l.split()
            alut.append(int(entry[0])-1)
    f.close()
    return alut
