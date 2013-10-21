from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
from clawpack.petclaw.io.petsc import write_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
#from matplotlib import rc
#rc('text', usetex=True)
import numpy as np

def get_max(frame,name):
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path+'_p/',file_prefix=file_prefix)
    sol_q=Solution(frame,file_format='petsc',read_aux=False,path=path)

    sig=sol.state.q[0,:,:]
    u=sol_q.state.q[1,:,:]

    y_index_max = sig.max(0).argmax()
    x_index_max =sig[:,y_index_max].argmax()

    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    my=len(y)

    file=open(name+'.txt','a')
    file.write(str(sol.state.t)+' '+str(x[x_index_max])+' '+str(sig[x_index_max,y_index_max])+' '+str(u[x_index_max,y_index_max])+'\n')
    file.close()

def record_max(from_frame,to_frame,name):
    t=np.zeros(to_frame-from_frame+1)
    x=np.zeros(to_frame-from_frame+1)
    file=open(name+'.txt','w')
    file.close()
    print '     ',
    for i in xrange(from_frame,to_frame+1):
        print str(i),
        get_max(i,name)
    print ''

if __name__== "__main__":    
    from_frame = 1470
    to_frame = 1570

    file_prefix='claw_p'

    print 'recording max points for sw1'
    path='./_output/sw1/'
    #record_max(from_frame,to_frame,'sw1')

    print 'recording max points for sw2'
    path='./_output/sw2/'
    #record_max(from_frame,to_frame,'sw2')

    print 'recording max points for sw3'
    path='./_output/sw3/'
    #record_max(from_frame,to_frame,'sw3')

    print 'recording max points for sw4'
    path='./_output/sw4/'
    #record_max(from_frame,to_frame,'sw4')

    print 'recording max points for sw5'
    path='./_output/sw5/'
    #record_max(from_frame,to_frame,'sw5')

    print 'recording max points for sw6'
    path='./_output/sw6/'
    #record_max(from_frame,to_frame,'sw6')

    from_frame = 1970
    to_frame = 2070
    
    print 'recording max points for sw7'
    path='./_output/sw7/'
    #record_max(from_frame,to_frame,'sw7')

    print 'recording max points for sw8'
    path='./_output/sw8/'
    record_max(from_frame,to_frame,'sw8')


