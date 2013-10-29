from clawpack.petclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os

def plot(frame,dirname,outname,field_name='stress',zlimits=[0.,1.],title='',xshift=0.):
    if not os.path.exists('./figures'):
        os.makedirs('./figures')

    path='./'+dirname
    if field_name=='stress':
        path=path+'/_p'
        sol=Solution(frame,file_format='petsc',path=path,read_aux=False,file_prefix='claw_p')
    else:
        sol=Solution(frame,file_format='petsc',path=path,read_aux=False)

    xx, yy = sol.state.grid.c_centers
    my=sol.state.grid.num_cells[1]
    x=sol.state.grid.x.centers

    field = sol.state.q[0,:,:]

    y_index_max = field.max(0).argmax()
    x_index_max = field[:,y_index_max].argmax()
    xm = x[x_index_max]


    plt.figure(figsize=(8,3))
    plt.pcolormesh(xx+xshift,yy,field,cmap=cm.Blues)
    plt.xticks(size=20); plt.yticks(size=20)
    cb = plt.colorbar();
    plt.clim(0,zlimits[1]);
    imaxes = plt.gca(); plt.axes(cb.ax)
    plt.yticks(fontsize=20); plt.axes(imaxes)
    plt.xlim([xm-50+xshift, xm+15+xshift])
    plt.savefig('./figures/'+outname+'_'+field_name+'.png')
    
    plt.figure(figsize=(8,3))
    plt.plot(x+xshift,field[:,3*my/4.],'-r',lw=2)
    plt.plot(x+xshift,field[:,my/4.],'--b',lw=2)
    plt.xticks(size=20); plt.yticks(size=20)
    plt.ylim([zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05]);
    plt.xlim([xm-50+xshift, xm+15+xshift])
    plt.savefig('./figures/'+outname+'_'+field_name+'_slices.eps')
    
    
def get_extremum(frame,dirname):
    path='./'+dirname
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path=path+'/_p',file_prefix='claw_p')
    sol_q=Solution(frame,file_format='petsc',read_aux=False,path=path)
    
    sig=sol_p.state.q[0,:]
    eps=sol_q.state.q[0,:]
    
    return [min(np.min(sig),np.min(eps)), max(np.max(sig),np.max(eps))]
