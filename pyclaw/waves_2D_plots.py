from clawpack.petclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
import matplotlib.cm as cm
import numpy as np
import os

def plot(frame,dirname,clim=None,axis_limits=None):
    if not os.path.exists('./figures'):
        os.makedirs('./figures')
        
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./saved_data/'+dirname+'/_p/',file_prefix='claw_p')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)
    
    mp=sol.state.num_eqn    
    yy,xx = np.meshgrid(y,x)

    p=sol.state.q[0,:,:]
    if clim is not None:
        pl.pcolormesh(xx,yy,p,cmap=cm.RdBu_r)
    else:
        pl.pcolormesh(xx,yy,p,cmap=cm.Reds)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    cb = pl.colorbar();

    if clim is not None:
        pl.clim(clim[0],clim[1]);
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    pl.axis('equal')
    if axis_limits is None:
        pl.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
    else:
        pl.axis([axis_limits[0],axis_limits[1],axis_limits[2],axis_limits[3]])
    pl.savefig('./figures/'+dirname+'.png')
    pl.close()
           
if __name__== "__main__":
    plot(65,'figure_14b',clim=[-0.3,0.3])
    plot(0,'figure_14a',axis_limits=[0.25, 10.25, 0.25, 10.25])
    plot(65,'figure_14c',clim=[-0.3,0.3])
    plot(65,'figure_14d',clim=[-0.3,0.3])
    
