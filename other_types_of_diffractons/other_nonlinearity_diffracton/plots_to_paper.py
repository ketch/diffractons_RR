from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import matplotlib.cm as cm
import numpy as np
import os

def plot_p(frame,file_prefix='claw_p',path='./_output/_p',plot_slices=True,plot_pcolor=True,slices_limits=None,xshift=0.0,name=''):
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    x=x+xshift
    mx=len(x); my=len(y)    
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "00"+str(frame)
    elif frame < 100:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)

    p=sol.state.q[0,:,:]

    if plot_pcolor:
        pl.figure(figsize=(8,3))
        #pl.pcolormesh(xx,yy,p,cmap=cm.OrRd)
        pl.pcolormesh(xx,yy,p,cmap=cm.Blues)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        pl.xticks([140, 150, 160, 170, 180])
        cb = pl.colorbar();
        #pl.clim(colorbar_min,colorbar_max);
        pl.clim(0,np.max(p));
        imaxes = pl.gca(); pl.axes(cb.ax)
        pl.yticks(fontsize=20); pl.axes(imaxes)
        pl.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
        pl.xlim([140,185])
        pl.savefig('./_plots_to_paper/other_nonlinearity_diffracton.png')
        #pl.show()                            
        pl.close()
    if plot_slices:
        pl.figure(figsize=(8,3))
        #pl.gcf().subplots_adjust(left=0.15)
        # plot solution of interaction
        pl.plot(x,p[:,3*my/4.],'-r',linewidth=2)
        pl.plot(x,p[:,my/4.],'-b',linewidth=2)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20)
        pl.ylabel('Stress',fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        if slices_limits is not None:
            pl.axis([slices_limits[0]+xshift,slices_limits[1]+xshift,slices_limits[2],slices_limits[3]])
        pl.xlim([140,185])
        pl.xticks([140, 150, 160, 170, 180])
        pl.savefig('./_plots_to_paper/other_nonlinearity_diffracton_slices.eps')
        pl.close()
   
if __name__== "__main__":
    if not os.path.exists('./_plots_to_paper/'): os.mkdir('./_plots_to_paper/')

    print('**********************')
    print('Plotting p ...')
    plot_p(frame=120,plot_pcolor=True,xshift=100)
    
