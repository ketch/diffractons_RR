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

def plot_p(frame,file_prefix='claw_p',path='./_output/_p',plot_slices=True,plot_pcolor=True,slices_limits=None,xshift=0.0,name='',title=True):
    sol_ref=Solution(frame+450,file_format='petsc',read_aux=False,path='_output/reference/_p/',file_prefix=file_prefix)
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
    p_ref=sol_ref.state.q[0,:,:]

    if plot_pcolor:
        pl.pcolormesh(xx,yy,p,cmap=cm.OrRd)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        cb = pl.colorbar();
        #pl.clim(colorbar_min,colorbar_max);
        imaxes = pl.gca(); pl.axes(cb.ax)
        pl.yticks(fontsize=20); pl.axes(imaxes)
        pl.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
        #pl.axis([0.25,60.25,0.25,60.25])
        pl.savefig('./_plots_to_paper/co-interaction_'+str_frame+name+'.png')
        #pl.show()                            
        pl.close()
    if plot_slices:
        pl.figure(figsize=(8,3))
        pl.gcf().subplots_adjust(left=0.10)
        # plot reference
        pl.plot(x,p_ref[:,my/4.],'--b',linewidth=1)
        pl.plot(x,p_ref[:,3*my/4.],'--r',linewidth=1)
        # plot solution of interaction
        pl.plot(x,p[:,3*my/4.],'-r',linewidth=2)
        pl.plot(x,p[:,my/4.],'-b',linewidth=2)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20)
        if title:
            pl.ylabel('Stress',fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        if slices_limits is not None:
            pl.axis([slices_limits[0]+xshift,slices_limits[1]+xshift,slices_limits[2],slices_limits[3]])
        pl.savefig('./_plots_to_paper/co-interaction_'+str_frame+name+'.eps')
        pl.close()

def get_extremum_p(frame):
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    sig=sol_p.state.q[0,:]
    return [np.min(sig), np.max(sig)]
   
if __name__== "__main__":
    if not os.path.exists('./_plots_to_paper'): os.mkdir('_plots_to_paper')   

    ylim = get_extremum_p(450)

    print('**********************')
    print('Plotting p ...')
    plot_p(frame=0,plot_pcolor=False,slices_limits=[15,45,ylim[0]-0.05*ylim[1],1.05*ylim[1]],xshift=0)
    plot_p(frame=200,plot_pcolor=False,slices_limits=[75,115,ylim[0]-0.05*ylim[1],1.05*ylim[1]],xshift=200)
    plot_p(frame=450,plot_pcolor=False,slices_limits=[20,50,ylim[0]-0.05*ylim[1],1.05*ylim[1]],xshift=600)
    plot_p(frame=450,plot_pcolor=False,slices_limits=[27.5,42.5,-0.01,0.1],xshift=600,name='_zoom',title=False)
