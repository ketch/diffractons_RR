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

def plot_p(frame,zlim=[0.,1.]):
    sol_ref=Solution(frame+450,file_format='petsc',read_aux=False,path='_output/reference/_p/',file_prefix='claw_p')
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers

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

    # plot pcolor
    pl.figure(figsize=(8,3))
    pl.pcolormesh(xx,yy,p,cmap=cm.OrRd)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    cb = pl.colorbar();
    pl.clim(zlim[0]-0.05*zlim[1],1.05*zlim[1]);
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    pl.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
    pl.savefig('./_plots_to_animation/pcolor/co-interaction_'+str_frame+'.png')
    pl.close()
    # Plot slices
    pl.figure(figsize=(8,3))
    pl.gcf().subplots_adjust(left=0.10)
    # plot reference
    pl.plot(x,p_ref[:,3*my/4.],'--r',linewidth=1)
    pl.plot(x,p_ref[:,my/4.],'--b',linewidth=1)
    # plot solution of interaction
    pl.plot(x,p[:,3*my/4.],'-r',linewidth=2)
    pl.plot(x,p[:,my/4.],'-b',linewidth=2)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    #pl.xlabel('x',fontsize=20)
    #pl.ylabel('Stress',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    pl.ylim(zlim[0]-0.05*zlim[1],1.05*zlim[1]);
    pl.savefig('./_plots_to_animation/slices/co-interaction_'+str_frame+'_slices.eps')
    pl.close()

def get_extremum_p(frame):
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    sig=sol_p.state.q[0,:]
    return [np.min(sig), np.max(sig)]
   
if __name__== "__main__":
    if not os.path.exists('./_plots_to_animation'): os.mkdir('_plots_to_animation')
    if not os.path.exists('./_plots_to_animation/pcolor'): os.mkdir('_plots_to_animation/pcolor')
    if not os.path.exists('./_plots_to_animation/slices'): os.mkdir('_plots_to_animation/slices')

    from_frame = 0
    to_frame   = 500

    zlim = get_extremum_p(450)

    print('**********************')
    print('Plotting p ...')
    for i in xrange(from_frame,to_frame+1):
        plot_p(frame=i,zlim=zlim)
        print ('frame '+str(i)+' plotted')

