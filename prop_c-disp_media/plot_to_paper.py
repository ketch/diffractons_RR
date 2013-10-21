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
import imp

def plot_p(frame,path,zlimits=[0.,1.],title='',xshift=0.):
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./'+path+'/_output/_p/',file_prefix='claw_p')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    x=x+xshift
    mx=len(x); my=len(y)

    p=sol.state.q[0,:,:]
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    else:
        str_frame = "0"+str(frame)

    y_index_max = p.max(0).argmax()
    x_index_max = p[:,y_index_max].argmax()
    xm = x[x_index_max]

    pl.figure(figsize=(8,3))
    #pl.pcolormesh(xx,yy,p,cmap=cm.OrRd)
    #pl.pcolormesh(xx,yy,p,cmap=cm.RdBu)
    pl.pcolormesh(xx,yy,p,cmap=cm.Blues)
    #pl.title("t= "+str(sol.state.t),fontsize=20)
    #pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    cb = pl.colorbar();
    #pl.clim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05);
    pl.clim(0,zlimits[1]);
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    pl.xlim([xm-50, xm+15])
    pl.savefig('./_plots_to_paper/'+path+'_stress.png')
    pl.close()
    
    pl.figure(figsize=(8,3))
    pl.plot(x,p[:,3*my/4.],'-r',lw=2)
    pl.plot(x,p[:,my/4.],'--b',lw=2)
    #pl.title("t= "+str(sol.state.t),fontsize=20)
    #pl.xlabel('x',fontsize=20); pl.ylabel('Stress',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    pl.ylim([zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05]);
    pl.xlim([xm-50, xm+15])
    pl.savefig('./_plots_to_paper/'+path+'_stress_slices.eps')
    pl.close()
    
def plot_q(frame,path,zlimits=[0.,1.],title='',xshift=0.):
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./'+path+'/_output/')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    x=x+xshift
    mx=len(x); my=len(y)

    eps=sol.state.q[0,:,:]
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    else:
        str_frame = "0"+str(frame)

    y_index_max = eps.max(0).argmax()
    x_index_max = eps[:,y_index_max].argmax()
    xm = x[x_index_max]

    pl.figure(figsize=(8,3))
    #pl.pcolormesh(xx,yy,eps,cmap=cm.OrRd)
    pl.pcolormesh(xx,yy,eps,cmap=cm.Blues)
    #pl.title("t= "+str(sol.state.t),fontsize=20)
    #pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    cb = pl.colorbar();
    #pl.clim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05);
    pl.clim(0,zlimits[1]);
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    pl.xlim([xm-50, xm+15])
    pl.savefig('./_plots_to_paper/'+path+'_strain.png')
    pl.close()
    
    pl.figure(figsize=(8,3))
    pl.plot(x,eps[:,3*my/4.],'-r',lw=2)
    pl.plot(x,eps[:,my/4.],'--b',lw=2)
    #pl.title("t= "+str(sol.state.t),fontsize=20)
    #pl.xlabel('x',fontsize=20); pl.ylabel('Stress',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    pl.ylim([zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05]);
    pl.xlim([xm-50, xm+15])
    pl.savefig('./_plots_to_paper/'+path+'_strain_slices.eps')
    pl.close()
    
def get_extremum(frame,path):
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path='./'+path+'/_output/_p/',file_prefix='claw_p')
    sol_q=Solution(frame,file_format='petsc',read_aux=False,path='./'+path+'/_output/')
    
    sig=sol_p.state.q[0,:]
    eps=sol_q.state.q[0,:]
    
    return [min(np.min(sig),np.min(eps)), max(np.max(sig),np.max(eps))]

if __name__== "__main__":
    if not os.path.exists('./_plots_to_paper'): os.mkdir('_plots_to_paper')

    print('**********************')
    print('Plotting solution ...')

    frame = 375

    xshift_matchC_linear = 200
    xshift_matchC_nonlinear = 200
    xshift_miss_matchC_linear = 400
    xshift_miss_matchC_nonlinear = 400
    #Match C
    path="het_matchC_linear"
    zlimits=get_extremum(frame,path)
    plot_p(frame=frame,path=path,zlimits=zlimits,xshift=xshift_matchC_linear)
    plot_q(frame=frame,path=path,zlimits=zlimits,xshift=xshift_matchC_linear)
    path="het_matchC_nonlinear"
    zlimits=get_extremum(frame,path)
    plot_p(frame=frame,path=path,zlimits=zlimits,xshift=xshift_matchC_nonlinear)
    plot_q(frame=frame,path=path,zlimits=zlimits,xshift=xshift_matchC_nonlinear)

    #Miss-Match C
    path="het_miss-matchC_linear"
    zlimits=get_extremum(frame,path)
    plot_p(frame=frame,path=path,zlimits=zlimits,xshift=xshift_miss_matchC_linear)
    plot_q(frame=frame,path=path,zlimits=zlimits,xshift=xshift_miss_matchC_linear)
    path="het_miss-matchC_nonlinear"
    zlimits=get_extremum(frame,path)
    plot_p(frame=frame,path=path,zlimits=zlimits,xshift=xshift_miss_matchC_nonlinear)
    plot_q(frame=frame,path=path,zlimits=zlimits,xshift=xshift_miss_matchC_nonlinear)
    


