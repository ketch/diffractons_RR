from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
import matplotlib.cm as cm
#rc('text', usetex=True)
import numpy as np
import os

def plot_p(frame,file_prefix='claw_p',path='./_output/_p',plot_slices=True,plot_pcolor=True,zlimits=[0.,1.]):
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
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

    if plot_pcolor:
        pl.figure(figsize=(8,3))
        pl.pcolormesh(xx,yy,p,cmap=cm.OrRd)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        cb = pl.colorbar();
        pl.clim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05);
        imaxes = pl.gca(); pl.axes(cb.ax)
        pl.yticks(fontsize=20); pl.axes(imaxes)
        pl.savefig('./_plots/no_track/stress/pcolor/stress'+str_frame+'.png')
        pl.xlim([xm-60, xm+20])
        pl.savefig('./_plots/track/stress/pcolor/stress'+str_frame+'.png')
        pl.close()
    if plot_slices:
        pl.figure(figsize=(8,3))
        pl.plot(x,p[:,3*my/4.],'-r',lw=2)
        pl.plot(x,p[:,my/4.],'--b',lw=2)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20); pl.ylabel('Stress',fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        pl.ylim([zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05]);
        pl.savefig('./_plots/no_track/stress/slices/stress'+str_frame+'_slices.eps')
        pl.xlim([xm-60, xm+20])
        pl.savefig('./_plots/track/stress/slices/stress'+str_frame+'_slices.eps')
        pl.close()

def plot_q(frame,file_prefix='claw',path='./_output/',plot_pcolor=True,plot_slices=True,zlimits=[0,1]):
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
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
    x_index_max =eps[:,y_index_max].argmax()
    xm = x[x_index_max]

    if plot_pcolor:
        pl.figure(figsize=(8,3))
        pl.pcolormesh(xx,yy,eps,cmap=cm.OrRd)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20); pl.ylabel('y',fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        cb = pl.colorbar(); 
        pl.clim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05);
        imaxes = pl.gca(); pl.axes(cb.ax)
        pl.yticks(fontsize=20); pl.axes(imaxes)
        pl.savefig('./_plots/no_track/strain/pcolor/strain'+str_frame+'.png')
        pl.xlim([xm-60, xm+20])
        pl.savefig('./_plots/track/strain/pcolor/strain'+str_frame+'.png')
        pl.close()

    if plot_slices:
        pl.figure(figsize=(8,3))
        pl.plot(x,eps[:,3*my/4.],'-r',lw=2)
        pl.plot(x,eps[:,my/4.],'-b',lw=2)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20)
        pl.ylabel("Strain",fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        pl.ylim([zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05]);
        pl.savefig('./_plots/no_track/strain/slices/strain'+str_frame+'_slices.eps')
        pl.xlim([xm-60, xm+20])
        pl.savefig('./_plots/track/strain/slices/strain'+str_frame+'_slices.eps')
        pl.close()

def get_max(frame):
    sol_q=Solution(frame,file_format='petsc',read_aux=False)
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p',file_prefix='claw_p')
    eps=sol_q.state.q[0,:,:]
    sig=sol_p.state.q[0,:,:]
    return np.max(np.max(eps)), np.max(np.max(sig))

def get_min(frame):
    sol_q=Solution(frame,file_format='petsc',read_aux=False)
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p',file_prefix='claw_p')
    eps=sol_q.state.q[0,:,:]
    sig=sol_p.state.q[0,:,:]
    return np.min(np.min(eps)), np.min(np.min(sig))

if __name__== "__main__":
    if not os.path.exists('./_plots'): os.mkdir('_plots')
    if not os.path.exists('./_plots/track'): os.mkdir('./_plots/track')
    if not os.path.exists('./_plots/track/stress'): os.mkdir('./_plots/track/stress')
    if not os.path.exists('./_plots/track/stress/pcolor'): os.mkdir('./_plots/track/stress/pcolor')
    if not os.path.exists('./_plots/track/stress/slices'): os.mkdir('./_plots/track/stress/slices')
    if not os.path.exists('./_plots/track/strain'): os.mkdir('./_plots/track/strain')
    if not os.path.exists('./_plots/track/strain/pcolor'): os.mkdir('./_plots/track/strain/pcolor')
    if not os.path.exists('./_plots/track/strain/slices'): os.mkdir('./_plots/track/strain/slices')
    if not os.path.exists('./_plots/no_track'): os.mkdir('./_plots/no_track')
    if not os.path.exists('./_plots/no_track/stress'): os.mkdir('./_plots/no_track/stress')
    if not os.path.exists('./_plots/no_track/stress/pcolor'): os.mkdir('./_plots/no_track/stress/pcolor')
    if not os.path.exists('./_plots/no_track/stress/slices'): os.mkdir('./_plots/no_track/stress/slices')
    if not os.path.exists('./_plots/no_track/strain'): os.mkdir('./_plots/no_track/strain')
    if not os.path.exists('./_plots/no_track/strain/pcolor'): os.mkdir('./_plots/no_track/strain/pcolor')
    if not os.path.exists('./_plots/no_track/strain/slices'): os.mkdir('./_plots/no_track/strain/slices')

    from_frame = 0
    to_frame   = 400

    eps_max, sig_max = get_max(10)
    eps_min, sig_min = get_min(10)

    plot_solution = True
    plot_derived_quantity = True

    if plot_solution:
        print('**********************')
        print('**********************')
        print('Plotting solution ...')
        for i in xrange(from_frame,to_frame+1):
            plot_q(frame=i,plot_pcolor=True,zlimits=[eps_min,eps_max])
            print ('frame '+str(i)+' plotted')

    if plot_derived_quantity:
        print('**********************')
        print('**********************')
        print('Plotting p ...')
        for i in xrange(from_frame,to_frame+1):
            plot_p(frame=i,plot_pcolor=True,zlimits=[sig_min,sig_max])
            print ('frame '+str(i)+' plotted')
