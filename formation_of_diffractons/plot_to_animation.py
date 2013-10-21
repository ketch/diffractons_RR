from clawpack.petclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import Axes3D as pl3
from matplotlib import rc
import matplotlib.cm as cm
import numpy as np
import os

def plot_p(frame,zlimits=None):
    if not os.path.exists('./_plots_animation/_p'): os.mkdir('./_plots_animation/_p')
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers; 
    mx=len(x); my=len(y)    
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    else:
        str_frame = "0"+str(frame)

    p=sol.state.q[0,:,:]
    # plot pcolor
    pl.figure(figsize=(8,3))
    pl.pcolormesh(xx,yy,p,cmap=cm.OrRd)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    cb = pl.colorbar();
    if zlimits is not None:
        pl.clim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05)
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    pl.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
    #pl.axis([0.25,60.25,0.25,60.25])
    if not os.path.exists('./_plots_animation/_p/pcolor'): os.mkdir('./_plots_animation/_p/pcolor')
    pl.savefig('./_plots_animation/_p/pcolor/p_frame'+str_frame+'.png')
    pl.close()
    # plot slices
    pl.figure(figsize=(8,3))
    pl.gcf().subplots_adjust(left=0.15)
    pl.plot(x,p[:,3*my/4.],'-r',linewidth=1)
    pl.plot(x,p[:,my/4.],'-b',linewidth=1)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xlabel('x',fontsize=20)
    pl.ylabel('Stress',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    if zlimits is not None:
        pl.ylim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05)
    if not os.path.exists('./_plots_animation/_p/slices'): os.mkdir('./_plots_animation/_p/slices')
    pl.savefig('./_plots_animation/_p/slices/p_frame'+str_frame+'_slices.eps')
    pl.close()

    #surface plots
    fig3 = pl.figure(figsize=(8,3))
    ax = fig3.gca(projection='3d')
    ax.plot_surface(xx,yy,p,linewidth=0)
    ax.set_axis_off()
    if zlimits is not None:
        ax.set_zlim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05)
    if not os.path.exists('./_plots_animation/_p/surface'): os.mkdir('./_plots_animation/_p/surface')
    pl.savefig('./_plots_animation/_p/surface/p_frame'+str_frame+'_surface.png')

def plot_q(frame,meqn,zlimits=None):
    if not os.path.exists('./_plots_animation/_q'+str(meqn)): os.mkdir('./_plots_animation/_q'+str(meqn))
    sol=Solution(frame,file_format='petsc',read_aux=False)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers; 
    mx=len(x); my=len(y)    
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    else:
        str_frame = "0"+str(frame)

    q=sol.state.q[meqn,:,:]
    # plot pcolor
    pl.figure(figsize=(8,3))
    pl.pcolormesh(xx,yy,q,cmap=cm.OrRd)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    cb = pl.colorbar();
    if zlimits is not None:
        pl.clim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05)
    imaxes = pl.gca(); pl.axes(cb.ax)
    pl.yticks(fontsize=20); pl.axes(imaxes)
    pl.axis([np.min(x),np.max(x),np.min(y),np.max(y)])
    #pl.axis([0.25,60.25,0.25,60.25])
    if not os.path.exists('./_plots_animation/_q'+str(meqn)+'/pcolor'): os.mkdir('./_plots_animation/_q'+str(meqn)+'/pcolor')
    pl.savefig('./_plots_animation/_q'+str(meqn)+'/pcolor/q'+str(meqn)+'_frame'+str_frame+'_pcolor.png')
    pl.close()
    # plot slices
    pl.figure(figsize=(8,3))
    pl.gcf().subplots_adjust(left=0.15)
    if meqn==2:
        pl.plot(x,q[:,0],'-r',linewidth=1)
        pl.plot(x,q[:,-1],'-b',linewidth=1)
    else:
        pl.plot(x,q[:,3*my/4.],'-r',linewidth=1)
        pl.plot(x,q[:,my/4.],'-b',linewidth=1)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xlabel('x',fontsize=20)
    pl.ylabel('q'+str(meqn),fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    if zlimits is not None:
        pl.ylim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05)
    if not os.path.exists('./_plots_animation/_q'+str(meqn)+'/slices'): os.mkdir('./_plots_animation/_q'+str(meqn)+'/slices')
    pl.savefig('./_plots_animation/_q'+str(meqn)+'/slices/q'+str(meqn)+'_frame'+str_frame+'_slices.eps')
    pl.close()

    #surface plots
    fig3 = pl.figure(figsize=(8,3))
    ax = fig3.gca(projection='3d')
    ax.plot_surface(xx,yy,q,linewidth=0)
    ax.set_axis_off()
    if zlimits is not None:
        ax.set_zlim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05)
    if not os.path.exists('./_plots_animation/_q'+str(meqn)+'/surface'): os.mkdir('./_plots_animation/_q'+str(meqn)+'/surface')
    pl.savefig('./_plots_animation/_q'+str(meqn)+'/surface/q'+str(meqn)+'_frame'+str_frame+'_surface.png')

def get_extremum_p(frame):
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    sig=sol_p.state.q[0,:]
    return [np.min(sig), np.max(sig)]

def get_extremum_q(frame,meqn):
    sol_q=Solution(frame,file_format='petsc',read_aux=False)
    q=sol_q.state.q[meqn,:]
    return [np.min(q), np.max(q)]

if __name__== "__main__":
    if not os.path.exists('./_plots_animation'): os.mkdir('./_plots_animation')
    
    from_frame = 0
    to_frame = 400

    zlim_p = get_extremum_p(400)
    zlim_q0 = get_extremum_q(400,0)
    zlim_q2 = get_extremum_q(400,2)

    print('**********************')
    print('Plotting ...')
    for i in xrange(from_frame,to_frame+1):
        plot_p(i,zlimits=zlim_p)
        plot_q(i,0,zlimits=zlim_q0)
        plot_q(i,2,zlimits=zlim_q2)
        print ('frame '+str(i)+' plotted')

