#!/usr/bin/env python
# encoding: utf-8

# Default parameters work for formation_of_diffractons/
import numpy as np
import clawpack.petclaw as pyclaw

# Medium structure
medium_type='piecewise-constant'

def qinit(state,A,x0,y0,varx,vary):
    r""" Set initial conditions:
         Gaussian stress, zero velocities."""
    x = state.grid.x.centers; y = state.grid.y.centers
    [yy,xx]=np.meshgrid(y,x)
    stress = A*np.exp(-(xx-x0)**2/(2*varx));

    #parameters from aux
    stress_rel=state.aux[2,:]
    K=state.aux[1,:]
    #initial condition
    state.q[0,:,:] = np.where(stress_rel==1,1,0) * stress/K \
                    +np.where(stress_rel==2,1,0) * np.log(stress+1)/K \
                    +np.where(stress_rel==3,1,0) * (np.sqrt(4*stress+1)-1)/(2*K)
    state.q[1,:,:]=0; state.q[2,:,:]=0

def setaux(x,y, KA, KB, rhoA, rhoB, stress_rel):
    r"""Return an array containing the values of the material
        coefficients.

        aux[0,i,j] = rho(x_i, y_j)              (material density)
        aux[1,i,j] = K(x_i, y_j)                (bulk modulus)
        aux[2,i,j] = stress-strain relation type at (x_i, y_j)
    """
    aux = np.empty((4,len(x),len(y)), order='F')
    if medium_type == 'piecewise-constant':
        yfrac = y - np.floor(y)
        xfrac = x - np.floor(x)
        # create a meshgrid out of xfrac and yfrac
        [yf,xf] = np.meshgrid(yfrac,xfrac)
        # density 
        aux[0,:,:] = rhoA*(yf<=0.5) + rhoB*(yf >0.5)
        #Young's modulus
        aux[1,:,:] = KA  *(yf<=0.5) + KB  *(yf >0.5)
        # Stress-strain relation
        aux[2,:,:] = stress_rel

    elif medium_type == 'sinusoidal' or medium_type == 'smooth_checkerboard':
        [yy,xx]=np.meshgrid(y,x)
        Amp_p=np.abs(rhoA-rhoB)/2; offset_p=(rhoA+rhoB)/2
        Amp_E=np.abs(KA-KB)/2; offset_E=(KA+KB)/2
        if medium_type == 'sinusoidal':
            frec_x=2*np.pi; frec_y=2*np.pi
            fun=np.sin(frec_x*xx)*np.sin(frec_y*yy)
        else:
            sharpness = 10
            fun_x=xx*0; fun_y=yy*0
            for i in xrange(0,1+int(np.ceil((x[-1]-x[0])/(0.5)))):
                fun_x=fun_x+(-1)**i*np.tanh(sharpness*(xx-i*0.5))
            for i in xrange(0,1+int(np.ceil((y[-1]-y[0])/(0.5)))):
                fun_y=fun_y+(-1)**i*np.tanh(sharpness*(yy-i*0.5))
            fun=fun_x*fun_y
        aux[0,:,:]=Amp_p*fun+offset_p
        aux[1,:,:]=Amp_E*fun+offset_E
        aux[2,:,:]=stress_rel

    return aux

def b4step(solver,state):
    r"""This routine does three things:
    
        1. Put in aux[3,:,:] the value of q[0,:,:] (eps). 
           This is required in rptpv.f.
           Only used by classic (not SharpClaw).

        2. Set the solution to zero in half of the domain at a specified time.

        3. Change the boundary conditions to periodic at a specified time.
    """
    state.aux[3,:,:] = state.q[0,:,:]

    # To set to 0 1st 1/2 of the domain. Used in rect domains with PBC in x
    if state.problem_data['turnZero_half_2D']==1:
        if state.t>=state.problem_data['t_turnZero'] and state.t<=state.problem_data['t_turnZero']+1:
            Y,X = np.meshgrid(state.grid.y.centers,state.grid.x.centers)
            state.q = state.q * (X>25)
            
    if state.problem_data['change_BCs']==1:
        if state.t>=state.problem_data['t_change_BCs']:
            solver.bc_lower[0]=pyclaw.BC.periodic
            solver.bc_upper[0]=pyclaw.BC.periodic
            solver.aux_bc_lower[0]=pyclaw.BC.periodic
            solver.aux_bc_upper[0]=pyclaw.BC.periodic

def compute_stress(state):
    """ Compute stress from strain and store in state.p."""
    K=state.aux[1,:,:]
    stress_rel=state.aux[2,:,:]
    eps=state.q[0,:,:]
    state.p[0,:,:] = np.where(stress_rel==1,1,0) * K*eps \
                    +np.where(stress_rel==2,1,0) * (np.exp(eps*K)-1) \
                    +np.where(stress_rel==3,1,0) * K*eps+K**2*eps**2

def total_energy(state):
    rho = state.aux[0,:,:]; K = state.aux[1,:,:]
    
    u = state.q[1,:,:]/rho
    v = state.q[2,:,:]/rho
    kinetic=rho * (u**2 + v**2)/2.

    eps = state.q[0,:,:]
    sigma = np.exp(K*eps) - 1.
    potential = (sigma-np.log(sigma+1.))/K

    dx=state.grid.delta[0]; dy=state.grid.delta[1]
    state.F[0,:,:] = (potential+kinetic)*dx*dy 

def gauge_stress(q,aux):
    p = np.exp(q[0]*aux[1])-1
    return [p]

def moving_wall_BC(state,dim,t,qbc,num_ghost):
    if dim.on_lower_boundary:
        qbc[0,:num_ghost,:]=qbc[0,num_ghost,:]
        qbc[2,:num_ghost,:]=qbc[2,num_ghost,:]
        t0=(t-10)/10
        a1=0.2;
        if abs(t0)<=1.: vwall = -a1/2.*(1.+np.cos(t0*np.pi))
        else: vwall=0.
        for ibc in xrange(num_ghost-1):
            qbc[1,num_ghost-ibc-1,:] = 2*vwall - qbc[1,num_ghost+ibc,:]


def setup(KA=5./8, KB=5./2, rhoA=8./5, rhoB=2./5, stress_rel=2, 
              oscillating_wall=False, square_domain=False,
              initial_amplitude=1, varx=5., Nx=32, Ny=128, tfinal=500,
              outdir='./_output',solver_type='sharpclaw'):
              
    """
    Solve the p-system in 2D with variable coefficients
    """
    # Domain
    x_lower=0.0; x_upper=200.00
    y_lower=0.0; y_upper=1.0

    if square_domain:
        x_upper = 100.
        y_upper = 100.
    else:
        x_upper = 200.
        y_upper = 1.
    # cells per layer
    mx=(x_upper-x_lower)*Nx
    my=(y_upper-y_lower)*Ny

    # Initial condition parameters
    x0=0.0 # Center of initial perturbation
    y0=0.25 # Center of initial perturbation
    vary=5.0 # Variance (in y) of initial Gaussian

    # Stress-strain relation:
    #   1: linear
    #   2: nonlinear (exponential)
    #   3: nonlinear (quadratic)

    # Boundary conditions
    bc_x_lower=pyclaw.BC.wall; bc_x_upper=pyclaw.BC.extrap
    bc_y_lower=pyclaw.BC.periodic; bc_y_upper=pyclaw.BC.periodic

    # Change x BCs to periodic at specified time:
    change_BCs=1
    t_change_BCs=50

    # Zero out 1st half of the domain at specified time:
    turnZero_half_2D=1 #flag
    t_turnZero=50
    
    num_output_times = tfinal

    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D()
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.cfl_max = 0.5
        solver.cfl_desired = 0.45
        solver.dimensional_split=False

    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()
        solver.cfl_max = 2.5
        solver.cfl_desired = 2.45


    if stress_rel < 3:
        from clawpack import riemann
        solver.rp = riemann.psystem_2D
    elif stress_rel==3:
        import psystem_quadratic_2D
        solver.rp = psystem_quadratic_2D

    solver.num_waves = 2

    solver.bc_lower     = [bc_x_lower, bc_y_lower]
    solver.bc_upper     = [bc_x_upper, bc_y_upper]
    solver.aux_bc_lower = [bc_x_lower, bc_y_lower]
    solver.aux_bc_upper = [bc_x_upper, bc_y_upper]

    if oscillating_wall:
        # This code assumes we'd never use an initial condition
        # and this boundary condition together.
        initial_amplitude = 0
        solver.user_bc_lower = moving_wall_BC
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.bc_lower[0] = pyclaw.BC.custom

    solver.fwave = True
    solver.before_step = b4step
        
    claw = pyclaw.Controller()
    claw.tfinal = tfinal
    claw.solver = solver
    claw.outdir = outdir
    claw.num_output_times = num_output_times

    # Domain
    x = pyclaw.Dimension('x',x_lower,x_upper,mx)
    y = pyclaw.Dimension('y',y_lower,y_upper,my)
    domain = pyclaw.Domain( [x,y] )
    num_eqn = 3
    num_aux = 4
    state = pyclaw.State(domain,num_eqn,num_aux)

    #Set global parameters
    state.problem_data = {}
    state.problem_data['turnZero_half_2D'] = turnZero_half_2D
    state.problem_data['t_turnZero']       = t_turnZero
    state.problem_data['change_BCs']       = change_BCs
    state.problem_data['t_change_BCs']     = t_change_BCs

    grid = state.grid
    state.aux = setaux(grid.x.centers,grid.y.centers, KA, KB, rhoA, rhoB, stress_rel)
    qinit(state,initial_amplitude ,x0,y0,varx,vary)

    claw.solution = pyclaw.Solution(state,domain)
    claw.num_output_times = num_output_times

    state.mp = 1
    state.mF = 1
    claw.compute_p = compute_stress
    claw.compute_F = total_energy

    claw.solution.state.grid.add_gauges([[25.0,0.75],[50.0,0.75],[75.0,0.75],[25.0,1.25],[50.0,1.25],[75.0,1.25]])
    solver.compute_gauge_values = gauge_stress
    # Do we need this?
    claw.solution.state.keep_gauges = True

    # This saves time on Shaheen, but otherwise one may wish to turn it on:
    claw.write_aux_init = False

    #Solve
    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
