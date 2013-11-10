""" Reproduce all figures from the paper that show results of Clawpack runs.
    Note that the default resolution (values of Nx, Ny) set here is much lower
    than that used in the paper.  This is to reduce the computational cost
    to something feasible on a good workstation in a few hours.  But for most
    cases, this is insufficient to accurately capture the dynamics.  To
    truly reproduce the figures in the paper, you should set Nx=32, Ny=128.
    But for this you will need (as of 2013) at least a large cluster.
    
    Figure names correspond to the arXiv version of the paper.

    For example, to reproduce figure 3a:

       > import reproduce
       > reproduce.reproduce_figure('3a',Nx=32, Ny=128)

    This puts the simulation results in a subdirectory of computed_data/,
    and the figures in the figures/ directory.

    You can also just plot the figure from the provided simulation results:

       > import reproduce
       > reproduce.figure_3a(use_saved_data=True)

    This uses data stored in a subdirectory of saved_data/.
"""
from psystem import setup
import os.path
import waves_2D_plots
import make_plots

fig_param = {

'3a' : {
    'description' : """Linear medium with constant sound speed.""",
    'frame' : 375,
    'xshift': 200,
    'simulation_params' : {
        'KA'    : 4., 'rhoA'  : 4., 'KB'    : 1., 'rhoB'  : 1.,
        'stress_rel' : 1,
        'oscillating_wall' : True,
    } },

'3b' : {
    'description' : """Linear medium with variable sound speed.""",
    'frame' : 375,
    'xshift': 400,
    'simulation_params' : {
        'KA'    : 5./8, 'rhoA'  : 8./5, 'KB'    : 5./2, 'rhoB'  : 2./5,
        'stress_rel' : 1,
        'oscillating_wall' : True,
    } },

'3c' : {
    'description' : """Nonlinear medium with constant sound speed.""",
    'frame' : 375,
    'xshift': 200,
    'simulation_params' : {
        'KA'    : 4., 'rhoA'  : 4., 'KB'    : 1., 'rhoB'  : 1.,
        'stress_rel' : 2,
        'oscillating_wall' : True,
    } },

'3d' : {
    'description' : """Nonlinear medium with variable sound speed.""",
    'frame' : 375,
    'xshift': 400,
    'simulation_params' : {
        'KA'    : 5./8, 'rhoA'  : 8./5, 'KB'    : 5./2, 'rhoB'  : 2./5,
        'stress_rel' : 2,
        'oscillating_wall' : True,
    } },

'4' : {
    'description' : """Diffractons for comparison with homogenized solution.""",
    'frame' : 120,
    'xshift' : 0,
    'simulation_params' : {
        'KA'    : 5./8, 'rhoA'  : 8./5, 'KB'    : 5./2, 'rhoB'  : 2./5,
        'stress_rel' : 2
    } },

'5' : {
    'description' : """Generation of a long train of diffractons.""",
    'frame' : 1470,
    'xshift': 1800,
    'simulation_params' : { 'varx' : 50. }
    },
'9' : {},
'10' : {},
'11' : {},

'13' : {
    'description' : """Formation of diffractons in a medium with quadratic nonlinearity.""",
    'frame' : 120,
    'xshift': 0,
    'simulation_params' : { 'stress_rel' : 3 }
    },

'14b' : {
    'description' : """Square domain with mismatched impedance.""",
    'frame' : 65,
    'xshift': 0,
    'simulation_params' : { 
        'KA'    : 4., 'rhoA'  : 4., 'KB'    : 1., 'rhoB'  : 1.,
        'square_domain' : True,
        'initial_amplitude' : 5
    } },

'14c' : {
    'description' : """Square domain with mismatched sound speed.""",
    'frame' : 65,
    'xshift': 0,
    'simulation_params' : { 
        'KA'    : 5./8, 'rhoA'  : 8./5, 'KB'    : 5./2, 'rhoB'  : 2./5,
        'square_domain' : True,
        'initial_amplitude' : 5
    } },

'14d' : {
    'description' : """Square domain with mismatched impedance and sound speed.""",
    'frame' : 65,
    'xshift': 0,
    'simulation_params' : { 
        'KA'    : 16., 'rhoA'  : 1., 'KB'    : 1., 'rhoB'  : 1.,
        'square_domain' : True,
        'initial_amplitude' : 5
    } }
}

def reproduce_figure(figname, Nx=2, Ny=8, plot=True, use_saved_data=True):
    """ If use_saved_data==False, re-run simulation and plot results.
        If use_saved_data==True, plot results from data provided."""
    params = fig_param[figname]
    frame = params['frame']
    filename = 'claw.pkl'+str(frame).zfill(4)
    if use_saved_data:
        dirname = 'saved_data/figure_'+figname
    else:
        dirname = 'computed_data/figure_'+figname
        if not os.path.exists('./'+dirname+'/'+filename):
            simulation_params = params['simulation_params']
            simulation_params['tfinal'] = frame
            simulation_params['Nx'] = Nx
            simulation_params['Ny'] = Ny
            simulation_params['outdir'] = './'+dirname

            claw = setup(**params['simulation_params'])
            claw.run()

    if plot:
        if (fig=='14a'):
            waves_2D_plots.plot(0,'figure_14a',axis_limits=[0.25, 10.25, 0.25, 10.25])
        elif (fig=='14b'):
            waves_2D_plots.plot(frame,'figure_14b',clim=[-0.3,0.3])
        elif (fig=='14c'):
            waves_2D_plots.plot(65,'figure_14c',clim=[-0.3,0.3])
        elif (fig=='14d'):
            waves_2D_plots.plot(65,'figure_14d',clim=[-0.3,0.3])
        else: 
            zlimits=make_plots.get_extremum(frame,dirname)
            for fields in ('stress','strain'):
                make_plots.plot(frame,dirname,'fig'+figname,fields,zlimits=zlimits,xshift=params['xshift'])
      
if __name__=="__main__":
    # Plotting turned off here so you can run in parallel from the command line, i.e.
    # mpiexec -np 8 python reproduce.py

    for fig in ('3a','3b','3c','3d','5','14b','14c','14d'):
        reproduce_figure(fig, plot = True, use_saved_data = True)
