#
#   generate_case
#                               
#   Script to generate geometry and input files for a range of different cases
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/generate_case.py casename"

# Import modules and functions
from routines import *

def gen_bend(casename):
    # Create a 90 degree bend with extended straight sections
 
    # Settings and boundary conditions
    av = default_settings(casename)
    av['p'] = 50e3; av['alpha'] = 45;

    # Initialise the dictionary to store the curves
    geom = {}

    # Geometrical parameters and resolution
    r_a = 0.8; r_b = 1; th_1 = np.deg2rad(135); th_2 = np.deg2rad(45); 
    d = np.atleast_1d(0.25); nn = 361;

    # Create the bend
    th = np.linspace(th_1,th_2,nn)
    x_a = r_a * np.cos(th); y_a = r_a * np.sin(th);
    x_b = r_b * np.cos(th); y_b = r_b * np.sin(th);

    # Add linear extensions
    dx_1 = d * np.cos(th_1 + np.pi/2); dy_1 = d * np.sin(th_1 + np.pi/2);
    dx_2 = d * np.cos(th_2 - np.pi/2); dy_2 = d * np.sin(th_2 - np.pi/2);
    geom['x_a'] = np.concatenate([x_a[0] + dx_1, x_a, x_a[-1] + dx_2])
    geom['y_a'] = np.concatenate([y_a[0] + dy_1, y_a, y_a[-1] + dy_2])
    geom['x_b'] = np.concatenate([x_b[0] + dx_1, x_b, x_b[-1] + dx_2])
    geom['y_b'] = np.concatenate([y_b[0] + dy_1, y_b, y_b[-1] + dy_2])

    return(av,geom)

################################################################################

def gen_bump(casename):
    # Create a duct with smooth bump protruding from lower wall

    # Settings and boundary conditions
    av = default_settings(casename)
    av['p'] = 85e3;

    # Initialise the dictionary to store the curves
    geom = {}

    # Geometrical parameters and resolution
    x = [0,0.2,0.7,1.5]; y = [0,0.06,0.4]; nn = 401;

    # Create straight line for upper wall
    geom['x_b'] = np.linspace(x[0],x[3],nn); geom['y_b'] = np.ones(nn) * y[2];

    # Create sinusoidal bump on lower wall 
    geom['x_a'] = geom['x_b']; geom['y_a'] = np.ones(nn) * y[0];
    q = (geom['x_a'] > x[1]) & (geom['x_a'] < x[2])
    dx = x[2] - x[1]; dy = y[1] - y[0]; 
    th = np.pi * (2 * (geom['x_a'][q] - x[1]) / dx - 1) 
    geom['y_a'][q] = 0.5 * dy * (1 + np.cos(th)) + y[0]

    # Copy lower bump to upper wall for symmetric shock case
    #geom['x_b'] = geom['x_a']; geom['y_b'] = -geom['y_a'] + y[2];

    return(av,geom)

################################################################################

def gen_tunnel(casename):
    # Create a supersonic wind tunnel with two throats, note that this function
    # does not currently specify time-varying boundary conditions, nor does the
    # write_settings function support this capability.

    # Settings and boundary conditions
    av = default_settings(casename)
    av['pstag'] = 200e3; av['p'] = 190e3;
    av['ni'] = 401; av['nj'] = 81; 
    av['sfac'] = 0.15; av['nsteps'] = 12000;

    # Initialise the dictionary to store the curves
    geom = {}

    # Final curve resolution
    nn = 2001;

    # First throat geometry from Vorst, Leon A., "Design and calibration of a 
    # Mach 2.5 two dimensional wind tunnel nozzle" (1970) 
    x = 0.1 * np.array([-4.86,-4.50,-4.00,-3.50,-3.00,-2.50,-2.00,-1.50,-1.00,
        -0.50,0.00,0.50,1.00,1.500,2.000,2.500,3.000,3.500,4.000,4.500,5.000,
        5.500,6.000,6.500,7.000,7.500,8.000,8.500,9.000,9.500,10.000,10.500,
        11.000,11.500,12.000,12.500,13.000,13.500,14.000,14.500,15.000,15.500,
        16.081])
    y = 0.1 * np.array([6.000,4.985,3.840,2.948,2.279,1.798,1.467,1.280,1.180,
        1.143,1.138,1.141,1.176,1.245,1.341,1.449,1.563,1.681,1.801,1.920,2.030,
        2.133,2.227,2.314,2.393,2.466,2.531,2.578,2.655,2.703,2.750,2.793,2.830,
        2.863,2.892,2.918,2.939,2.958,2.972,2.983,2.992,2.996,3.000])

    # Add on inlet contraction
    dx_inlet = np.array([-0.8,-0.12]); dy_inlet = 0.15;
    x = np.concatenate([dx_inlet + x[0],x])
    y = np.concatenate([np.array([1,1]) * dy_inlet + y[0],y])
    x_curv = [x[0],x[-1]]

    # Add on test section    
    dx_test = np.atleast_1d(1.2);
    x = np.concatenate([x,dx_test + x[-1]])
    y = np.concatenate([y,np.atleast_1d(y[-1])])

    # Add on exit contraction and diffuser
    dx_exit = np.array([0.6,2.4,3]); dy_exit = np.array([-0.07,0.25,0.25]);
    x = np.concatenate([x,dx_exit + x[-1]])
    y = np.concatenate([y,dy_exit + y[-1]])

    # Construct two interpolants, pchip for curved and linear for straight
    q = (x >= x_curv[0]) & (x <= x_curv[1])
    f_curv = interp.PchipInterpolator(x[q],y[q])
    f_lin = interp.interp1d(x,y)

    # Reinterpolate to consistent and increased resolution
    x_interp = np.linspace(x[0],x[-1],nn)
    y_interp = f_lin(x_interp) 
    q = (x_interp >= x_curv[0]) & (x_interp <= x_curv[1])
    y_interp[q] = f_curv(x_interp[q]) 

    # Assemble the curves
    geom['x_a'] = x[[0,-1]]; geom['y_a'] = np.zeros(2);
    geom['x_a'] = x_interp; geom['y_a'] = -y_interp;
    geom['x_b'] = x_interp; geom['y_b'] = y_interp;

    return(av,geom)

################################################################################

def gen_corner(casename):
    # Create a 90 degree turning duct to investigate unsteady shock structures

    # Settings and boundary conditions
    av = default_settings(casename)
    av['pstag'] = 314.64e3; av['tstag'] = 402.88; av['p'] = 100e3;
    av['ni'] = 301; av['nj'] = 101; 
    av['sfac'] = 0.15; av['nsteps'] = 12000;

    # Initialise the dictionary to store the curves
    geom = {}

    # Geometrical parameters to define the geometry
    L = 0.2;

    # Inner boundary of the corner
    geom['x_a'] = np.array([-L, L, L])
    geom['y_a'] = np.array([0, 0, -2 * L])

    # Outer boundary of the corner
    geom['x_b'] = np.array([-L, 2 * L, 2 * L]); 
    geom['y_b'] = np.array([L, L, -2 * L]);

    return(av,geom)

################################################################################

def gen_waves(casename):
    # Create a duct sharp turning to perform supersonic compression with perfect
    # cancellation

    # Settings and boundary conditions
    av = default_settings(casename)

    # Set static pressure as ratio of inlet stagnation - note this is static
    # pressure at the inlet plane for a supersonic inlet boundary condition
    av['p'] = av['pstag'] * 0.1278;

    # Initialise the dictionary to store the curves
    geom = {}

    # Geometrical parameters from 3A3 - Two-dimensional compressible flow - 
    # Examples Paper 2 - Q3 
    theta = np.deg2rad(10); beta_1 = np.deg2rad(39.314); 
    beta_2 = np.deg2rad(49.429); dy_1 = 0.1; dx_bc = 0.2; 

    # Dimensions of corners and intersection
    dx_1 = dy_1 / np.tan(beta_1)
    dx_2 = (dy_1 - dx_1 * np.tan(theta)) / (np.tan(theta) + np.tan(beta_2 - theta))
    dy_2 = dx_2 * np.tan(beta_2 - theta)

    # Create straight lower wall with coordinate point at shock intersection
    geom['x_a'] = np.array([-dx_bc, dx_1 , dx_1 + dx_2 + dx_bc])
    geom['y_a'] = np.array([0, 0, 0])

    # Create upper wall with two corners
    geom['x_b'] = np.array([-dx_bc, 0, dx_1 + dx_2, dx_1 + dx_2 + dx_bc]); 
    geom['y_b'] = np.array([dy_1, dy_1, dy_2, dy_2]);

    return(av,geom)

################################################################################

def gen_tube(casename):
    # Create a simple straight section of duct, this currently assumes the inlet
    # plane is at i = 1 and the outlet is at i = ni as per the basic test cases
 
    # Settings and boundary conditions
    av = default_settings(casename)
    av['ni'] = 101; av['nj'] = 5;

    # Initialise the dictionary to store the curves
    geom = {}

    # Create the geometry directly in numpy arrays
    geom['x_a'] = np.array([0.0,1.0]); geom['x_b'] = geom['x_a'];
    geom['y_a'] = np.zeros(2); geom['y_b'] = np.ones(2) * 0.05;

    return(av,geom)

################################################################################

def gen_naca(casename):
    # Create a single cambered NACA aerofoil in a straight wind tunnel

    # Optional elliptic grid smoothing
    nsmooth = 10
 
    # Settings and boundary conditions
    av = default_settings(casename)
    av['ni'] = -1; av['nj'] = -1;
    av['p'] = 90000;

    # Initialise the dictionary to store the curves
    geom = {};

    # Specify NACA 2412 coefficients and surface resolution
    m = 0.02; p = 0.4; t = 0.12; ns = 101; 

    # Construct the camberline
    x = polyspace(0,1,0.005/float(ns),1.0/float(ns),ns); yc = np.zeros(ns);  
    q = x < p; yc[q] = m * (2 * p * x[q] - x[q]**2) / p**2;
    q = x >= p; yc[q] = m * ((1 - 2 * p) + 2 * p * x[q] - x[q]**2) / (1 - p)**2;

    # Gradient of the camberline to aid construction later
    dycdx = np.zeros(ns)
    q = x < p; dycdx[q] = 2 * m * (p - x[q]) / p**2; 
    q = x >= p; dycdx[q] = 2 * m * (p - x[q]) / (1 - p)**2;
    th = np.arctan(dycdx)
   
    # Construct the thickness - modified for zero trailing edge thickness
    yt = 5 * t * (0.2969 * x**0.5 - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 
        0.1036 * x**4)

    # Create the upper and lower surfaces by adding thickness to the camberline
    geom['x_a'] = x - yt * np.sin(th); geom['x_b'] = x + yt * np.sin(th);
    geom['y_a'] = yc + yt * np.cos(th); geom['y_b'] = yc - yt * np.cos(th);

    # Domain boundary curves
    geom['x_c'] = np.array([-4.0,5.0]); geom['x_d'] = geom['x_c'];
    geom['y_c'] = -2.5 * np.ones(2); geom['y_d'] = np.ones(2) * 2.5;
  
    # Create a 5-block mesh around the aerofoil
    g = mesh_aerofoil(geom,25,93,25,21,0.2,[-0.8,0])

    # Set walls on the aerofoil and the tunnel
    for n in range(g['nn']):
        g['wall'].append(np.zeros([g['ni'][n],g['nj'][n]],dtype=int))
    g['wall'][0][:,0] = 1;
    g['wall'][1][:,0] = 1; g['wall'][3][:,0] = 1; g['wall'][4][:,0] = 1;
    g['wall'][2][:,-1] = 1; g['wall'][3][:,-1] = 1; g['wall'][4][:,-1] = 1;

    # Smooth the grid
    g = elliptic_smooth(g,nsmooth)

    return(av,geom,g)

################################################################################

def gen_turbine_c(casename):
    # Create a turbine cascade with multi-block mesh including blade C-mesh

    # Optional elliptic grid smoothing
    nsmooth = 10
 
    # Settings and boundary conditions
    av = default_settings(casename)
    av['ni'] = -1; av['nj'] = -1;
    av['p'] = 85000;

    # Initialise the dictionary to store the curves
    geom = {};

    # Read G9 turbine cascade from file, used in 4A3 wind tunnel experiment
    arr = np.loadtxt('g9.raw',skiprows=1); nn = np.shape(arr)[0];

    # Extract upper and lower coordinates from the array
    geom['x_a'] = arr[:,0]; geom['x_b'] = arr[:,2];
    geom['y_a'] = arr[:,1]; geom['y_b'] = arr[:,3];

    # Reinterpolate to desired resolution in streamwise direction
    ns = 101;
    si_a = dist(geom['x_a'],geom['y_a'],True)
    si_b = dist(geom['x_b'],geom['y_b'],True)
    si = polyspace(0,1,0.1/float(ns),1.0/float(ns),ns);
    geom['x_a'] = interp.pchip_interpolate(si_a,geom['x_a'],si)
    geom['y_a'] = interp.pchip_interpolate(si_a,geom['y_a'],si)
    geom['x_b'] = interp.pchip_interpolate(si_b,geom['x_b'],si)
    geom['y_b'] = interp.pchip_interpolate(si_b,geom['y_b'],si)

    # Domain boundary curves
    geom['x_c'] = np.array([-2.0,0.0,1.0,2.0,3.0,5.0])
    geom['y_c'] = np.array([1.5,1.5,1.2,0.5,-1.0,-3.0])
    geom['x_d'] = geom['x_c']; geom['y_d'] = geom['y_c'] + 2.99;
  
    # Create a 5-block mesh around the aerofoil
    g = mesh_aerofoil(geom,25,93,25,21,0.1,[-0.1,0.3])

    # Set walls on the turbine blade surface
    for n in range(g['nn']):
        g['wall'].append(np.zeros([g['ni'][n],g['nj'][n]],dtype=int))
    g['wall'][0][:,0] = 1;

    # Add periodic matching patches on inlet, upper and lower domain blocks and 
    # the outlet
    g['match'].append(match(4,4,np.arange(1,g['ni'][3]+1),
        np.array(1),np.arange(1,g['ni'][3]+1),np.array(g['nj'][3])))
    g['match'].append(match(5,5,np.arange(1,g['ni'][4]+1),
        np.array(1),np.arange(1,g['ni'][4]+1),np.array(g['nj'][4])))
    g['match'].append(match(2,3,np.arange(1,g['ni'][1]+1),
        np.array(1),np.arange(1,g['ni'][2]+1),np.array(g['nj'][2])))
    g['nm'] = g['nm'] + 3

    # Smooth the grid
    g = elliptic_smooth(g,nsmooth)

    return(av,geom,g)

################################################################################

def gen_turbine_h(casename):
    # Create a turbine cascade with single block sheared-H mesh

    # Optional elliptic grid smoothing
    nsmooth = 0
 
    # Settings and boundary conditions
    av = default_settings(casename)
    av['ni'] = -1; av['nj'] = -1;
    av['p'] = 85000;

    # Initialise the dictionaries to store the curves and the mesh
    geom = {}; g = {}; g['nn'] = 1;
    for var in ['x','y','ni','nj','match','n_in','n_out','wall']:
        g[var] = []; 

    # Read G9 turbine cascade from file, used in 4A3 wind tunnel experiment
    arr = np.loadtxt('g9.raw',skiprows=1); nn = np.shape(arr)[0];

    # Extract upper and lower coordinates from the array and offset by a pitch
    dy = 2.99
    geom['x_a'] = arr[:,0]; geom['x_b'] = arr[:,2];
    geom['y_a'] = arr[:,1]; geom['y_b'] = arr[:,3] + dy;

    # Salient coordinates
    x_le = geom['x_a'][0]; x_te = geom['x_a'][-1]; x_in = -2.0; x_out = 5.0; 
    dx_le = 0.001 * (x_te - x_le); dx_in = 0.025 * (x_le - x_in);
    dx_te = 0.008 * (x_te - x_le); dx_out = 0.025 * (x_le - x_in);
    y_a_out = np.atleast_1d(-2.5); y_b_out = y_a_out + dy;

    # Extend the curves up and downstream in axial direction
    geom['x_a'] = np.concatenate([np.ones(1)*x_in,geom['x_a'],np.ones(1)*x_out])
    geom['x_b'] = np.concatenate([np.ones(1)*x_in,geom['x_b'],np.ones(1)*x_out])
    geom['y_a'] = np.concatenate([geom['y_a'][[0]],geom['y_a'],y_a_out])
    geom['y_b'] = np.concatenate([geom['y_b'][[0]],geom['y_b'],y_b_out])

    # Resolution and spacing in i-direction for inlet, blade and outlet sections
    ni = [41,101,41]
    x = np.unique(np.concatenate([polyspace(x_in,x_le,dx_in,dx_le,ni[0]), 
        polyspace(x_le,x_te,dx_le,dx_te,ni[1]),
        polyspace(x_te,x_out,dx_te,dx_out,ni[2])])).reshape([-1,1]);

    # Resolution and spacing in j-direction
    nj = 73
    sj = polyspace(0,1,0.2/float(nj),0.2/float(nj),nj).reshape([1,nj])

    # Reinterpolate to desired resolution in x-coordinate
    geom['y_a'] = interp.pchip_interpolate(geom['x_a'],geom['y_a'],x)
    geom['y_b'] = interp.pchip_interpolate(geom['x_b'],geom['y_b'],x)
    geom['x_a'] = x; geom['x_b'] = x;

    # Force periodic boundaries to be periodic
    geom['y_b'][0:ni[0]] = geom['y_a'][0:ni[0]] + dy
    geom['y_b'][-ni[2]:,:] = geom['y_a'][-ni[2]:,:] + dy

    # Use j-direction spacing to create the intermediate mesh
    g['x'].append(np.tile(x.reshape([-1,1]),[1,nj]))  
    g['y'].append((1-sj).reshape([1,nj]) * geom['y_a'] + \
        sj.reshape([1,nj]) * geom['y_b'])  

    # Store the shape of the block
    g['ni'].append(np.shape(g['x'][0])[0]); g['nj'].append(nj);

    # Set walls on the turbine blade surface
    g['wall'].append(np.zeros([g['ni'][0],g['nj'][0]],dtype=int))
    g['wall'][0][ni[0]-1:-ni[2]+1,[0,-1]] = 1;

    # Add periodic matching patches upstream of the LE and downstream of TE
    g['match'].append(match(1,1,np.arange(1,ni[0]+1),
        np.array(1),np.arange(1,ni[0]+1),np.array(nj)))
    g['match'].append(match(1,1,np.arange(ni[0]+ni[1]-1,g['ni'][0]+1),
        np.array(1),np.arange(ni[0]+ni[1]-1,g['ni'][0]+1),np.array(nj)))
    g['nm'] = 2

    # Add inlet and outlet patches on the only block
    g['n_in'] = 1; g['n_out'] = 1;

    # Smooth the grid
    g = elliptic_smooth(g,nsmooth)

    return(av,geom,g)

################################################################################

def gen_multi(casename):
    # Create the entire mesh for a four element wing
 
    # Settings and boundary conditions
    av = default_settings(casename)

    # Aerofoil geometry from Suddhoo, A. and Hall, I. M. "Test cases for the 
    # plane potential flow past multi-element aerofoils" (1985). Note this paper
    # also has the exact inviscid static pressure coefficients that can be read
    # from the input file
    arr = np.loadtxt('sudhoo.raw',skiprows=1); nn = np.shape(arr)[0];
    x = arr[:,:4]; y = arr[:,4:8]; 

    # Initialise the dictionary to store the curves
    geom = {}

    # Split into upper and lower surfaces for plotting in Python
    n_a = np.arange(round(nn/2),-1,-1); n_b = np.arange(round(nn/2),nn);
    geom['x_a'] = x[n_a,:]; geom['y_a'] = y[n_a,:]; 
    geom['x_b'] = x[n_b,:]; geom['y_b'] = y[n_b,:];

    # Domain boundary curves
    geom['x_c'] = np.array([-6.0,8.0]); geom['x_d'] = geom['x_c'];
    geom['y_c'] = -4.5 * np.ones(2); geom['y_d'] = np.ones(2) * 3.5;

    # Create the H-mesh between each aerofoil element

    return(av,geom)

################################################################################

def main():
 
    # Get the casename from command line input
    casename = sys.argv[-1]

    # Create the curves for the desired case and set the boundary conditions
    if casename == 'bend':
        av,geom = gen_bend(casename)
    elif casename == 'bump':
        av,geom = gen_bump(casename)
    elif casename == 'tunnel':
        av,geom = gen_tunnel(casename)
    elif casename == 'corner':
        av,geom = gen_corner(casename)
    elif casename == 'waves':
        av,geom = gen_waves(casename)
    elif casename == 'tube':
        av,geom = gen_tube(casename)
    elif casename == 'naca':
        av,geom,g = gen_naca(casename)
    elif casename == 'turbine_c':
        av,geom,g = gen_turbine_c(casename)
    elif casename == 'turbine_h':
        av,geom,g = gen_turbine_h(casename)
    elif casename == 'multi':
        av,geom = gen_multi(casename)

    # Save the settings and curves to their input files
    write_settings(av)
    write_geom(av,geom)

    # Open figure window to plot the curves
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('x / m'); ax.set_ylabel('y / m');
    ax.set_aspect('equal',adjustable='box'); ax.tick_params(direction='in')
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)

    # Plot the geometry curves
    ax.plot(geom['x_a'],geom['y_a'],'.-',color=cols[0,:])
    ax.plot(geom['x_b'],geom['y_b'],'.-',color=cols[1,:])

    # Plot the domain curves if present
    if 'x_c' in geom:
        ax.plot(geom['x_c'],geom['y_c'],'.-',color=cols[2,:])
        ax.plot(geom['x_d'],geom['y_d'],'.-',color=cols[3,:])
    
    # Plot the mesh coordinates and patches if present
    if 'g' in locals():
        
        # Open a new figure window
        plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
        ax.set_xlabel('x / m'); ax.set_ylabel('y / m');
        ax.set_aspect('equal',adjustable='box'); ax.tick_params(direction='in')

        # Plot the mesh coordinates and walls
        for n in range(g['nn']):

            # Mesh in both i and j-directions
            ax.plot(g['x'][n],g['y'][n],color=cols[n,:],linewidth=0.5)
            ax.plot(np.transpose(g['x'][n]),np.transpose(g['y'][n]),
                color=cols[n,:],linewidth=0.5)

            # Extract the block and plot the walls
            plot_wall(ax,cut_block(g,n),False);

        # Plot the patch matching data
        for m in range(g['nm']):

            # Pull out the index lists and convert to Python indexing
            p = g['match'][m]; n_1 = p['n_1']-1; n_2 = p['n_2']-1;
            i_1 = p['i_1']-1; i_2 = p['i_2']-1; 
            j_1 = p['j_1']-1; j_2 = p['j_2']-1;

            # Extract coordinates
            x_1 = g['x'][n_1][i_1,j_1]; y_1 = g['y'][n_1][i_1,j_1];
            x_2 = g['x'][n_2][i_2,j_2]; y_2 = g['y'][n_2][i_2,j_2];

            # Plot both sides of the patch with different symbols
            ax.plot(x_1,y_1,'+',color=cols[m,:])
            ax.plot(x_2,y_2,'x',color=cols[m,:])

            # Calculate and print the error
            d_max = np.max(np.hypot(x_2 - x_1,y_2 - y_1))
            print('Patch', m+1, 'd_max =', d_max)
            
        # Plot the inlet and outlets, these are assumed to be at i = 1 and i = ni
        x_in = np.squeeze(g['x'][g['n_in']-1][0,:])
        y_in = np.squeeze(g['y'][g['n_in']-1][0,:])
        ax.plot(x_in,y_in,color=[0.8,0.8,0.8])
        x_out = np.squeeze(g['x'][g['n_out']-1][-1,:])
        y_out = np.squeeze(g['y'][g['n_out']-1][-1,:])
        ax.plot(x_out,y_out,color=[0.8,0.8,0.8])

        # Write the entire grid definition to file
        write_mesh(av,g)

    # Show all the plots
    plt.show()

    
main()


