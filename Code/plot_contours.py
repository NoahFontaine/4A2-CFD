#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)

    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    g = calc_secondary(av,g)    

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    # INSERT
    cut = cut_i(g,0)
    mass_av = np.mean(cut['ro'] * cut['v'] * cut['area'])
    area_av = np.mean(cut['area'])
    g['p_ref'] = mass_av**2 / (av['gam'] * av['rgas'] * area_av**2)
    g['pstag_ref'] = g['p_ref'] + 0.5 * mass_av**2 / area_av**2
    g['cp'] = (g['p'] - g['p_ref']) / (g['pstag_ref'] - g['p_ref'])
    g['mach'] = g['v'] / (av['gam'] * av['rgas'] * g['t'])**0.5
    # END INSERT

    # Specify the parameters to plot
    fieldnames = ['cp', 'mach']; 
    colnames = ['Static pressure coefficient','Mach number']

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')
 
        # Plot filled contour levels
        hc = ax.pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Add Mach = 1 contours
        if name == 'mach':
            ax.contour(g['x'],g['y'],g['mach'],[1.0],colors='w',
                linewidths=0.5) 

        # Draw the walls of the block
        plot_wall(ax,g)

    # Show all the plots
    plt.savefig('contour_plots.png', dpi=300, bbox_inches='tight')

    
main()


