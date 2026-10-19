#
#   plot_guess                       
#                               
#   Script to plot an initial flowfield guess created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_guess.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the guess data
    filename = 'out_guess_' + sys.argv[-1] + '.bin'

    # Read the case from file
    g = read_case(filename)

    # Open figure window and open four subplots
    fig,ax = plt.subplots(2,2,sharex=True,sharey=True,figsize=[14.4,7.2]); 
    fig.tight_layout()

    # Set subplot aspect ratios as equal and remove axes labels
    ax = ax.flatten()
    for a in ax:
        a.set_aspect('equal',adjustable='box'); a.axis('off')

    # Plot the primary flow variables to show the guess
    fieldnames = ['ro','roe','rovx','rovy']
    for n,name in enumerate(fieldnames):
 
        # Plot filled contour levels
        hc = ax[n].pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

  	# Add colorbar with variable name
        colorbar(hc,name)

        # Draw the walls of the block
        plot_wall(ax[n],g)

    # Show all the plots
    plt.savefig(f'guess_plot_{sys.argv[-1]}.png', dpi=300, bbox_inches='tight')

    
main()


