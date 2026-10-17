#
#   plot_coord
#                               
#   Script to plot a mesh created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_coord.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filename to read the grid data
    filename = 'out_coord_' + sys.argv[-1] + '.bin'

    # Read the case from file
    g = read_case(filename)

    # Open figure window and set the axes to be equal
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('x / m'); ax.set_ylabel('y / m');
    ax.set_aspect('equal',adjustable='box'); ax.tick_params(direction='in')

    # Plot the mesh coordinates to show the cells
    ax.plot(g['x'],g['y'],color=cols[0,:],linewidth=0.5)
    ax.plot(np.transpose(g['x']),np.transpose(g['y']),color=cols[0,:],
        linewidth=0.5)

    # Draw the boundary of the block
    plot_bound(ax,g)

    # Show all the plots
    plt.savefig('mesh_plot_' + sys.argv[-1] + '.png')

    
main()


