#
#   plot_advection
#                               
#   Script to plot the flowfield calculated using the example Fortran advection
#   program
#
#   Change to the directory containing the output file and execute the script 
#   with "python path_to_script/plot_advection.py"

# Import other functions, including those written for this project and numpy etc
from routines import *

# Put the bulk of your program in a function called "main"
def main():

    # Dictionaries, like the Fortran derived data types, can be used to store
    # anything you like. They also simplify the syntax when you want to pass
    # a group of multiple variables into a function. First it is initialised as
    # empty.
    g = {}

    # Read the output from the advection solver, "loadtxt" will read everything
    # into a single 2D array in one line
    arr = np.loadtxt('advection_output.txt')

    # The data is arranged into the dictionary to make the code more readable
    # for anything else we might like to do. In Python you can also put multiple 
    # statements on one line if they are terminated with a semicolon. Yay! Just 
    # like Fortran!
    g['x'] = arr[0,:]; g['phi'] = arr[1,:];

    # Matplotlib has everything we need to inspect our 1D and 2D solutions, in
    # this script we will plot the value of "phi" and how it varies in the
    # x-direction. First a figure window is opened, one set of axes is drawn and
    # its handle is stored in "ax".
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); 

    # To draw pretty graphs we have to add axes labels, choose a nice colour 
    # scheme and add the gridlines
    ax.set_xlabel('x'); ax.set_ylabel('phi'); cols = gen_cols();
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.tick_params(direction='in')

    # Plot the data as a single line with dots for markers to show the
    # resolution. You can continue across multiple lines in Python if you have
    # brackets open, if not use a "\" character.
    ax.plot(g['x'],g['phi'],'.-',color=cols[0,:],linewidth=0.5,
        label='Upwind advection solution')

    # Add a legend to keep track of multiple sets of data for your reports,
    # adding the "label" when you create the plot object is useful as it gets
    # used automatically in the legend.
    ax.legend()

    # Save the plot to a PNG file.
    plt.savefig('advection_plots.png')

    # Oh no! The dissipation is awful for this first order scheme isn't it, 
    # where has the crisp step function gone?


# This is where "main" actually gets executed
main()


