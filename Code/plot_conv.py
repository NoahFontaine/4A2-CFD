#
#   plot_conv                      
#                               
#   Script to plot convergence history of a run executed using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_conv.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the guess data
    filename = 'conv_' + sys.argv[-1] + '.csv'

    # Read the history from file
    l = read_conv(filename)

    # Open figure window for all residual data
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('Iteration'); ax.set_ylabel('Residual');

    # Set y-axis as a log scale and turn on the gridlines
    ax.set_yscale('log'); ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.grid(linestyle='-',color=[0.8,0.8,0.8],linewidth=0.5,which='minor',
        axis='y')

    # 2x4 nested list of residuals and names for legend
    fieldnames = [['dro_max', 'droe_max', 'drovx_max', 'drovy_max'],
        ['dro_avg', 'droe_avg', 'drovx_avg', 'drovy_avg']]
    colnames = ['dro', 'droe', 'drovx', 'drovy']

    # Plot the residuals changing with time
    for m in range(len(fieldnames)):
        for n,name in enumerate(fieldnames[m]):
            ax.plot(l['nstep'],l[name],color=cols[n,:])

    # Add the legend for colours which are shared for max and average residuals
    ax.legend(colnames)

    # Show all the plots
    plt.savefig(f'convergence_plot_{sys.argv[-1]}.png', dpi=300, bbox_inches='tight')

    
main()


