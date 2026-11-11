# Sensitivity analysis script

cfls=(1e-4 1e-3 1e-2 1e-1)

# for each cfl, edit input_bend.txt and run the simulation, then plot the convergence history to bend_cfl.csv
cfl_old=1
declare -i n=1
for cfl in "${cfls[@]}"; do
    # Edit input_bend.txt to set the CFL number
    echo "Running simulation with CFL = $cfl" 2>&1
    sed -i "s/${cfl_old} 0\.300000 0\.000100/${cfl} 0\.300000 0\.000100/" input_bend.txt

    # Run the simulation
    ../../Code/solver.x < input_bend.txt > log_$n.txt
    mv conv_bend.csv conv_bend_$n.csv

    cfl_old=$cfl
    n+=1
done