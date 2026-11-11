# Sensitivity analysis script

nis=(15 40 75 150 300 600)
njs=(10 25 50 100 200 400)

# for each ni, edit input_bend.txt and run the simulation, then plot the convergence history to bend_cfl.csv
ni_old=600
nj_old=400
declare -i n=1
for ni in "${nis[@]}"; do
    # Edit input_bend.txt to set the CFL number
    echo "Running simulation with ni = $ni" 2>&1
    sed -i "s/${ni_old} ${nj_old}/${ni} ${njs[$((n-1))]}/" input_bend.txt
    
    # Run the simulation
    ../../Code/solver.x < input_bend.txt > log_ni.txt
    mv conv_bend.csv conv_bend_$n.csv

    ni_old=$ni
    nj_old=${njs[$((n-1))]}
    n+=1
done