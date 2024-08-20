# TradeOff package
Starting a new project on how trade-offs effect ecosystem assembly. Including work on how immgration effects carbon use efficiency.

## Installing the TradeOff packaged in Julia
1. Enter Pkg mode in Julia and run "activate . "
2. Run "instantiate"
3. Run "dev . "

These three commands lets Julia know that a non-registered packaged is being used.

## Full procedure used to produce all data and figures for Immigration thesis.

## Full procedure for generating the data for the paper

## Full summary of the scripts in this repository
### vars_over_t.jl
This script takes in -- as input, averages over the community and generates a file per replicate simulation.

### traj_stats_.jl
This script uses the output files of vars_over_t.jl, calculates averages across replicate simulations and produces one file.