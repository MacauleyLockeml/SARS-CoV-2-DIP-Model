Description of the files in this folder:

parameters.txt - the reference set of the calibrated model parameters

det_model.jl - numerical solution of the calibrated deterministic model with varying initial conditions
stoch_model.jl - numerical solution of the stochastic MC-based model (simulation of one trajectory)

run_timecourse.jl - simulation of an ensemble of stochastic trajectories, saving the evolution in time of the ensemble statistics for each variable to csv-files
run_densities.jl - simulation of an ensemble of stochastic trajectories, saving the histograms of the selected variables at the end timepoint
run_dose_response.jl - simulation of an ensemble of stochastic trajectories, saving the mean and median estimates of the number of WT virions and DIPs produced present 24 hours post infection, as well as the probability of productive infection

draw_timecourse.jl, draw_densities.jl, draw_dose_response.jl scripts can be used to visualized the results of running the scripts above
