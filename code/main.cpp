//
//  main.cpp
//
//  Created by Emiliano Méndez Salinas.
//  Copyright © 2019 Emiliano Méndez Salinas. All rights reserved.
//  Code for the paper "Neural network models for the evolution of associative learning"
//  by Emiliano Méndez Salinas, Franz J. Weissing and Magdalena Kozielska
//

#include <iostream>
#include <vector>
#include <array>


#include <cmath>
#include <random>
#include <sys/stat.h>
#include <chrono>

#include <sstream>
#include <fstream>




#include "Parameters.h"
#include "Networks.h"
#include "Simulation.h"
#include "Pop.h"




//// Declaration of random number engines and distributions that are used throughout the code. Notice that
//// mutational and environmental processes have different random numbers. This allows for example to have
//// simulations where the environment is exactly the same throughout evolution but the population is initialized
//// differently and undergoes different mutations throughout evolution. The third random number engine is only used
//// for the random values needed for the RW_Plots and thus allows full replicability of the other two regardless of
//// whether or not and how the rescorla_wagner_plot function is called.
std::mt19937_64 rngMutations;
std::mt19937_64 rngEnvironment;
std::mt19937_64 rngPlots;
std::uniform_real_distribution<double> Uniform0to1(0.0, 1.0);
std::uniform_real_distribution<double> UniformNeg1to1(-1.0, 1.0);






//// If called (when type of third input specified in the parameters file is not valid), stops simulation and
//// prints error message.
void third_input_type_error()
{
    std::cout << "ERROR: UNKNOWN PARAMETER PROCESSING_TF" << std::endl;
    exit(EXIT_FAILURE);
}


//// Second outermost wrapper for most of the things that happen in the simulation. It has the function of selecting
//// the type of third input to be used during the simulation (based on Parameters) and call the function
//// void main_loop (the one where everything really happens) with the appropiate template argument. Similarly,
//// after execution of main_loop function finishes, it also calls the function void rescorla_wagner_plot with the
//// same set of template arguments.
template<class N, class P>
void do_stuff(int & argc, P& Pop, const Parameters& Par)
{
    
    if(argc <= 2){
        Par.LEARNING_AGENT_TYPE == "Averaging_Bayesian_Rule" ? main_loop<N, P, simple_N>(Pop, Par) :
        Par.THIRD_INPUT_TYPE == "one_over_N" ? main_loop<N, P, one_over_N>(Pop, Par) :
        Par.THIRD_INPUT_TYPE == "Delta_estimate" ? main_loop<N, P, Delta_estimate>(Pop, Par):
        Par.THIRD_INPUT_TYPE == "simple_N" ? main_loop<N, P, simple_N>(Pop, Par):
        Par.THIRD_INPUT_TYPE == "one_over_N_plus_two" ? main_loop<N, P, one_over_N_plus_two>(Pop, Par):
        third_input_type_error();
    }
    std::ofstream RW_Plot(get_folder() + Par.SIMULATION_NAME + "__RW_Plot.txt");
    RW_Plot << "Generation,Individual,Difference,Delta" << '\n';
    
    const size_t RW_PLOT_REPS = 100;
    
    Par.LEARNING_AGENT_TYPE == "Averaging_Bayesian_Rule" ? rescorla_wagner_plot<N, P, simple_N>(Pop, Par, Par.GENERATIONS,RW_PLOT_REPS, RW_Plot) :
    Par.THIRD_INPUT_TYPE == "one_over_N" ? rescorla_wagner_plot<N, P, one_over_N>(Pop, Par, Par.GENERATIONS,RW_PLOT_REPS, RW_Plot) :
    Par.THIRD_INPUT_TYPE == "Delta_estimate" ? rescorla_wagner_plot<N, P, Delta_estimate>(Pop, Par, Par.GENERATIONS,RW_PLOT_REPS, RW_Plot):
    Par.THIRD_INPUT_TYPE == "simple_N" ? rescorla_wagner_plot<N, P, simple_N>(Pop, Par, Par.GENERATIONS,RW_PLOT_REPS, RW_Plot):
    Par.THIRD_INPUT_TYPE == "one_over_N_plus_two" ? rescorla_wagner_plot<N, P, one_over_N_plus_two>(Pop, Par, Par.GENERATIONS,RW_PLOT_REPS, RW_Plot):
    third_input_type_error();

    
}




//// Prints simulation runtime information to an output file.
void running_time_data(std::ofstream & ofs, const std::chrono::milliseconds & time)

{
    if(ofs.is_open())
    {
        ofs << '\n' << " Execution time:  "
        << time.count() << "  milliseconds.\n\n                = "
        << ( time.count() / 1000.0 ) / 60.0 << "  minutes.\n\n                = "
        << ( time.count() / 1000.0 ) / 3600 << "  hours.";
    }
    ofs.close();
}


int main(int argc, char* argv[]) {
    
    //// Create clock object to store simulation begin time.
    auto Sim_beginning = std::chrono::high_resolution_clock::now();
    
    
    
    //// Reads the parameter file in.
    std::ifstream Parameters_ifs = get_input_file(argc, argv);
    
    
    //// Creates a Parameter object from the read parameter file.
    const Parameters Par = readParameters(Parameters_ifs);
    
    //// Saves a copy of the parameter file as explained in the definition of the function in Parameters.h
    copy_and_save_parameter_file(argc, argv, Par);
    
    
    //// IMPORTANT rngPlots should be seeded here too!!!
    rngMutations.seed(Par.SEED_MUTATIONS);
    rngEnvironment.seed(Par.SEED_ENVIRONMENTS);
    
    
    //// As explained in Pop.h, outermost wrapper for most of the things that should happen in the simulation.
    //// Next outermost wrapper is create_population (only used when learning_agents are Networks) and afterwards
    //// do_stuff, where the explicit code of what happens lies.
    start_everything(argc, argv, Par);
    
    //// Create clock object to store simulation end time.
    auto Sim_ending = std::chrono::high_resolution_clock::now();
    
    //// Estimate simulation run time using the two clocks.
    auto Sim_execution_time =
    std::chrono::duration_cast<std::chrono::milliseconds>(Sim_ending - Sim_beginning);
    
    //// Create output file for simulation runtime information and print to it.
    std::ofstream Run_time_data(get_folder() + Par.SIMULATION_NAME + "__Run_time_data.txt");
    running_time_data(Run_time_data, Sim_execution_time);
    
    
    
    return 0;
}




