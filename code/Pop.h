//
//  Pop.h
//
//  Created by Emiliano Méndez Salinas.
//  Copyright © 2020 Emiliano Méndez Salinas. All rights reserved.
//

#ifndef Pop_h
#define Pop_h





//// Stop simulation and print error if either processing or output transfer functions specified in the
//// parameters are not valid.
void output_tf_error()
{
    std::cout << "ERROR: UNKNOWN PARAMETER OUTPUT_TF" << std::endl;
    exit(EXIT_FAILURE);
}

void processing_tf_error()
{
    std::cout << "ERROR: UNKNOWN PARAMETER PROCESSING_TF" << std::endl;
    exit(EXIT_FAILURE);
}





//// As seen in next function, this one is called only when the learning agent is a type of Network. It instantiates
//// the four static member variables of class Network and then creates the population of networks and calls the
//// do_stuff function which is the second outermost wrapper for most of the things that happen in the simulation.
//// Read forward for a clearer explanation.

template <class Processing_layer, class Output_node>
void create_population(int & argc, char* argv[], const Parameters &Par)
{
    //// The following four static variables and functions called to initialize them are explained in Networks.h
    Network<Processing_layer,Output_node>::weights_number = calculate_weights_number(Par);
    Network<Processing_layer,Output_node>::temp_values_size = calculate_temp_values_size(Par);
    
    Network<Processing_layer,Output_node>::old_values =
    std::vector<double> (Network<Processing_layer,Output_node>::temp_values_size,0);
    Network<Processing_layer,Output_node>::new_values =
    std::vector<double> (Network<Processing_layer,Output_node>::temp_values_size,0);
    
    //// The corresponding population is created according to the template arguments passed to the create_population
    //// template function, and then the do_stuff function is called to continue with the simulation.
    std::vector< learning_Agent< Network<Processing_layer, Output_node> > >Population;
    for (int i = 0; i < Par.POPULATION_SIZE; ++i)
    {
        Population.push_back(learning_Agent<Network<Processing_layer, Output_node> >(Par));
    }
    do_stuff<learning_Agent<Network<Processing_layer, Output_node>>, std::vector< learning_Agent<Network<Processing_layer, Output_node> > > >(argc,Population, Par);
}



//// Other than functions related to importing and managing parameters and measuring run time, this is the only
//// function that is called in the main function in main.cpp. This means that most things of the simulation will
//// happen within this function, by calling many other functions, so it works as the outermost wrapper for most
//// of the action. Its specific function is to create the appropiate type of populations of learning agents
//// depending on the parameters, and by been made a separate function, it is more organized and concise.
void start_everything(int & argc, char* argv[], const Parameters &Par)
{
    //// When the learning agent is a Network, types for transfer functions of processing layer and output layer
    //// have to be specified (acording to parameter values) and with that information given as template arguments,
    //// for the sake of more ordered code, another function has to be called where the objects that constitute the
    //// population are actually created.
    if (Par.LEARNING_AGENT_TYPE == "Network")
    {
        if(Par.PROCESSING_TF == "tanh")
        {
            Par.OUTPUT_TF == "none" ? create_population<tanh_tf,none>(argc,argv,Par) :
            Par.OUTPUT_TF == "log" ? create_population<tanh_tf,log_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "tanh" ? create_population<tanh_tf,tanh_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "relu" ? create_population<tanh_tf,relu_tf>(argc,argv,Par) :
            output_tf_error();
        }
        else if(Par.PROCESSING_TF == "none")
        {
            Par.OUTPUT_TF == "none" ? create_population<none,none>(argc,argv,Par) :
            Par.OUTPUT_TF == "log" ? create_population<none,log_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "tanh" ? create_population<none,tanh_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "relu" ? create_population<none,relu_tf>(argc,argv,Par) :
            output_tf_error();
        }
        else if(Par.PROCESSING_TF == "log")
        {
            Par.OUTPUT_TF == "none" ? create_population<log_tf,none>(argc,argv,Par) :
            Par.OUTPUT_TF == "log" ? create_population<log_tf,log_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "tanh" ? create_population<log_tf,tanh_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "relu" ? create_population<log_tf,relu_tf>(argc,argv,Par) :
            output_tf_error();
        }
        else if(Par.PROCESSING_TF == "relu")
        {
            Par.OUTPUT_TF == "none" ? create_population<relu_tf,none>(argc,argv,Par) :
            Par.OUTPUT_TF == "log" ? create_population<relu_tf,log_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "tanh" ? create_population<relu_tf,tanh_tf>(argc,argv,Par) :
            Par.OUTPUT_TF == "relu" ? create_population<relu_tf,relu_tf>(argc,argv,Par) :
            output_tf_error();
        }
        else
            processing_tf_error();
    }
    //// When the learning agent is any of the rules, the population is created right away according to the type of
    //// rule specified in the parameters and the "do_stuff" function, the second outermost wrapper for most of the
    //// action, is called.
    else if (Par.LEARNING_AGENT_TYPE == "Rescorla_Wagner_Rule")
    {
        std::vector<learning_Agent<Rescorla_Wagner_Rule> > Population;
        
        for(int i = 0; i < Par.POPULATION_SIZE; ++i)
        {
            Population.push_back(Par);
        }
        do_stuff<learning_Agent<Rescorla_Wagner_Rule>,std::vector<learning_Agent<Rescorla_Wagner_Rule> > >(argc,Population,Par);
    }
    else if (Par.LEARNING_AGENT_TYPE == "Optimal_Performing_Rule")
    {
        std::vector<learning_Agent<Optimal_Performing_Rule> > Population;
        for(int i = 0; i < Par.POPULATION_SIZE; ++i)
        {
            Population.push_back(Par);
        }
        do_stuff<learning_Agent<Optimal_Performing_Rule>,std::vector<learning_Agent<Optimal_Performing_Rule> > >(argc,Population,Par);
    }
    else if (Par.LEARNING_AGENT_TYPE == "Averaging_Bayesian_Rule")
    {
        std::vector<learning_Agent<Averaging_Bayesian_Rule> > Population;
        for(int i = 0; i < Par.POPULATION_SIZE; ++i)
        {
            Population.push_back(Par);
        }
        do_stuff<learning_Agent<Averaging_Bayesian_Rule>,std::vector<learning_Agent<Averaging_Bayesian_Rule> > >(argc,Population,Par);
    }
    else
    {
        std::cout << "ERROR: UNKNOWN PARAMETER LEARNING_AGENT_TYPE" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    
    
    
}













#endif /* Pop_h */
