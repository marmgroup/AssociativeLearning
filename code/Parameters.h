//
//  Parameters.h
//
//  Created by Emiliano Méndez Salinas.
//  Copyright © 2019 Emiliano Méndez Salinas. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include <filesystem>


//===================Parameter related stuff============================================================



std::vector<int> string_to_vector(const std::string& str, char delimiter)
{
    std::vector<int> vec;
    std::string value;
    std::istringstream strStream(str);
    while (std::getline(strStream, value, delimiter))
    {
        vec.push_back(std::stoi(value));
    }
    return vec;
}


//// This class contains all the parameters that are read from an input file and available throughtout the code
//// to be used when appropiate.
class Parameters
{
public:
    const std::string SIMULATION_NAME; //// Any name should work
    const std::string MUT_DISTRIBUTION; //// Only "normal" or "cauchy" allowed for now
    const std::string INITIAL_WEIGHTS; //// Only "random" or "0" allowed for now
    const std::string INITIAL_LEARNING_RATE; //// Only "fixed" allowed, but parameter not used anymore
    const double LOW_BOUND; //// Any real number
    const double UPP_BOUND; //// Any real number
    const double MUTATION_PROBABILITY; //// Real number between 0 an 1
    const double MUTATION_DISPERSION; //// Real number between 0 and 1
    const std::string LEARNING_AGENT_TYPE; //// Only "Network", "Rescorla_Wagner_Rule", "Optimal_Performing_Rule" or "Averaging_Bayesian_Rule" allowed for now
    const int UPDATES_INPUT; //// Only "0" or "1" allowed
    const std::string THIRD_INPUT_TYPE; //// Only "simple_N","one_over_N","Delta_estimate" and "one_over_N_plus_two"
    const bool PROCESSING_BIAS; //// Only "0" or "1" allowed
    const std::string PROCESSING_TF; //// Only "tanh", "relu", "log", "none" allowed for now
    const bool OUTPUT_BIAS; //// Only "0" or "1" allowed
    const std::string OUTPUT_TF; //// Only "tanh", "relu", "log", "none" allowed for now
    const double SEED_MUTATIONS; //// Any number larger than 1
    const double SEED_ENVIRONMENTS; //// Any number larger than 1
    const double SEED_PLOTS; //// Any number larger than 1
    const int POPULATION_SIZE; //// Any discrete number larger than 0
    const int GENERATIONS; //// Any discrete number larger than 0
    
    const int UPDATES_NUMBER; //// Any discrete number larger than 0
    const int TESTS_NUMBER; //// Any discrete number larger than 0
    const double INITIAL_ESTIMATE; //// Any real number
    const std::string INITIAL_ESTIMATE_TYPE; //// Only "evolving", "fixed" or "random" allowed for now
    const std::string UPDATES_NUMBER_TYPE; //// Only "random" or "fixed" allowed for now
    const int MAX_UPDATES_NUMBER; //// Any discrete number larger than 1
    const double STOP_PROBABILITY; //// Real number between 0 an 1
    const std::string ERROR_CALCULATION; //// Only "squared" or "absolute" allowed for now
    const std::vector<int> NODES_PER_LAYER; //// One or more positive integers larger than 0 separated by commas
    
//// Constructor
    Parameters(
               const std::string & p_simulation_name,
               const std::string & p_mut_distribution,
               const std::string & p_initial_weights,
               const std::string & p_initial_learning_rate,
               const double & p_low_bound,
               const double & p_upp_bound,
               const double & p_mutation_probability,
               const double & p_mutation_dispersion,
               const std::string & p_learning_agent_type,
               const int & p_updates_input,
               const std::string & p_third_input_type,
               const bool & p_processing_bias,
               const std::string & p_processing_tf,
               const bool & p_output_bias,
               const std::string & p_output_tf,
               const double & p_seed_mutations,
               const double & p_seed_environments,
               const double & p_seed_plots,
               const int & p_population_size,
               const int & p_generations,
               
               const int & p_updates_number,
               const int & p_tests_number,
               const double & p_initial_estimate,
               const std::string & p_initial_estimate_type,
               const std::string & p_updates_number_type,
               const int & p_max_updates_number,
               const double & p_stop_probability,
               const std::string & p_error_calculation,
               
               
               const std::vector<int> & p_nodes_per_layer
               ) :
    SIMULATION_NAME(p_simulation_name),
    MUT_DISTRIBUTION(p_mut_distribution),
    INITIAL_WEIGHTS(p_initial_weights),
    INITIAL_LEARNING_RATE(p_initial_learning_rate),
    LOW_BOUND(p_low_bound),
    UPP_BOUND(p_upp_bound),
    MUTATION_PROBABILITY(p_mutation_probability),
    MUTATION_DISPERSION(p_mutation_dispersion),
    LEARNING_AGENT_TYPE(p_learning_agent_type),
    UPDATES_INPUT(p_updates_input),
    THIRD_INPUT_TYPE(p_third_input_type),
    PROCESSING_BIAS(p_processing_bias),
    PROCESSING_TF(p_processing_tf),
    OUTPUT_BIAS(p_output_bias),
    OUTPUT_TF(p_output_tf),
    SEED_MUTATIONS(p_seed_mutations),
    SEED_ENVIRONMENTS(p_seed_environments),
    SEED_PLOTS(p_seed_plots),
    POPULATION_SIZE(p_population_size),
    GENERATIONS(p_generations),
    
    UPDATES_NUMBER(p_updates_number),
    TESTS_NUMBER(p_tests_number),
    INITIAL_ESTIMATE(p_initial_estimate),
    INITIAL_ESTIMATE_TYPE(p_initial_estimate_type),
    UPDATES_NUMBER_TYPE(p_updates_number_type),
    MAX_UPDATES_NUMBER(p_max_updates_number),
    STOP_PROBABILITY(p_stop_probability),
    ERROR_CALCULATION(p_error_calculation),
    NODES_PER_LAYER(p_nodes_per_layer)
    {
        
    }
};

//// Reads the parameter file and returns a Parameter object constructed from the values contained in the
//// parameter file. This is the object that will be used throughout the simulation to have access to the
//// parameter values wherever neeeded.
Parameters readParameters(std::ifstream &Parameters_ifs)
{
    
    if(!Parameters_ifs.is_open())
    {
        std::cout << "ERROR: UNABLE TO OPEN PARAMETERS FILE" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //// Variables that will store the values read from the parametr file.
    std::string p_SIMULATION_NAME;
    std::string p_MUT_DISTRIBUTION;
    std::string p_INITIAL_WEIGHTS;
    std::string p_INITIAL_LEARNING_RATE;
    double p_LOW_BOUND;
    double p_UPP_BOUND;
    double p_MUTATION_PROBABILITY;
    double p_MUTATION_DISPERSION;
    std::string p_LEARNING_AGENT_TYPE;
    int p_UPDATES_INPUT;
    std::string p_THIRD_INPUT_TYPE = "empty"; //// The only parameter that is pre-initialized in case it is not
    bool p_PROCESSING_BIAS;                 //// found in the (old type) parameter file. This should be removed
    std::string p_PROCESSING_TF;            //// for the final version of the program.
    bool p_OUTPUT_BIAS;
    std::string p_OUTPUT_TF;
    double p_SEED_MUTATIONS;
    double p_SEED_ENVIRONMENTS;
    double p_SEED_PLOTS;
    int p_POPULATION_SIZE;
    int p_GENERATIONS;
    
    int p_UPDATES_NUMBER;
    int p_TESTS_NUMBER;
    double p_INITIAL_ESTIMATE;
    std::string p_INITIAL_ESTIMATE_TYPE;
    std::string p_UPDATES_NUMBER_TYPE;
    int p_MAX_UPDATES_NUMBER;
    double p_STOP_PROBABILITY;
    std::string p_ERROR_CALCULATION;
    
    std::string p_NODES_PER_LAYER;
    
    //// The parameter file is read and the values of the parameters in the file are assigned to the respective
    //// variable. Order should not matter here.
    for(;;)
    {
        std::string parameter_line;
        Parameters_ifs >> parameter_line;
        if(Parameters_ifs.good()){
            
            if(parameter_line == "SIMULATION_NAME")
                Parameters_ifs >> p_SIMULATION_NAME;
            if(parameter_line == "MUT_DISTRIBUTION")
                Parameters_ifs >> p_MUT_DISTRIBUTION;
            if(parameter_line == "INITIAL_WEIGHTS")
                Parameters_ifs >> p_INITIAL_WEIGHTS;
            if(parameter_line == "INITIAL_LEARNING_RATE")
                Parameters_ifs >> p_INITIAL_LEARNING_RATE;
            if(parameter_line == "LOW_BOUND")
                Parameters_ifs >> p_LOW_BOUND;
            if(parameter_line == "UPP_BOUND")
                Parameters_ifs >> p_UPP_BOUND;
            if(parameter_line == "MUTATION_PROBABILITY")
                Parameters_ifs >> p_MUTATION_PROBABILITY;
            if(parameter_line == "MUTATION_DISPERSION")
                Parameters_ifs >> p_MUTATION_DISPERSION;
            if(parameter_line == "LEARNING_AGENT_TYPE")
                Parameters_ifs >> p_LEARNING_AGENT_TYPE;
            if(parameter_line == "UPDATES_INPUT")
                Parameters_ifs >> p_UPDATES_INPUT;
            if(parameter_line == "THIRD_INPUT_TYPE")
                Parameters_ifs >> p_THIRD_INPUT_TYPE;
            if(parameter_line == "PROCESSING_BIAS")
                Parameters_ifs >> p_PROCESSING_BIAS;
            if(parameter_line == "PROCESSING_TF")
                Parameters_ifs >> p_PROCESSING_TF;
            if(parameter_line == "OUTPUT_BIAS")
                Parameters_ifs >> p_OUTPUT_BIAS;
            if(parameter_line == "OUTPUT_TF")
                Parameters_ifs >> p_OUTPUT_TF;
            if(parameter_line == "SEED_MUTATIONS")
                Parameters_ifs >> p_SEED_MUTATIONS;
            if(parameter_line == "SEED_ENVIRONMENTS")
                Parameters_ifs >> p_SEED_ENVIRONMENTS;
            if(parameter_line == "SEED_PLOTS")
                Parameters_ifs >> p_SEED_PLOTS;
            if(parameter_line == "POPULATION_SIZE")
                Parameters_ifs >> p_POPULATION_SIZE;
            if(parameter_line == "GENERATIONS")
                Parameters_ifs >> p_GENERATIONS;
            
            if(parameter_line == "UPDATES_NUMBER")
                Parameters_ifs >> p_UPDATES_NUMBER;
            if(parameter_line == "TESTS_NUMBER")
                Parameters_ifs >> p_TESTS_NUMBER;
            if(parameter_line == "INITIAL_ESTIMATE")
                Parameters_ifs >> p_INITIAL_ESTIMATE;
            if(parameter_line == "INITIAL_ESTIMATE_TYPE")
                Parameters_ifs >> p_INITIAL_ESTIMATE_TYPE;
            if(parameter_line == "UPDATES_NUMBER_TYPE")
                Parameters_ifs >> p_UPDATES_NUMBER_TYPE;
            if(parameter_line == "MAX_UPDATES_NUMBER")
                Parameters_ifs >> p_MAX_UPDATES_NUMBER;
            if(parameter_line == "STOP_PROBABILITY")
                Parameters_ifs >> p_STOP_PROBABILITY;
            if(parameter_line == "ERROR_CALCULATION")
                Parameters_ifs >> p_ERROR_CALCULATION;
            
            
            if(parameter_line == "NODES_PER_LAYER")
                Parameters_ifs >> p_NODES_PER_LAYER;
            
        }
        else
            break;
    }
    
    Parameters_ifs.close();
    
    
    //// Instantiation and return of the Parameter object using the values read from the parameter file and stored
    //// in the variables used as arguments in the constructor.
    return Parameters ( p_SIMULATION_NAME,
                       p_MUT_DISTRIBUTION,
                       p_INITIAL_WEIGHTS,
                       p_INITIAL_LEARNING_RATE,
                       p_LOW_BOUND,
                       p_UPP_BOUND,
                       p_MUTATION_PROBABILITY,
                       p_MUTATION_DISPERSION,
                       p_LEARNING_AGENT_TYPE,
                       p_UPDATES_INPUT,
                       p_THIRD_INPUT_TYPE,
                       p_PROCESSING_BIAS,
                       p_PROCESSING_TF,
                       p_OUTPUT_BIAS,
                       p_OUTPUT_TF,
                       p_SEED_MUTATIONS,
                       p_SEED_ENVIRONMENTS,
                       p_SEED_PLOTS,
                       p_POPULATION_SIZE,
                       p_GENERATIONS,
                       
                       p_UPDATES_NUMBER,
                       p_TESTS_NUMBER,
                       p_INITIAL_ESTIMATE,
                       p_INITIAL_ESTIMATE_TYPE,
                       p_UPDATES_NUMBER_TYPE,
                       p_MAX_UPDATES_NUMBER,
                       p_STOP_PROBABILITY,
                       p_ERROR_CALCULATION,
                       
                       string_to_vector(p_NODES_PER_LAYER,',') );
}








//// Function to check if a file or folder exists. Used to select the folder where the output files will be saved
//// and to stop simulations with the same name as another simulation found in the same folder, that is, to prevent
//// overwritting of output files.
bool folder_or_file_Exists(const std::string& folder_or_file) {
    struct stat buf;
    return (stat(folder_or_file.c_str(), &buf) == 0);

    //// Suggested by Hanno but not working in my OS. Maybe add it later.
    //return std::filesystem::exists(folder_or_file);
}


//  You can specify where to save the simulations, otherwise it is saved in a current directory
std::string get_folder()
{

    if(folder_or_file_Exists("C:/simulation/"))
        return "C:/simulation/";
    else
    {
        return (std::filesystem::current_path().string() + "/");
    }
}


//// Imports the parameter file, which might or not be given as an argument to main function, and saves it a
//// ifsgream object.
std::ifstream get_input_file(int & argc, char* argv[])
{
    if(argc > 1)
    {
        //// The name of the parameter file can be arbitrary and should be located in the same path and provided
        //// together, after the executable in the command line.
        std::ifstream Parameters_ifs(argv[1]);
        return Parameters_ifs;
    }
    else
    {
        //// If not running from the command line or if no argument for parameter file is provided in the command
        //// line, the default file must be named "myParams.txt" and be located in the output folder.
        std::ifstream Parameters_ifs(get_folder() + "myParams.txt");
        return Parameters_ifs;
    }
    
}




//// Actually duplicates a text file.
void save_as_txt(const std::string & inputFileName, const std::string & outputFileName)
{
    std::ifstream infile(inputFileName);
    std::ofstream outfile(outputFileName);
    std::string content = "";
    int i;
    
    for (i = 0; infile.eof() != true; i++) //// Get content of infile
        content += infile.get();
    
    i--;
    content.erase(content.end() - 1);     //// Erase last character
    
    infile.close();
    
    outfile << content;                 //// Output
    outfile.close();
}


//// A copy of the parameter file is generated and saved together with the other ouput files. This is mostly useful
//// to always save the parameters together with the other output files so that it is easy to consult what were the
//// simulation parameters that gave such results. Also allows to have a folder with parameter files that can be
//// used for several simulations without the need to manually copy or move them. Note that this function takes
//// a Parameters object as function argument, which means that this happens when the parameter file has already
//// been read and used to create a Parameters object. This can be clearly seen in main.cpp
void copy_and_save_parameter_file(int & argc, char* argv[], const Parameters &Par)
{
    
    //std::string Run_time_data_file = get_folder() + Par.SIMULATION_NAME + "__Run_time_data.txt";
    
    //// Create the names for the "_Parameters.txt" and "__RW_Plot.txt" output files, which are later used to check
    //// whether results of a simulation with the same name are saved in the same directory, in which case the
    //// new simulation is force stopped sending an error messsage.
    std::string Parameter_file_name = get_folder() + Par.SIMULATION_NAME + "_Parameters.txt";
    std::string RWPlot_file_name = get_folder() + Par.SIMULATION_NAME + "__RW_Plot.txt";
    if(folder_or_file_Exists(Parameter_file_name) )
    {
        std::cout << "ERROR: A SIMULATION WITH THE SAME NAME IS SAVED IN THE SAME DIRECTORY. PLEASE USE A DIFFERENT NAME" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(folder_or_file_Exists(RWPlot_file_name) )
    {
        std::cout << "ERROR: A RW_PLOT FILE WITH THE SAME NAME IS SAVED IN THE SAME DIRECTORY. PLEASE CHECK" << std::endl;
        exit(EXIT_FAILURE);
    }
    //// If there is no danger of overwritting simulation output, create and save the copy of the parameter file.
    else
    {
        //// If parameter file is given as argument to the main function, the value of the parameter
        //// "SIMULATION_NAME" is used as the name for the copy of the parameter file. Since seldom the name
        //// actual name of the file and the value of the parameter will be the same, both files should be in
        //// different folders, for which is useful to have a folder containing exclusively executable(s) and
        //// "original" (the copy is also original in the sense that they are identical files) parameter files
        //// in a folder that is different to that where output files (including copied parameter files) will
        //// be saved.
        if(argc > 1)
        {
            std::string New_File_Name = argv[1];
            save_as_txt(New_File_Name, Parameter_file_name);
        }
        //// Otherwise, the parameters are taken from the standardized "myParams.txt" file that should be
        //// located in the output folder. Next a copy of this same file is saved with a different name, that
        //// contained in the parameter file itself as the parameter "SIMULATION_NAME", which is, by the way, also
        //// already "stored" in the Parameters object.
        else
        {
            save_as_txt(get_folder() + "myParams.txt", Parameter_file_name);
        }
    }
}





#endif /* Parameters_h */
