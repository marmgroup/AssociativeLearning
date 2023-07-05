//
//  Simulation.h
//
//  Created by Emiliano Méndez Salinas.
//  Copyright © 2019 Emiliano Méndez Salinas. All rights reserved.
//

#ifndef Simulation_h
#define Simulation_h

extern std::uniform_real_distribution<double> Uniform0to1;
extern std::mt19937_64 rngEnvironment;
extern std::mt19937_64 rngPlots;








//// All individuals that take part of the simulation are learning agents. The argument to this template class
//// determines which kind of learning agent they are and thus what kind of learning mechanism they possess and how
//// the update of the "state of the world" takes place. All kinds of organism have the same member functions, but
//// each with its particular implementation, so that they can all be called from the container learning_Agent.
//// One kind of learning_Agent an artificial nerual network (Network), which is itself a template class and can
//// be of different kind depending on the type of transfer functions it has. Network class is the only one
//// defined in its own Networks.h file. The other kinds are all rules defined later in this file. learning_Agent
//// class contains all the variables and functions needed for the organisms to interact with the environment and
//// to reproduce. As their name implies, learning_Agent objects have to undergo several learning experiences
//// throughout their life, in which they have to estimate de probability of getting a reward (1) or no reward (0)
//// in the next experience, based on the outcome (1 or 0) of previous experiences. In this way, every new
//// experience, they will update their estimate of the probability of getting a reward next time. The fitness of
//// this individuals is (inversely) determined by the cumulative error of their estimates with respect to the true
//// probabilities in the environment and they reproduce in proportion to their fitness.
template<class T>
class learning_Agent
{
    double initial_estimate; //// The "innate" expected probability of reward, prior to any experience.
    double estimate; //// Recent most expected probability of reward.
    double error; //// Cumulative difference of the estimated probability with respect to the true probability.
    double Delta_estimate; //// Difference between recent most and immediately previous estimate.
    T learn_mechanism; //// A Network or a lerning rule as defined elsewhere. Each with their own way of learning.
public:
    //// Constructor. Apart from itself calling the constructor of the learn_mechanism, initializes the other four
    //// member variables of the class to defaulr and Parameter values.
    learning_Agent(const Parameters &Par) : learn_mechanism(Par), initial_estimate(Par.INITIAL_ESTIMATE), estimate(0.0), error(0.0), Delta_estimate(1.0) {}
    
    //// Calls the homonym function in learn_mechanism and is used only during reproduction to decide when and where
    //// mutations happen and what the new values of the genes are.
    void mutate_all_genome(const Parameters &Par)
    {
        learn_mechanism.mutate_all_genome(Par);
    }
    
    //// Sets values of all member variables back to initial conditions, with exception of initial_estimate that
    //// takes its new value from an argument to the function. This will be useful for the scenarios where this
    //// member variable is subject to evolution.
    void reset_values(const double & init_estimate)
    {
        initial_estimate = init_estimate;
        estimate = 0.0;
        error = 0.0;
        Delta_estimate = 1.0;
    }
    
    //// Calls the activate function of learn_mechanism which is where the learning takes place according to the
    //// specific implementation of each mechanism and then updates the values of the estimate and Delta_estimate.
    void update_estimate2(const Input & input, const Parameters & Par) {
        double new_estimate = learn_mechanism.activate(input,Par);
        Delta_estimate = fabs(estimate - new_estimate);
        estimate = new_estimate;
    }
    
    //// Both commented later in their definition.
    void select_first_estimate(const Parameters &Par);
    void update_error(const double & val, const Parameters &Par);

    //// Return the values of the four member variables. 
    double get_initial_estimate() {return initial_estimate;}
    double get_estimate() {return estimate;}
    double get_error() {return error;}
    double get_Delta_estimate() { return Delta_estimate; }
    
    //double set_estimate(const double &val) {estimate = val;}
    
    //// Call the homonym functions in learn_mechanism to export genome of the learning_Agent or the labels that
    //// correspond to and identify each element of the genome.
    void print_labels(std::ofstream &ofs, const Parameters & Par) {learn_mechanism.print_labels(ofs,Par);}
    void print_values(std::ofstream &ofs) {learn_mechanism.print_values(ofs);}
    
};




//// Implements the "Rescorla_Wagner_Rule" as described in Trimmer et al. 2012.
class Rescorla_Wagner_Rule
{
private:
    double learning_rate; //// The only (evolvable) element of the genome of this class.
public:
    //// Constructor, initializes the only member variable of the class.
    Rescorla_Wagner_Rule(const Parameters &Par) : learning_rate(value_after_mutation(Par,0.0) ) {}
    
    //// Only needs to check if learning_rate should mutate and if so, what should its new value be.
    void mutate_all_genome(const Parameters & Par)
    {
        learning_rate = value_after_mutation(Par, learning_rate);
    }
    
    //// Solves the learning rule's equation
    double activate(const Input & input, const Parameters & Par)
    {
        return input.estimate() + learning_rate * (input.reward() - input.estimate() ) ;
    }
    
    //// Export the only member of the class and its respective label.
    void print_values(std::ofstream & ofs) const { ofs << "," << learning_rate; }
    void print_labels(std::ofstream & ofs, const Parameters & Par) const { ofs << "," << "learning_rate"; }
    
};

//// Implements the "Optimal_Performing_Rule" as described in Trimmer et al. 2012.
class Optimal_Performing_Rule
{
private:
    double learning_rate; //// The only (evolvable) element of the genome of this class.
public:
    //// Constructor, initializes the only member variable of the class.
    Optimal_Performing_Rule(const Parameters &Par) : learning_rate(value_after_mutation(Par,0.0) ) {}
    
    //// Only needs to check if learning_rate should mutate and if so, what should its new value be.
    void mutate_all_genome(const Parameters & Par)
    {
        learning_rate = value_after_mutation(Par, learning_rate);
    }
    
    //// Solves the learning rule's equation
    double activate(const Input & input, const Parameters & Par)
    {
        return input.estimate() + learning_rate * (input.reward() - 0.5 );
    }
    
    //// Export the only member of the class and its respective label.
    void print_values(std::ofstream & ofs) const { ofs << "," << learning_rate; }
    void print_labels(std::ofstream & ofs, const Parameters & Par) const { ofs << "," << "learning_rate"; }
    
};

//// Implements the "Averaging_Bayesian_Rule" as described in Trimmer et al. 2012.
class Averaging_Bayesian_Rule
{
public:
    //// Nothing to initialize in constructor.
    Averaging_Bayesian_Rule(const Parameters &Par) {}
    
    //// Nothing to mutate so function is empty, but still present because eventually it has to be called
    //// and if not included the program would not compile.
    void mutate_all_genome(const Parameters & Par) {}

    //// Solves the learning rule's equation
    double activate(const Input & input, const Parameters & Par)
    {
        return input.estimate() + (1.0 /  ( input.number() + 2.0) ) * (input.reward() - input.estimate() );
    }
    
    //// No member variables to export but functions still present because eventually they have to
    //// be called and if not included the program would not compile.

    void print_values(std::ofstream & ofs) const {}
    void print_labels(std::ofstream & ofs, const Parameters & Par) const {}

};



//// Allows to calculate the error (difference between estimated and real probability) in two different ways,
//// as squared error or as absolute error, depending on the value of the parameter, and add it to the cumulative
//// error.
template<class T>
void learning_Agent<T>::update_error(const double & val, const Parameters &Par)
{
    if(Par.ERROR_CALCULATION == "squared")
    {
        error += ( (val * val) / static_cast<double>(Par.TESTS_NUMBER) );
    }
    else if(Par.ERROR_CALCULATION == "absolute")
    {
        error += ( fabs(val) / static_cast<double>(Par.TESTS_NUMBER) );
    }
    else
    {
        std::cout << "ERROR: UNKNOWN PARAMETER ERROR_CALCULATION" << std::endl;
        exit(EXIT_FAILURE);
    }
}


//// Assigns an initial value to the estimate. This will be done every time the agent "arrives" to a new environment
//// (a new probability to be estimated), prior to any experience. This initial value can be a random number, or
//// the value of the member variable initial_estimate, depending on the value of the respective Parameter. Also an
//// initial value for Delt_estimate is assigned.
template<class T>
void learning_Agent<T>::select_first_estimate(const Parameters &Par)
{
    Delta_estimate = 1.0;
    if(Par.INITIAL_ESTIMATE_TYPE == "random")
    {
        estimate = Uniform0to1(rngEnvironment);
    }
    else if(Par.INITIAL_ESTIMATE_TYPE == "fixed" || Par.INITIAL_ESTIMATE_TYPE == "evolving")
    {
        estimate = initial_estimate;
    }
    else
    {
        std::cout << "ERROR: UNKNOWN PARAMETER INITIAL_ESTIMATE_TYPE" << std::endl;
        exit(EXIT_FAILURE);
    }
}


//// Sets the number of experiences the learning_Agent will have to learn and update its estimate of the probability
//// before it is compared agains the real one to calculate the error. Depending on the respective Parameters, it
//// can be a fixed and predetermined number of times or a random number of times (within certain bounds) where
//// each time elapsed has a certain probability to be the last one.
int calculate_updates_number(const Parameters &Par)
{
    if(Par.UPDATES_NUMBER_TYPE == "fixed")
    {
        return Par.UPDATES_NUMBER;
    }
    else if (Par.UPDATES_NUMBER_TYPE == "random")
    {
        for(int i = 1; i <= Par.MAX_UPDATES_NUMBER; ++i)
        {
            if(Uniform0to1(rngEnvironment) < Par.STOP_PROBABILITY)
            {
                return i;
            }
        }
    }
    else
    {
        std::cout << "ERROR: UNKNOWN PARAMETER UPDATES_NUMBER_TYPE" << std::endl;
        exit(EXIT_FAILURE);
    }
    return Par.MAX_UPDATES_NUMBER;
}





//// Useful for the scenarios where the initial estimate is evolving. In such cases, the offspring individual will
//// inherit the (possibly) mutated initial estimate of the parental individual. So this function will be called
//// when reproduction is taking place. If the initial estimate is fixed it is always the same and it is a value
//// predetermined in the Parameters. If it is random, the value returned does not really matter beacause it will
//// be reassigned anyway in each environment by the function select_first_estimate defined earlier.
double transmit_initial_estimate(const double &val, const Parameters &Par)
{
    if(Par.INITIAL_ESTIMATE_TYPE == "evolving")
        return value_after_mutation(Par, val);
    else if(Par.INITIAL_ESTIMATE_TYPE == "fixed" || Par.INITIAL_ESTIMATE_TYPE == "random")
        return Par.INITIAL_ESTIMATE;
    else
    {
        std::cout << "ERROR: UNKNOWN PARAMETER INITIAL_ESTIMATE_TYPE" << std::endl;
        exit(EXIT_FAILURE);
    }
}



//// Creates a new offspring population, of the same size, out of the parental population. Each offspring individual
//// inherits its genome and initial estimate from a parental individual randomly selected with a probability that
//// is proportional to its fitness. This function estimates fitness values for each individual (from their
//// life time cumulative error), selects organisms that will leave offspring depending on fitness values, and
//// copies the genomes and initial estimates (with mutations) from parents to offspring.
template<class N, class P>
void reproduce(P & Population, const Parameters & Par)
{
    //// Create the container for the offspring population
    P next_generation;
    next_generation.reserve(Par.POPULATION_SIZE);
    
    //// Check if at least one individual has different error.
    bool all_equal = true;
    for(int i = 1; i < Par.POPULATION_SIZE; ++i)
    {
        if(Population[i - 1].get_error() != Population[i].get_error() )
        {
            all_equal = false;
            break;
        }
    }
    
    //// If all individuals have the same error, the population is just replicated as it is now and then
    //// the values of the member variables of each individual are reset to its initial values and finally
    //// the genome undergoes mutations. Notice that the call to transmit_initial_estimate will take care of
    //// mutation of the initial_estimate in the scenarios where this member variable should evolve.
 
    if(all_equal)
    {
        for(int i = 0; i < Par.POPULATION_SIZE; ++i)
        {
            next_generation.push_back(Population[i] );
            next_generation[i].reset_values(transmit_initial_estimate(Population[i].get_initial_estimate(), Par) );
            next_generation[i].mutate_all_genome(Par);
        }
        Population.swap(next_generation);
    }
    
    //// If at least one individual has different fitness, then chances of leaving offspring are not uniform
    //// throughout the population and reproduction chance should be in some way proportional to fitness...
    else
    {
        //// ...therefore the first thing to do is calculate fitness for each individual. This can be done in
        //// different ways, we chose the following: The individual with the largest lifetime cumulative error
        //// is considered to have a fitness value of 0. For the rest of the individuals, their fitness is the
        //// difference between their own error and the largest error...
        std::vector<double> fitness;
        double largest_error = 0.0;
        for(int i = 0; i < Par.POPULATION_SIZE; ++i) //// Obtain largest error.
        {
            fitness.push_back(Population[i].get_error());
            if(fitness[i] > largest_error)
                largest_error = fitness[i];
        }
        
        for(auto & ind: fitness) ind = largest_error - ind; //// Obtain fitness for each individual by subtracting
                                                            //// its lifetime error from the largest error.
        
        //// ...then a probability distribution is created where each individual probability is weighted by its
        //// fitness...
        std::discrete_distribution<int> weightedReproduction(fitness.begin(), fitness.end());
        
        //// ...and finally individuals that reproduce are drawn randomly from this "weighted lottery" and then
        //// the values of the member variables of each individual are reset to its initial values and the
        //// genome undergoes mutations. Notice that the call to transmit_initial_estimate will take care of
        //// mutation of the initial_estimate in the scenarios where this member variable should evolve.
        for(int i = 0; i < Par.POPULATION_SIZE; ++i)
        {
            int Reproductor = weightedReproduction(rngMutations);
            next_generation.push_back(Population[Reproductor] );
            next_generation[i].reset_values(transmit_initial_estimate(Population[Reproductor].get_initial_estimate(), Par) );
            next_generation[i].mutate_all_genome(Par);
        }
        Population.swap(next_generation);
    }
}


//// This function generates data to be used to create "Rescorla Wagner Plots" which depict the updating behaviour
//// of networks and rules by plotting change in the reward probability estimate as a function of the difference
//// between reward value and old estimate. In these plots, the Rescorla Wagner rule produces a straight line with
//// slope equal to its learning rate. If a learning agent of Network type produced a similar line, then that would
//// mean that experimentally observable behaviour would be the same for both the Rescorla Wagner rule and the
//// Network. Thus, this provides an intuitive way to compare how different agents work compared to each other.
template<class N, class P, class T>
void rescorla_wagner_plot(P & Population, const Parameters & Par, const int &generation, const int &samples, std::ofstream &RW_Plot)
{
    //// Create a bernoulli distribution to obtain rewar or non-reward where both outcomes are equally likely.
    std::bernoulli_distribution Reward_or_not(0.5);
    
    //// Instantiate object to calculate third input.
    T third_input_type;
    
    
    for (int i = 0; i < Population.size(); ++i) //// Loop through individuals in the population.
    {
        for (int j = 0; j < samples; ++j) //// Loop to generate a number of data points per individual.
        {
            //// The difference between the two following variables will generate a single x coordinate within the
            //// range -1 and 1.

            //// 1 or 0 obtained with equal probability.
            double reward = static_cast<double>(Reward_or_not(rngPlots));
            //// An arbitrary value for previous_estimate is taken randomly from a uniform distribution
            //// between 0 and 1.
            double previous_estimate = Uniform0to1(rngPlots);
            
            
            //// To take into account the effect of having experienced n learning chances on the updating of the
            //// probability of reward estimate. This value is relevant only for some networks and rules, and in
            //// others is simply ignored. The range for possible values of n differs between different states of
            //// the Parameter UPDATES_NUMBER_TYPE.
            int n;
            if (Par.UPDATES_NUMBER_TYPE == "random")
            {
                std::uniform_int_distribution<> single_n(1, Par.MAX_UPDATES_NUMBER);
                n = single_n(rngPlots);
            }
            else
            {
                std::uniform_int_distribution<> single_n(1, Par.UPDATES_NUMBER);
                n = single_n(rngPlots);
            }
            
            //// Create and initialize input.
            Input input({ 0.0,0.0,0.0 });
            
            //// The input object that will be fed to the learning agent consists of the reward, the
            //// previous_estimate and the n just generated. The n might be further transformed depending on the
            //// class of the third_input_type object, and, as already mentioned might be used or ignored depending
            //// on the type of the learning agent.
            input.reassign_values
            (
             reward,
             previous_estimate,
             third_input_type.get_third_input(n,Population[i].get_Delta_estimate() ) //// The last element in input is safely ignored by networks and rules that don't take updates into account.
             );
            
            //// The learning agent is fed with the input to make use of the learning mechanism and updates the
            //// probability of reward estimate according to it.
            Population[i].update_estimate2(input,Par);
            
            //// Data exported to generate the Rescorla Wagner Plots. For the x coordinate the difference between
            //// the reward and the previous_estimate, for the y coordinate the difference between the
            //// previous_estimate and the new estimate, updated after using the learning mechanism. in the previous
            //// line.
            RW_Plot << generation << "," << i << "," <<
            reward - previous_estimate << "," <<
            Population[i].get_estimate() - previous_estimate << '\n';
            
        } //// End of loop for number of data points.
    } //// End of loop through individuals in the population.
} //// End of rescorla_wagner_plot function.



//// Classes that implement four different ways of determning the value that should be fed as third input
//// when applicable. These are used in the function do_stuff as template arguments for the call to function
//// void main_loop.
class simple_N
{
public:
    double get_third_input(const int n, const double D) { return static_cast<double>(n);}
};
class one_over_N
{
public:
    double get_third_input(const int n, const double D) {return 1.f / static_cast<double>(n);}
};
class Delta_estimate
{
public:
    double get_third_input(const int n, const double D) {return D;}
};
class one_over_N_plus_two
{
public:
    double get_third_input(const int n, const double D) {return 1.f / ( static_cast<double>(n) + 2.f);}
};



//// This function contains the core of the simulation, where the individuals "visit" different environments, learn
//// from their rewarded or non rewarded experiences, and reproduce proportionally to their fitness, for a number of
//// generations.
template<class N, class P, class T>
void main_loop(P & Population, const Parameters & Par)
{
    //// Print to the console a bar that will serve as reference for the progress bar.
    std::cout << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| - Progress Reference"
    << std::endl << std::endl << std::endl;
    
    //// Output file to export the current state of the whole population for each generation.
    std::ofstream Data(get_folder() + Par.SIMULATION_NAME + "__main_Data.txt");

    //// Output file to export the main error of the population for each generation.
    std::ofstream Error(get_folder() + Par.SIMULATION_NAME + "__error_Data.txt");
    
    //// Output file to export Rescorla Wagner Plots for each generation.
    std::ofstream RW_Plot(get_folder() + Par.SIMULATION_NAME + "__RW_Plot_Long.txt");
    
    
    //// Export headers for the output files.
    RW_Plot << "Generation,Individual,Difference,Delta" << '\n';
    Data << "Generation,Individual,Error,Initial_estimate"; //// For Data it is done in two steps because the first
    Population[0].print_labels(Data,Par);          //// elements are standard regardless of the type of learning
    Data << '\n';                                  //// agent while the others correspond to the genome labels
    Error << "Generation,Error" << '\n';               //// which are specific to each type of agent and are generated
                                                   //// by the print_labels function.
    
    
    
    
    //// Create and initialize input.
    Input input({0.0,0.0,0.0});
    
    //// Instantiate object to calculate third input.
    T third_input_type;
    
    
    
    
       for(int generation = 0; generation <= Par.GENERATIONS; ++generation) //// Loop through generations.
    //for(int generation = 0; generation <= 300; ++generation) //TESTING
    {
        for(int test = 1; test <= Par.TESTS_NUMBER; ++test) //// Loop through environments.
        {
            //// The true probability of reward of the current envioronment is drawn from a uniform distribution.
            double Probability = Uniform0to1(rngEnvironment);
            
            //// Create a bernoulli distribution to obtain rewar or non-reward according to current environment
            //// true probability of reward.
            std::bernoulli_distribution Reward_or_not(Probability);
            
            for(int pop = 0; pop < Par.POPULATION_SIZE; ++pop) //// Loop through individuals in the population.
            {
                //// See description in function calculate_updates_number's definition.
                int updates = calculate_updates_number(Par);
                
                //// See description in function select_first_estimate's definition.
                Population[pop].select_first_estimate(Par);
                
                
                for(int up = 1; up <= updates; ++up) //// Loop through learning experiences.
                {
                    //// The input object that will be fed to the learning agent consists of: 1) reward or non
                    //// reward (1 and 0 respectively) obtained from the afore mentioned distribution, 2) the last
                    //// estimate of probability of reward by the current agent, and 3) a third input that is
                    //// calculated in different ways depending on the class of the third_input_type object, and
                    //// which may be used or ignored after all, depending on the value of the Parameter
                    //// UPDATES_INPUT and the type of learning agent.
                    input.reassign_values
                    (
                     static_cast<double>(Reward_or_not(rngEnvironment)),
                     Population[pop].get_estimate(),
                     third_input_type.get_third_input(up,Population[pop].get_Delta_estimate() ) //// The last element in input is safely ignored by networks and rules that don't take updates into account.
                     );
                    
                    //// Learning takes place here. See further description in function update_estimate2's
                    //// definition.
                    Population[pop].update_estimate2(input,Par);
                    
                    
                } //// End of loop through learning experiences.
                
                //// After all learning experiences in the current enviornment have taken place, the cumulative
                //// life time error of the individual is updated. See function update_error's definition for
                //// more details.
                double Difference = Population[pop].get_estimate() - Probability;
                Population[pop].update_error(Difference,Par);
                
            } //// End of loop through individuals in the population.
        } //// End of loop through environments.
        

        //// Data of the current state of the population and its behaviour (rescorla wagner plot) is exported.
        //// Special conditions are set so that not every generation these data is exported to prevent having
        //// very large output files.
        if(generation % 2000 == 0 || generation < 100 || generation > Par.GENERATIONS - 50 )
        {
            for(int ind = 0; ind < Population.size(); ++ind)
            {
                Data << generation << "," << ind << "," << Population[ind].get_error() << "," << Population[ind].get_initial_estimate();
                Population[ind].print_values(Data);
                Data << '\n';
            }
            rescorla_wagner_plot<N, P, T>(Population, Par, generation, 1, RW_Plot);
        }
  
        //// Mean error of the population is estimated and exported.
        double mean_error = 0;
        for(int ind = 0; ind < Population.size(); ++ind)
        {
            mean_error += Population[ind].get_error();
        }
        mean_error /= static_cast<double>(Population.size());
        
        Error << generation << "," << mean_error << '\n';
        
        
        
        //// Dots are printed to the console to show progress bar of the simulation.
        //if(generation % 3 == 0) //TESTING
            if(generation % ( Par.GENERATIONS/(100) ) == 0)
            std::cout << "." << std::flush;
        
        //// Reproduction of the population takes place to go to next generation. See function
        //// reproduce's definition for full description.
        reproduce<N,P>(Population,Par);
        
    } //// End of loop through generations.
} //// End of main_loop function.









#endif /* Simulation_h */
