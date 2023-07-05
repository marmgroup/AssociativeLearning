//
//  Networks.h

//
//  Created by Emiliano Méndez Salinas.
//  Copyright © 2020 Emiliano Méndez Salinas. All rights reserved.
//

#ifndef Networks_h
#define Networks_h



//// Defined in main.cpp
template<class N, class P>
void do_stuff(int & argc, P & Pop, const Parameters & Par);


extern std::mt19937_64 rngMutations;










//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//// ++++++++++ CLASS DEFINITIONS DIRECTLY RELATED TO ANNs +++++++++++++++++++++

//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//// Input class objects are used as the input for a ANNs (or learning rules). Although all Input objects have
//// the third (number) Value specified, whether this is used or not during the operation of ANNs depends on
//// the parameter (UPDATES_INPUT) as determined by the activate(const Input, const Parameters) member function
//// of class Network.

class Input
{
    template <class N> friend class Node;
private:
    std::vector<double> Values;
public:
    Input(const std::vector<double> in_values) : Values(in_values) {}

    //// Get respective values
    double reward() const {return Values.at(0);}              //// Use of at() to prevent out of bound access.
    double estimate() const {return Values.at(1);}
    double number() const {return Values.at(2);}
    double get_val(const int &i) const {return Values[i];}

    //// Set respective values
    void reward(const double &val) {Values[0] = val;}
    void estimate(const double &val) {Values[1] = val;}
    void number(const double &val) {Values[2] = val;}
    void reassign_values(const double & reward, const double &estimate, const double &number)
    {
        Values[0] = reward;
        Values[1] = estimate;
        Values[2] = number;
    }
};




//// Transfer function classes. These are the classes used to instantiate templated Networks. 
//// Each class has its respective definition of how the data that comes out of a node will be transformed. 
//// Hence, the two template classes used to instantiate a Network determine, respectively, which kind of transfer functions will the networks use in a) all processing nodes, and b) output nodes.

class none
{
    template<class P, class O> friend class Network;
    double transfer(const double & weighted_sum) { return weighted_sum; }
    std::string print() {return "none";}
};

class tanh_tf
{
    template<class P, class O> friend class Network;
    double transfer(const double & weighted_sum) { return tanh(weighted_sum); }
    std::string print() {return "tanh";}
};

class log_tf
{
    template<class P, class O> friend class Network;
    double transfer(const double & weighted_sum) { return 1.0 / (1.0 + exp(-(weighted_sum) ) ); }
    std::string print() {return "log";}
};

class relu_tf
{
    template<class P, class O> friend class Network;
    double transfer(const double & weighted_sum) { return weighted_sum > 0 ? weighted_sum : 0; }
    std::string print() {return "relu";}
};







//// Template class. The classes used to instantiate Network objects must be any of the transfer function classes defined beforehand. 
//// The first class determines the type and transfer function for all processing nodes and the second class does the same for all the output nodes.

template<class P, class O>
class Network
{
    template <class Processing_layer, class Output_node>
    friend void create_population(int & argc, char* argv[], const Parameters &Par);
private:
    
    //// Static variables, equal for, and shared among, all the Network objects within a simulation.
    
    static size_t weights_number; //// Size of the genome,: number of weights + biases the ANN contains
    static size_t temp_values_size; //// Size of the vectors old_values and new_values
    static std::vector<double> old_values; //// These two containers temporarilly store the weighted sums of each
    static std::vector<double> new_values; //// layer and are overwritten with the output of each subsequent layer.
    
    
    
    //// Genome containing the values of all the weights and biases that the network consists of.
    //// Since this is linear (a simple vector of doubles), the identity of each value (what weight/bias
    //// of the network it is) can only be deduced from the following Network information contained in the
    //// Parameter file: a) number and size of each layer (NODES_PER_LAYER), b) presence/lack of bias in
    //// processing (PROCESSING_BIAS) and output nodes (OUTPUT_BIAS), and c) presence/lack of the third value of the
    //// input (UPDATES_INPUT). Therefore, member function activate(const Input, const Parameters) makes use of
    //// all this parameters information in order to correctly carry on the operation of the network. And similarly,
    //// member function print_labels(std::ofstream, const Parameters) makes use of the same parameters information
    //// to print the headers of the Data ofstream object ("...main_data.txt") so that the human user can understand
    //// which weight/bias each output value corresponds to.
    std::vector<double> Connection_weights;
    
    
    //// Instantiation of transfer function objects based on the classes passed to the Network template. These
    //// determine which transfer functins will the network use in processing and output nodes, respectively.
    P proc_tf;
    O out_tf;
    
    //// Returns values that will fill member Connection_weights. More descritpion below in the function definition.
    std::vector<double> create_weights(const Parameters & Par);
    
    
public:
    //// Constructor. Only gives values to member Connection_weights, for which function
    //// create_weights(Par) is called.
    Network (const Parameters & Par) : Connection_weights(create_weights(Par) )
    {}
    
    //// Single call to mutate the genome of a Network during/after reproduction. It itself uses other functions
    //// that are in charge of actually changing the individual values within the genome.
    void mutate_all_genome(const Parameters &Par)
    {
        for(int i = 0; i < Connection_weights.size(); ++i) {
            Connection_weights[i] = value_after_mutation(Par,Connection_weights[i]); }
    }

    //// Makes the network "react" to an input and return an output. Full description in the declaration below.
    double activate(const Input & in_vec, const Parameters & Par);
 
    //// Makes use of the parameters information to print the headers of the Data ofstream object
    //// ("...main_data.txt") so that the human user can understand which weight/bias each output value
    //// corresponds to.
    void print_labels(std::ofstream &ofs, const Parameters & Par) const;

    //// The print_labels(std::ofstream, const Parameters) function does the difficult job of matching genome values
    //// to their actual role/position within the network. This one only prints the values in the order they are
    //// found in the genome Connection_weights.
    void print_values(std::ofstream &ofs) const
    {
        for(int i = 0; i < Connection_weights.size(); ++i) {
            ofs << "," << Connection_weights[i]; }
    }
    
};



//// Static members declaration outside of the class Network, otherwise the compiler complains. The initialization
//// of these variables is done in function void create_population in Pop.h
template<class P, class O>
size_t Network<P,O>::weights_number;
template<class P, class O>
size_t Network<P,O>::temp_values_size;
template<class P, class O>
std::vector<double> Network<P,O>::old_values;
template<class P, class O>
std::vector<double> Network<P,O>::new_values;









//// Returns the whole vector of weights/biases used by the constructor of Network. Here the identity of each value
//// (which weight/bias within the network it is) is not important, but the number of values has to be correct.
//// Therefore the function relies on the static member "weights_number". The only other action it accomplishes
//// is choosing between initializing the weights of the network as zeros or random values (within certain bounds)
//// all based on the Parameters of the simulation.

template<class P, class O>
std::vector<double> Network<P,O>::create_weights(const Parameters & Par)
{
    std::vector<double> values;
    if(Par.INITIAL_WEIGHTS == "random")
    {
        std::uniform_real_distribution<double> Uniform(Par.LOW_BOUND, Par.UPP_BOUND);
        
        for(int i = 0; i < weights_number; ++i)
            values.push_back( Uniform(rngMutations) );
        return values;
    }
    else if(Par.INITIAL_WEIGHTS == "0")
    {
        for(int i = 0; i < weights_number; ++i)
            values.push_back(0);
        return values;
    }
    else
    {
        std::cout << "ERROR: UNKNOWN PARAMETER InitialWeights" << std::endl;
        exit(EXIT_FAILURE);
    }
}








//// This is the "blackboxed" and "arbitrary" function that makes the Network work properly. Apart from actually
//// computing all the calculations/operations performed by the network, by determining the order of operations
//// based on parameter values, in effect it interprets the identity of each value in the genome Connection_weights,
//// so that operands in each calculation are the correct ones. To temporarily store values that are transferred by
//// each layer, it makes use of static vectors "old_values" and "new_values".
template<class P, class O>
double Network<P,O>::activate(const Input & in_vec, const Parameters & Par)
{
    int i = 0;
    //// old_values is always used to give an input to the next layer. Input to the first layer is exactly the
    //// Input to the Network, so here network Input is copied to old_values.
    for(; i < (2 + Par.UPDATES_INPUT); ++i) //// If parameter UPDATES_INPUT is 0, the third value of Network input
        old_values[i] = in_vec.get_val(i);  //// is not copied and thus ignored. (See comments on class Input)
    
    //// Contents of old_values will be used as one set of factors to calculate the weighted sum for each node (the
    //// other set of factor being the weights). Thus the last value of old_values vector corresponds to 1 if
    //// processing nodes have a bias. Otherwise, although it takes value 0, this is irrelevant because it will be
    //// ignored anyway as explained in comment below for "size_t weights_n". 

    old_values[i] = Par.PROCESSING_BIAS;
    
    size_t pos_w = 0; //// Counter. Will serve as index to weights in genome so that they are used sequentially
    size_t pos_n = 0; //// Counter. Will be reset to 0 to serve as index for the weighted sum of each node.
    size_t weights_n = 2 + Par.UPDATES_INPUT + Par.PROCESSING_BIAS; //// weights_n controls how many products
    //// will be used for the weighted sum of each node (since all nodes in the same layer receive same input,
    //// weights_n needs to be recalculated only every new layer). It's initial value must be between 2 and 4
    //// depending on whether the network has processing bias and whether there should be third input to the ANN.
    //// If weights_n is not correct, pos_w will get out of phase with respect to pos_n and the Network will not
    //// compute the calculations it is supposed to, giving a wrong output in the end. Depending on other features
    //// of the Network, this might also cause an "out of bounds" error, but this kind of error will not always
    //// be the case!
    
    for(int layer = 0; layer < Par.NODES_PER_LAYER.size() - 1; ++layer) //// If the Network has only output layer
    {                                                                   //// this loop will not execute.
        for(int node = 0; node < Par.NODES_PER_LAYER[layer]; ++node)
        {
            for(int values = 0; values < weights_n; ++values) //// See comment above
            {
                new_values[node] += Connection_weights[pos_w] * old_values[pos_n]; //// Compute the weighted sum
                ++pos_w; //// move on to next weight
                ++pos_n; //// move on to next node
            }
            pos_n = 0; //// pos_n must be reset so that the same input values are used again and
        }              //// again by every node of the same layer.
        for(int i = 0; i < new_values.size(); ++i)               //// The weighted sum of every node is transformed
        {                                                        //// using the indicated transfer function and
            old_values[i] = proc_tf.transfer(new_values[i]);     //// the resultant values (one per node) are
            new_values[i] = 0;                                   //// copied to the container that will be used as
        }                                                        //// input for the next layer.
        
        //// weights_n must be recalculated since the number of products for the weighted sum of next layer depends
        //// on both, the number of nodes in the current layer, and whether processing bias is present or not.
        weights_n = Par.NODES_PER_LAYER[layer] + Par.PROCESSING_BIAS;
        
        //// Partially repeated comment. The last value of old_values vector corresponds to 1 if processing nodes
        //// have a bias. Otherwise, although it takes value 0, this is irrelevant because it will be ignored
        //// anyway as explained above in comment for "size_t weights_n". 
        old_values[Par.NODES_PER_LAYER[layer]] = Par.PROCESSING_BIAS;
    }
    
    
    size_t ind; //// Will store the size of the input that is fed to the output layer.
    
    //// pos_w would be larger than 0 only if there are processing layers, in which case
    //// Par.NODES_PER_LAYER.end()[-2] indicates the size of the last processing layer, which is itself the size
    //// of the input to output layer. Otherwise, there are no processing layers and the input to the network is
    //// the input to output layer. Additionally, in the next line, it must be considered if output nodes have a
    //// bias, which will be used in the weighted sum. Again... The last value of old_values vector corresponds to 1
    //// if processing nodes have a bias. Otherwise, although it takes value 0, this is irrelevant because it will
    //// be ignored anyway as explained above in comment for "size_t weights_n".

    pos_w > 0 ? ind = Par.NODES_PER_LAYER.end()[-2] : ind = 2 + Par.UPDATES_INPUT;
    old_values[ind] = Par.OUTPUT_BIAS;
    
    
    weights_n -= Par.PROCESSING_BIAS; /// since this will set the size of the input to the output layer,
    weights_n += Par.OUTPUT_BIAS;     /// processing bias should be neglected and output bias considered.
    
    //// Loop for the output layer. If there were no processing layers, execution skips previous loop and enters
    //// directly into this one.
    for(int node = 0; node < Par.NODES_PER_LAYER.back(); ++node)
    {
        for(int values = 0; values < weights_n; ++values)
        {
            new_values[node] += Connection_weights[pos_w] * old_values[pos_n]; //// Compute the weighted sum
            ++pos_w; //// move on to next weight
            ++pos_n; //// move on to next node
        }
        pos_n = 0; //// pos_n must be reset so that the same input values are used again and
    }              //// again by every node of the same layer.
    for(int i = 0; i < new_values.size(); ++i)                   //// The weighted sum of every node is transformed
    {                                                            //// using the indicated transfer function and
        old_values[i] = out_tf.transfer(new_values[i]);          //// the resultant values (one per node) are
        new_values[i] = 0;                                       //// copied to the container that will be returned
    }                                                            //// by the function in the next code line.
    
    //// IMPORTANT: This code line (and the type of the function return) will eventually have to be changed if more
    //// than one output node are needed.
    return old_values[0];
    
} //// End of activate function.

//// Makes use of the parameters information to print the headers of the Data ofstream object
//// ("...main_data.txt") so that the human user can understand which weight/bias each output value
//// corresponds to.
template<class P, class O>
void Network<P,O>::print_labels(std::ofstream &ofs, const Parameters & Par) const
{
    //// This variable will indicate how many incomming weights each and all nodes of each layer have. To begin
    //// with, the number of weights corresponds to the number of inputs to the whole network. This is true
    //// regardless of whether there are processing layers or not.
    size_t previous = 2 + Par.UPDATES_INPUT;
    
    //// If no processing layer exists, this loop will be skipped altogether.
    for(int i = 0; i < Par.NODES_PER_LAYER.size() - 1; ++i) //// Loop through processing layers.
    {
        //// In the following lines, a unique composite identifier string is created for each weight in
        //// processing layers...
        
        //// First the number of the processing layer (P+number)...
        std::string layer_prefix = "P" + std::to_string(i) + "_";
        for(int j = 0; j < Par.NODES_PER_LAYER[i]; ++j) //// Loop through nodes in the current layer.
        {
            //// Then the number of the node and its type of transfer function (N+number+processing_tf)...
            std::string node_prefix = "N" + std::to_string(j) + "_" + Par.PROCESSING_TF + "_";
            for(int w = 0; w < previous; ++w) //// Loop through weigths in the current node.
            {
                //// Finally the number of the weight (W+number).
                std::string weight_label = "W" + std::to_string(w);
                //// Each element of the unique composite identifier is printed to the output file.
                ofs << "," << layer_prefix << node_prefix << weight_label;
            } //// End of loop through weights.
            
            //// If the processing nodes have a bias, one more composite identifier is printed which corresponds to
            //// the bias to current node, therefore it uses "B" instead of "W"+number.
            if(Par.PROCESSING_BIAS)
                ofs << "," << layer_prefix << node_prefix << "B";
        } //// End of loop through nodes.
        
        //// The number of incomming weights into all and each node of the next layer is the number of nodes in the
        //// current layer.
        previous = Par.NODES_PER_LAYER[i];
    } //// End of loop through layers.

    
    //// Again, a unique composite identifier string is created for each weight, but now in the output layer...

    //// Since there is only one output layer, only "O" and no number is needed...
    std::string layer_prefix = "O";
    for(int j = 0; j < Par.NODES_PER_LAYER.back(); ++j) //// Loop through nodes in the output (always last) layer.
    {
        //// Then the number of the node and its type of transfer function (N+number+output_tf)...
        std::string node_prefix = std::to_string(j) + "_" + Par.OUTPUT_TF + "_";
        for(int w = 0; w < previous; ++w) //// Loop through weigths in the current node.
        {
            //// Finally the number of the weight (W+number).
            std::string weight_label = "W" + std::to_string(w);
            //// Each element of the unique composite identifier is printed to the output file.
            ofs << "," << layer_prefix << node_prefix << weight_label;
        }  //// End of loop through weights.
        
        //// If the output nodes have a bias, one more composite identifier is printed which corresponds to the
        //// bias to the current node, therefore it uses "B" instead of "W"+number.
        if(Par.OUTPUT_BIAS)
            ofs << "," << layer_prefix << node_prefix << "B";
    } //// End of loop through nodes.

} //// End of print_labels function.




//// Uses the information of the parameters to calculate and return the total number of weights (including biases)
//// that the genome of the network must contain. Called to initialize the static member variable
//// "static size_t weights_number".
size_t calculate_weights_number(const Parameters & Par)
{
    size_t previous = 2 + Par.UPDATES_INPUT; //// Number of weights into the first (possibly only) layer.
    size_t count = 0;    //// Initialize counter of number of processing nodes to 0.
    size_t weights = 0;     //// Initialize counter of number of weights to 0.
    
    //// Loop through processing layers. Skipped if there are no processing layers.
    for(int i = 0; i < Par.NODES_PER_LAYER.size() - 1; ++i)
    {
        count += Par.NODES_PER_LAYER[i]; //// Add number of nodes in the current layer.
        weights += Par.NODES_PER_LAYER[i] * previous; //// Add number of weights in the current layer (incomming
                                                      //// weights plus number of nodes).
        previous = Par.NODES_PER_LAYER[i]; //// Reset value to number of weights incomming into next layer.
    }
    
    weights += count * Par.PROCESSING_BIAS;    //// If processing layers exist and have a bias, one additional
                                               //// weight is added per processing node.
    weights += Par.NODES_PER_LAYER.back() * previous; //// Add the number of weights incomming into the output layer
                                                      //// (nodes in previous layer, or input, times nodes in output
                                                      //// layer).
    weights += Par.NODES_PER_LAYER.back() * Par.OUTPUT_BIAS; //// If output nodes have a bias, one additional weight
                                                             //// is added per output node.
    return weights;
    
} 



//// The output of this function will be used to determine the size of static member variable vector old_values and
//// new_values. Since these are used within the activate function to temporarily store the weighted sums of each
//// layer (that is, the values that will be input fed to subsequent layer) the output of this function should
//// match the size of the largest among all layers and input vector.
size_t calculate_temp_values_size(const Parameters & Par)
{
    //// Start with the size of the input and replace it if any layer has larger size than that.
    size_t initial = 2 + Par.UPDATES_INPUT;
    for(int i = 0; i < Par.NODES_PER_LAYER.size(); ++i)
        if(Par.NODES_PER_LAYER[i] > initial)
            initial = Par.NODES_PER_LAYER[i];
    ++initial; //// Instruct to add an additionl space to the vector in case there is a bias. See comments to
               //// activate function to get a clearer picture.
    return initial;
}




//// Template function that allows automatic implementation of a probability distribution as specified by the
//// template argument. Then returns a double value drawn from such distribution with given mean and deviation
//// values.
template<class Probability_distribution>
double mutate(double mean, double dispersion)
{
    Probability_distribution mut(mean, dispersion);
    return mut(rngMutations);
}

//// It serves to determine whether a mutation takes place or not. If it does, this works as a wrapper to the
//// previous function and actually outputs the output of the previous function, by calling it with the template
//// argument according to type of mutation distribution specified in Parameters. Notice that in that case, the old
//// weight is not modified itself but used to create a new random value around this weight as mean. That is why the
//// wheight argument can be const here. If no mutation occurs, the weight argument of the function is returned
//// without modification.
double value_after_mutation(const Parameters & Par, const double & weight)
{
    std::bernoulli_distribution Mutation_happens(Par.MUTATION_PROBABILITY);
    if(Mutation_happens(rngMutations) )
    {
        if(Par.MUT_DISTRIBUTION == "normal")
        {
            return mutate<std::normal_distribution<double> >(weight,Par.MUTATION_DISPERSION);
        }
        else
            if(Par.MUT_DISTRIBUTION == "cauchy")
            {
                return mutate<std::cauchy_distribution<double> >(weight,Par.MUTATION_DISPERSION);
            }
            else
                std::cout << "ERROR: UNKNOWN PARAMETER Mutation Distribution" << std::endl;
        exit(EXIT_FAILURE);
    }
    else
        return weight;
}







#endif /* Networks_h */
