#### This is the script that will decide whether a substitution is accepted or not by applying the fitness function ####

#### The difference between this script and the first version is that this one can implement
#### stabilizing selection by means of gamma distributions.

library('optparse')

# It receives the following arguments:

#### Argument list ####

# $1 = The path to the file that contains the deltaGs for the candidate substitution
# $2 = The path to the file that contains the deltaGs for the previous substitution
# $3 = The path to the output file that will contain the information for this run
# $4 = The fitness program to use
# $5 = The path to the reference deltaG values for the starting complex
# $6 = The value of N for the population size
# $7 = The kind of fitness function (0 for exponential, 1 for a gamma distribution)
# $8 = first value to use (beta if the function will be exponential, shape (alpha) if gamma)
# $9 = when using stabilizing selection, the value for scale (theta)

#### Assigning arguments with commandArgs ####

#args <- commandArgs(trailingOnly = TRUE)
#candidate_file <- args[1]
#previous_file <- args[2]
#outfile <- args[3]
#fitness_program <- args[4]
#reference_file <- args[5]
#N <- args[6]
#fitness_function <- args[7]

# Based on my test subject: substitution 28 for replicate 2 for my homodimers run for 1a82
# These test values give me an acceptance probability of around 11.17%
# test_input <- c(8.15936, 4.90969, -20.4032)
# test_previous <- c(8.15935, 4.90968, -20.4033)
# test_thresholds <- c(9, 4.5, -21)


#### Assigning arguments with optparse ####

option_list <- list(
  make_option(c('-c', '--cand'), type = 'character', default = NULL,
              help = 'File with values for the candidate substitution', metavar = 'character'),
  make_option(c('-p', '--prev'), type = 'character', default = NULL,
              help = 'File with values for the previous structure', metavar = 'character'),
  make_option(c('-o', '--output'), type = 'character', default = NULL,
              help = 'Output file for the results', metavar = 'character'),
  make_option(c('-m', '--prog'), type = 'numeric', default = NULL,
              help = 'The fitness program to use: 0 (low stab), 1 (non-bound), 2 (wildtype)', metavar = 'numeric'),
  make_option(c('-r', '--ref'), type = 'character', default = NULL,
              help = 'File with the reference deltaG values'),
  make_option(c('-N', '--pop_size'), type = 'numeric', default = NULL,
              help = 'The population size number for the Paccept function', metavar = 'numeric'),
  make_option(c('-f', '--fit_func'), type = 'numeric', default = NULL,
              help = 'The fitness function to use: 0 (exponential), 1 (gamma distributed - stabilizing selection)', metavar = 'numeric'),
  make_option(c('-b', '--beta'), type = 'numeric', default = NULL,
              help = 'The beta value for the shape of the exponential curve (if the exponential function is selected)', metavar = 'numeric'),
  make_option(c('-k', '--shape'), type = 'numeric', default = NULL,
              help = 'The k (shape) parameter for the gamma distribution (if selected)', metavar = 'numeric'),
  make_option(c('-t', '--scale'), type = 'numeric', default = NULL,
              help = 'The theta (scale) parameter for the gamma distribution (if selected)', metavar = 'numeric'),
  make_option(c('-l', '--plateau_length'), type = 'numeric', default = NULL,
              help = 'The length of the plateau at the top of the gamma distribution (if selected) on each side of the central value', metavar = 'numeric')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

candidate_file <- opt$cand
previous_file <- opt$prev
outfile_prefix <- opt$output
fitness_program <- as.numeric(opt$prog)
reference_file <- opt$ref
N <- as.numeric(opt$pop_size)
fitness_function <- as.numeric(opt$fit_func)

if(fitness_function == 0){
  beta <- as.numeric(opt$beta)
}else{
  shape <- as.numeric(opt$shape)
  scale <- as.numeric(opt$scale)
  plateau_radius <- as.numeric(opt$plateau_length)
}

candidate_deltaG <- read.table(candidate_file, h = T)
previous_deltaG <- read.table(previous_file, h = T)
reference_deltaG <- read.table(reference_file, h = T)

# A function to calculate the fitness component for each term
fitness_term <- function(deltaG, deltaGThreshold, fitness_function){
  deltaG <- as.numeric(deltaG)
  deltaGThreshold <- as.numeric(deltaGThreshold)
  # The case for the original exponential function
  if(fitness_function == 0){
    return(-log(exp(beta*(deltaG - deltaGThreshold))+1))
  }else if(fitness_function == 1){
    # The case for stabilizing selection
    # A gamma probability distribution that is centered around the wildtype value and has our parameters,
    # as well as a faster fitness descent on the side of chain and complex destabilization than on the
    # side of overstabilization. The maximum value will be set to 1.
    
    # I must convert my input value so that the wildtype becomes this maximum value
    # I can take the maximum value to any value if I subtract that value and the max_val
    # to all the points in the distribution. I can just shift my values accordingly.
    max_val <- (shape - 1)*scale
    shift <- max_val - deltaGThreshold
    
    
    # Save some variables for easier reading
    # To have the max_val be 1
    norm_factor <- dgamma(max_val, shape = shape, scale = scale)
    # To move the observed deltaG values to the corresponding region under the distribution
    shifted_deltaG <- deltaG + shift 
    # To invert the distribution. For example, a Gamma distribution with parameters
    # shape = 2 and scale = 2 decays more slowly to the right (destabilization), which
    # is not what I want. Using this inversion factor (basically 2*<distance to the threshold>),
    # I will get a mirrored distribution.
    inversion_factor <- 2*(deltaG - deltaGThreshold)
    
    # A shift for the plateau, pretty much the plateau radius with a sign that depends on the position of
    # deltaG with respect to the deltaGThreshold
    plateau_shift = (deltaG - deltaGThreshold)*plateau_radius/abs((deltaG - deltaGThreshold))
    
    # Use everything to get the fitness value
    # return(dgamma(shifted_deltaG - inversion_factor, shape = shape, scale = scale)/norm_factor)
    
    # Using a plateau
    return(ifelse(abs(deltaG - deltaGThreshold) < plateau_radius, # Test if we are in the plateau region
      1, # Take a 1 if we are inside
      dgamma(shifted_deltaG - inversion_factor + plateau_shift, shape = shape, scale = scale)/norm_factor # Sample from the distribution if we are outside
      )
    )
  }
}

# I will use the reference values and the selected program to decide the values for the threshold
# fitness_program will indicate how the thresholds for fitness will be used on homodimers
# 0 for selection based only on binding energy (Kachroo's low stability)
# 1 for selection based only on protein stability (Kachroo's non-bound)
# 2 for selection based on both protein stability and binding energy (Kachroo's wildtype)
if(fitness_program == 0){
  binding_energy_threshold <- reference_deltaG$Binding_energy*(1 - 1/N)
  threshold_deltaG <- c(0, 0, binding_energy_threshold)
}else if(fitness_program == 1){
  stability_threshold_A <- reference_deltaG$ChainA_stability*(1 - 1/N)
  stability_threshold_B <- reference_deltaG$ChainB_stability*(1 - 1/N)
  # Kachroo set the deltaG for binding energy to infinity. I will use 1e10 as a proxy for that.
  threshold_deltaG <- c(stability_threshold_A, stability_threshold_B, 1e10)
}else{
  stability_threshold_A <- reference_deltaG$ChainA_stability*(1 - 1/N)
  stability_threshold_B <- reference_deltaG$ChainB_stability*(1 - 1/N)
  binding_energy_threshold <- reference_deltaG$Binding_energy*(1 - 1/N) 
  threshold_deltaG <- c(stability_threshold_A, stability_threshold_B, binding_energy_threshold)
}

fitness_terms_candidate <- fitness_term(candidate_deltaG, threshold_deltaG, fitness_function)
fitness_terms_previous <- fitness_term(previous_deltaG, threshold_deltaG, fitness_function)


fitness_candidate <- sum(fitness_terms_candidate)
fitness_previous <- sum(fitness_terms_previous)

# I will write the fitness contributions for each of those terms
sink(paste(outfile_prefix, '_verdict.txt', sep=''))
cat('Step', 'Stability_A', 'Stability_B', 'Binding_energy', 'Total_fitness\n')
cat('Previous_step', fitness_terms_previous, fitness_previous, '\n')
cat('Candidate', fitness_terms_candidate, fitness_candidate, '\n')

# If the fitness is higher, the substitution is accepted right away
if(fitness_candidate > fitness_previous){
  cat('ACCEPTED\n')
  }else{
  # Otherwise, I need to calculate a probability of acceptance between 0 and 1
  Paccept <- exp(-2*N*(fitness_previous - fitness_candidate)) 

  # I will draw a random value from an uniform distribution to decide whether it should be accepted or not
  # If Paccept = 0.70, I want to accept the mutation 70% of the time. As 70% of the values between 0 and 1 are
  # less or equal to 0.70, I will only accept it if the drawn value is less or equal to 0.70
  # I tried this in a loop for 1000 draws several times, and I get close to the expected number of accepted substitutions
  val <- runif(1, min = 0, max = 1)
  if(val <= Paccept){
    cat('ACCEPTED\n')
  }else{
    cat('REJECTED\n')
  }
}

# Stop printing to the outfile
sink()


