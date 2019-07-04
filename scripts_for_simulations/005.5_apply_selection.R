#### This is the script that will decide whether a substitution is accepted or not by applying the fitness function ####

library('optparse')

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
  make_option(c('-l', '--plateau_length'), type = 'numeric', default = NULL,
              help = 'The length of the plateau at the top of the gamma distribution (if selected) on each side of the central value', metavar = 'numeric'),
  make_option(c('-s', '--threshold'), type = 'numeric', default = NULL,
              help = 'The factor (between 0 and 1) that will be multiplied times the starting values to obtain the threshold values', metavar = 'numeric')
              
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
threshold_factor <- as.numeric(opt$threshold)

if(fitness_function == 0){
  beta <- as.numeric(opt$beta)
}else if(fitness_function == 1){
  beta <- as.numeric(opt$beta)
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
    # This is the case for stabilizing selection
    y = rep(0, length(deltaG))
    
    for(i in 1:length(deltaG)){
      if((deltaG[i] - deltaGThreshold[i]) <= 0){
        # The case in which this parameter is more stable than the threshold
        new_threshold <- deltaGThreshold[i] - plateau_radius
        y[i] = -log(exp(-(beta*(deltaG[i] - new_threshold)/2))+1)
      }else if((deltaG[i] - deltaGThreshold[i]) > 0){
        # The case in which this parameter is less stable than the threshold
        new_threshold <- deltaGThreshold[i] + plateau_radius
        y[i] = -log(exp(beta*(deltaG[i] - new_threshold))+1)
      }
    }
    out_y = sapply(y, function(x) min(x, -0.1))
    return(out_y)
    
  }
}

# fitness_program will indicate which variables will be considered for selection
# 0 for selection based only on binding energy (Kachroo's low stability)
# 1 for selection based only on protein stability (Kachroo's non-bound)
# 2 for selection based on both protein stability and binding energy (Kachroo's wildtype)
if(fitness_program == 0){
  binding_energy_threshold <- reference_deltaG$Binding_energy*(threshold_factor)
  threshold_deltaG <- c(0, 0, binding_energy_threshold)
}else if(fitness_program == 1){
  stability_threshold_A <- reference_deltaG$ChainA_stability*(threshold_factor)
  stability_threshold_B <- reference_deltaG$ChainB_stability*(threshold_factor)
  # Set the threshold to a very high value to simulate Kachroo's 
  # non-bound scenario
  threshold_deltaG <- c(stability_threshold_A, stability_threshold_B, 1e10)
}else{
  stability_threshold_A <- reference_deltaG$ChainA_stability*(threshold_factor)
  stability_threshold_B <- reference_deltaG$ChainB_stability*(threshold_factor)
  binding_energy_threshold <- reference_deltaG$Binding_energy*(threshold_factor) 
  threshold_deltaG <- c(stability_threshold_A, stability_threshold_B, binding_energy_threshold)
}

fitness_terms_candidate <- fitness_term(candidate_deltaG, threshold_deltaG, fitness_function)
fitness_terms_previous <- fitness_term(previous_deltaG, threshold_deltaG, fitness_function)


fitness_candidate <- sum(fitness_terms_candidate)
fitness_previous <- sum(fitness_terms_previous)

# Write the fitness contributions for each of those terms
sink(paste(outfile_prefix, '_verdict.txt', sep=''))
cat('Step', 'Stability_A', 'Stability_B', 'Binding_energy', 'Total_fitness\n')
cat('Previous_step', fitness_terms_previous, fitness_previous, '\n')
cat('Candidate', fitness_terms_candidate, fitness_candidate, '\n')

# If the fitness is higher, the substitution is accepted right away
if(fitness_candidate > fitness_previous){
  cat('ACCEPTED\n')
  }else{
  # Otherwise, calculate a probability of acceptance between 0 and 1
  Paccept <- exp(-2*N*(fitness_previous - fitness_candidate)) 

  # Draw a random number between 0 and 1 to compare against the 
  # probability of acceptance
  val <- runif(1, min = 0, max = 1)
  if(is.nan(Paccept)){
    cat('REJECTED')
  }else{
    if(val <= Paccept){
      cat('ACCEPTED\n')
    }else{
      cat('REJECTED\n')
    }
  }
}

# Stop printing to the outfile
sink()


