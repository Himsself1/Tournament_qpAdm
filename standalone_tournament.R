# * Libraries

list_of_packages <- c( 
 "ggplot2", "devtools",
  "argparse", "stringr",
  "Cairo", "tibble",
 "reshape", "dplyr",
 "yaml"
)

for (i in list_of_packages) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = T)
  }
}

devtools::install_github("uqrmaie1/admixtools")
library(admixtools)


# * Command line arguments

args<-commandArgs(TRUE)
yaml_input <- args[1]

## yaml_input <- "./qpadm_inference/apoikia_tournament_example.yaml"

# ** Verify arguments

if( file.exists( yaml_input ) ){
  input_params <- read_yaml( yaml_input )
}else{
  cat(sprintf( "YAML file: %s not found. \n", yaml_input ))
  stop("Stop 0")
}

# *** Population file names

if( !file.exists( input_params$file_with_ind_names ) ){
  cat(sprintf( "File %s does not exist. \n", input_params$file_with_ind_names ))
  stop( "Stop 1" )
}else{
  ind_names <- read.table(input_params$file_with_ind_names)[,1]
}

if( !file.exists( input_params$file_with_pop_names ) ){
  cat(sprintf( "File %s does not exist. \n", input_params$file_with_pop_names ))
  stop("Stop 2")
}else{
  pop_names <- read.table(input_params$file_with_pop_names)[,1]
}

if( !dir.exists(input_params$path_to_eigenstrat) ){
  cat(sprintf( "Input folder: %s does not exist. \n",
          input_params$path_to_eigenstrat
          ))
  stop("Stop 3")
}

eig_files <- list.files(path = input_params$path_to_eigenstrat,
           pattern = input_params$prefix)
if( length(eig_files) < 3 ){
  cat(sprintf("EIGENSTRAT files missing or incomplete. \n"))
  stop("Stop 4")
}else{
  input_prefix <- paste0(c(input_params$path_to_eigenstrat, input_params$prefix), collapse = '/')
}

# *** Check suitability of targets, left and right populations

if( length(input_params$competing_pops) < 2 ){
  cat(sprintf( "Please select populations to compete \n" ))
  stop("Stop 5")
}else{
  competing_pops <- input_params$competing_pops
  if( sum(!(competing_pops %in% pop_names)) ){
    cat(sprintf( "Populations  \n %s are not in \n %s \n",
                paste0(competing_pops[!(competing_pops %in% pop_names)], collapse = ' '),
                input_params$file_with_pop_names ))
    stop("Stop 6")
  }
}

if( input_params$target_pop_name == '' ){
  cat(sprintf( "Please select a target population.\n" ))
  stop( "Stop 6.5" )
}else{
  target_pop_name <- input_params$target_pop_name
  if( sum(!(target_pop_name %in% pop_names)) ){
    cat(sprintf( "Populations \n %s are not in \n %s \n",
                paste0(target_pop_name[!(target_pop_name %in% pop_names)], collapse = ' '),
                input_params$file_with_pop_names ))
    stop("Stop 7")
  }
}

other_left_pops <- input_params$other_left_pops
if( sum(!(other_left_pops %in% pop_names)) ){
  cat(sprintf( "Population(s) \n %s are not in \n %s \n",
              paste0(target_pop_namee[!(other_left_pops %in% pop_names)], collapse = ' '),
              input_params$file_with_pop_names ))
  stop("Stop 8")
}

if( target_pop_name %in% competing_pops ){
  cat(sprintf( "Remove target population: %s from competing populations. \n",
              target_pop_name))
}

if( target_pop_name %in% other_left_pops ){
  cat(sprintf( "Remove target population: %s from left populations. \n",
              target_pop_name))
}

# *** Create output folder

out_dir_full_name <- paste0( c(input_params$output_folder, input_params$run_name), collapse = '/')
out_dir_for_plots <- paste0(out_dir_full_name, "/plots/", collapse = '')
out_dir_for_stats <- paste0(out_dir_full_name, "/stats/", collapse = '')
dir.create( out_dir_for_plots, recursive = TRUE )
dir.create( out_dir_for_stats, recursive = TRUE )

# * Building Models

all_pops <- unique(pop_names)

## Populations that will never be left or target are going to always be right
always_right <- all_pops[!(all_pops %in% c(target_pop_name, competing_pops, other_left_pops))] 

# Need to make lists of `left`, `right`, `target` and `excluded` populations.
# Each line is going to correspond to a seperate model.
left_all = list()
right_all = list()
exclude = c()

## For each pair (A,B) of competing pops we need 2 models:
## A = left and B in right | A = left and B NOT in right
for(i in 1:length(competing_pops)){
  for( j in 1:length(competing_pops) ){
    if( i == j ) { next }
    left_all <- c(left_all, list(c(competing_pops[i], other_left_pops)))
    left_all <- c(left_all, list(c(competing_pops[i], other_left_pops)))
    right_all <- c(right_all, list(c(competing_pops[-i], always_right)))
    right_all <- c(right_all, list(c(competing_pops[c(-i,-j)], always_right)))
    exclude <- c(exclude, c("None", competing_pops[j]))
  }
}

target_all <- rep(target_pop_name, length(left_all))

qp_models <- tibble(
  left = left_all,
  right = right_all,
  target = target_all
)

# * Runnign qpAdm

f2_blocks_for_all_pops <- f2_from_geno(
  pref = input_prefix,
  pops = pop_names,
  inds = ind_names,
  adjust_pseudohaploid = FALSE  
)

all_qpadms <- qpadm_multi(
  data = f2_blocks_for_all_pops,
  models = qp_models
  )

# ** Data processing

# *** Information for all models

feasibility <- c()
weights <- list()
p_values <- c()
for (i in 1:length(all_qpadms)) {
  feasibility <- c(feasibility, all_qpadms[[i]]$popdrop$feasible[1])
  weights[[i]] <- unname(all_qpadms[[i]]$popdrop[1, c(7, 8)])
  p_values <- c(p_values, all_qpadms[[i]]$popdrop[1, 5])
}
  
## `data_v1` will be the data frame with all the information tidied up
## In the final data frame, `p_value_all` and `feasible_all` show the p_value and
## feasibility of the model that doesn't exclude the respective population.
data_v1 <- data.frame(
  left_1 = do.call(rbind, left_all)[, 1],
  weights_1 = as.data.frame(do.call(rbind, weights))[, 1],
  exclude = exclude,
  p_values = unname(unlist(p_values)),
  p_values_all = unname(unlist(p_values[rep(which(exclude == "None"), each = 2)])),
  feasible = feasibility,
  feasible_all = feasibility[rep(which(exclude == "None"), each = 2)]
)

write.table(
  x = data_v1, quote = FALSE, sep = '\t',
  file = paste0(c(out_dir_for_stats, "all_models.tsv"), collapse = '')
)

# ** Preparing scoring functions

data_for_scoring_function <- as.data.frame(
    data_v1[exclude != "None", -2 ],
    row.names = 1:sum(exclude != "None")
  )

# *** Calculating Pairwise `Resilience`

## Returns 1 if resillient, -1 if not resillient and 0 if a decision cannot be made.
score_stef <- c()
for (i in 1:nrow(data_for_scoring_function)) {
  if ((data_for_scoring_function$p_values[i] < 0.05) & (data_for_scoring_function$p_values_all[i] < 0.05)) {
    ## Both models are bad
    score_stef <- c(score_stef, 0)
    next 
  } else if ((data_for_scoring_function$p_values[i] < 0.05) & (data_for_scoring_function$p_values_all[i] > 0.05)) {
    ## Include is better, so I need to check if it is also feasible.
    if (data_for_scoring_function$feasible_all[i] == TRUE) {
      score_stef <- c(score_stef, 0)
      next
    } else {
      score_stef <- c(score_stef, 0)
      next
    }
  } else if ((data_for_scoring_function$p_values[i] > 0.05) & (data_for_scoring_function$p_values_all[i] < 0.05)) {
    ## Exclude is better, so I need to check if it is also feasible.
    if (data_for_scoring_function$feasible[i] == TRUE) {
      score_stef <- c(score_stef, -1)
      next
    } else {
      score_stef <- c(score_stef, 0)
      next
    }
  } else {
    ## Both P_include and P_exclude are > 0.05. Need to check if they are also feasible.
    if ((data_for_scoring_function$feasible[i] == FALSE) & (data_for_scoring_function$feasible_all[i] == FALSE)) {
      ## Both are not feasible => tie.
      score_stef <- c(score_stef, 0)
      next
    } else if ((data_for_scoring_function$feasible[i] == FALSE) & (data_for_scoring_function$feasible_all[i] == TRUE)) {
      score_stef <- c(score_stef, 0) # Need to discuss this. Include > Exclude
      next
    } else if ((data_for_scoring_function$feasible[i] == TRUE) & (data_for_scoring_function$feasible_all[i] == FALSE)) {
      score_stef <- c(score_stef, -1) # Exclude > Include
      next
    } else {
      ## Both P-values are > 0.05 and both are feasible => Inference is good
      score_stef <- c(score_stef, 1)
      next
    }
  }
}
data_for_scoring_function$score_stef <- score_stef

#### TODO: This needs testing with a real dataset #### DONE (04/03)

# *** Calculating Lazaridis score

## If A is resillient to B and B is not resillient to A then A vs B = 1
## If A is resillient to B and B is resillient to A then A vs B = 0
## If A is not resillient to B and B is resillient to A then A vs B = -1
## The sign of the score denotes which population is favored in the comparison.
data_versus <- data_for_scoring_function[, c(1:2, 7)]
results_versus <- data_frame()
for (l1 in 1:(length(competing_pops)-1)) {
  for (ex in (l1+1):(length(competing_pops))) {
    ## print( paste0(c(l1, ex), collapse = ' '))
    if (l1 == ex) {
      next
    }
    index_model_A <- which(
    (data_versus$left_1 == competing_pops[l1]) &
      (data_versus$exclude == competing_pops[ex])
    )
    index_model_B <- which(
    (data_versus$left_1 == competing_pops[ex]) &
      (data_versus$exclude == competing_pops[l1])
    )
    ## cat(sprintf( "Scores %d \t %d \n", data_versus$score_stef[index_model_A], data_versus$score_stef[index_model_B]))
    if (data_versus$score_stef[index_model_A] > data_versus$score_stef[index_model_B]) {
      temp_result_versus <- c(competing_pops[l1], competing_pops[ex], 1)
    } else if (data_versus$score_stef[index_model_A] < data_versus$score_stef[index_model_B]) {
      temp_result_versus <- c(competing_pops[l1], competing_pops[ex], -1)
    } else {
      temp_result_versus <- c(competing_pops[l1], competing_pops[ex], 0)
    }
    results_versus <- rbind(results_versus, temp_result_versus)
  }
}
names(results_versus) <- c("left", "exclude", "score")
results_versus$left <- factor(results_versus$left, levels = competing_pops)
results_versus$exclude <- factor(results_versus$exclude, levels = competing_pops)

# ** Plotting
