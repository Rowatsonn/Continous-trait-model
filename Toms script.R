## Load in required packages ##
library(tidyverse)
library(ggplot2)
##

###### Save useful functions ##
prob_to_log <- function(p){
  return(log(p / (1 - p)))
}

log_to_prob <- function(logit){
  return(1 / (1 + exp(-logit)))
}

## Model function ##
generations = 150
N = 500
P = 0.5
run = NA
cost = 0.2 # The cost for all of the traits 
selection_strength = 1
learning_rate = 0.1
min_lifespan = 20
max_lifespan = 20
mu_lifespan = 0
min_baseline = -3
max_baseline = 1
min_learningspeed = 0
max_learningspeed = 1
mu_learningspeed = 0.01
min_allocation = 0
max_allocation = 1
mu_allocation = 0.01
min_generalisation = 0
max_generalisation = 1
mu_generalisation = 0.01
min_brainsize = 1
max_brainsize = 1
mu_brainsize = 0
min_verticaltransmission = 0
max_verticaltransmission = 0
mu_verticaltransmission = 0
vertical_transmission_rate = 0.5
  
output <- tibble(
  total_payoff = NA,
  learning_speed = NA,
  allocation = NA,
  generalisation = NA,
  brain_size = NA,
  vertical_transmission = NA,
  Lifespan = NA, 
  Baseline = NA, 
  prob_A = NA,
  prob_B = NA,
  max_A = NA,
  max_B = NA,
  run = run,
  generation = NA,
  P = NA,
)

Difference = max_baseline - min_baseline # What is the difference between the minimum and maximum baselines?

start_time <- Sys.time()
print(paste0("Starting simulations where P = ", P))

# Create ranges to sample from

learning_speeds <- seq(from = min_learningspeed, to = max_learningspeed, by = 0.01)
allocations <- seq(from = min_allocation, to = max_allocation, by = 0.01)
generalisations <- seq(from = min_generalisation, to = max_generalisation, by = 0.01)
brain_sizes <- seq(from = min_brainsize, to = max_brainsize, by = 0.01)
vertical_transmissions <- seq(from = min_verticaltransmission, to = max_verticaltransmission, by = 0.01)
lifespans <- seq(from = min_lifespan, to = max_lifespan+1, by = 1) # This workaround is to solve a bug in the sampling

# 1) Sample values for a fresh generation
current_population <- tibble(total_payoff = 0,
                             learning_speed = sample(learning_speeds, size = N, replace = T),
                             allocation = sample(allocations, size = N, replace = T),
                             generalisation = sample(generalisations, size = N, replace = T),
                             brain_size = sample(brain_sizes, size = N, replace = T),
                             vertical_transmission = sample(vertical_transmissions, size = N, replace = T)) %>%
  mutate(Lifespan = sample(lifespans, size = N, replace = T),
         Lifespan = ifelse(Lifespan > max_lifespan, Lifespan - 1, Lifespan), # This is part of solving the bug
         Baseline = 1 - (4 * learning_speed),
         prob_A = Baseline * allocation,
         prob_B = Baseline * (1 - allocation),
         max_A = ((brain_size * 10) * allocation) + 1,
         max_B = ((brain_size * 10) * (1 - allocation) + 1),
         max_A = 2000, ## NOTE. If we want to allow allocation strategy to also increase learning capacity, delete these lines.
         max_B = 2000)

# 2) Start the model for the current generation

for(gen in 1:generations){
  
  # A_prob <- current_population$prob_A
  # B_prob <- current_population$prob_B
  # A_max <- current_population$max_A
  # B_max <- current_population$max_B
  
  for(agent in 1:N) { # For each agent in the population.
    
    # Pull out their information
    A_prob <- current_population$prob_A[agent]
    B_prob <- current_population$prob_B[agent]
    A_max <- current_population$max_A[agent]
    B_max <- current_population$max_B[agent]
    generalisation <- current_population$generalisation[agent]
    learning_speed <- current_population$learning_speed[agent]
    allocation <- current_population$allocation[agent]
    L <- current_population$Lifespan[agent]
    
    # Reset payoffs
    total_payoff <- 0
    Tasks <- sample(c("A", "B"), size = L, replace = T, prob = c(P, 1-P))
    
    for(task in Tasks){
      Success <- ifelse(task == "A", rbinom(1, size = 1, prob = log_to_prob(A_prob)), rbinom(1, size = 1, prob = log_to_prob(B_prob)))
      total_payoff <- total_payoff + Success
      
      # Learning about the tasks
      if(task == "A"){
        A_prob <- A_prob + (learning_rate * learning_speed * allocation)
        B_prob <- B_prob + (learning_rate * learning_speed * allocation * generalisation)
      } else if(task == "B"){
        B_prob <- B_prob + (learning_rate * learning_speed * (1 - allocation))
        A_prob <- A_prob + (learning_rate * learning_speed * (1 - allocation) * generalisation)
      }
      
      # Constrain probabilities based on brain size allocation
      A_prob <- min(A_prob, A_max)
      B_prob <- min(B_prob, B_max)  
    } # End of agent's lifespan
    
    # Record the agent's payoff and update their A and B probs 
    current_population$total_payoff[agent] <- current_population$total_payoff[agent] + total_payoff
    current_population$prob_A[agent] <- A_prob
    current_population$prob_B[agent] <- B_prob
    
  } # End of agent loop
  
  # 4) Apply the fitness penalties
  # NOTE - Need to include brain size if we allow this to evolve also
  
  current_population$total_payoff <- current_population$total_payoff - (current_population$learning_speed * cost + current_population$generalisation * cost + current_population$vertical_transmission * cost) * current_population$Lifespan
  
  # Add extra fitness in case of negative fitness values
  current_population$total_payoff <- (current_population$total_payoff + 10) * selection_strength
  
  end_population <- current_population %>%
    mutate(run = run, generation = gen, P = P)
  
  output <- rbind(output, end_population)
  
  # 5) Sample for the next generation
  
  population_payoff <- sum(current_population$total_payoff)
  
  new_population <- current_population[
    sample(
      nrow(current_population),
      N, # Sample size
      replace = T,
      prob = current_population$total_payoff / population_payoff
    ),
  ]
  
  new_population$total_payoff <- 0
  
  ### Vertical transmission ### Figure out the bonus inheritance to add later.
  
  A_bonus <- new_population$prob_A * (new_population$vertical_transmission * vertical_transmission_rate)
  B_bonus <- new_population$prob_B * (new_population$vertical_transmission * vertical_transmission_rate)
  
  # 6) Mutation. (maybe could be optimised somehow, not sure how)
  mutation_learning <- runif(N, 0, 1) < mu_learningspeed
  mutation_allocation <- rep(sample(c(TRUE, FALSE), N, prob = c(mu_allocation, 1-mu_allocation), replace = T)) 
  mutation_generalisation <- rep(sample(c(TRUE, FALSE), N, prob = c(mu_generalisation, 1-mu_generalisation), replace = T))
  mutation_brain <- rep(sample(c(TRUE, FALSE), N, prob = c(mu_brainsize, 1-mu_brainsize), replace = T))
  mutation_vertical <- rep(sample(c(TRUE, FALSE), N, prob = c(mu_verticaltransmission, 1-mu_verticaltransmission), replace = T))      
  mutation_lifespan <- rep(sample(c(TRUE, FALSE), N, prob = c(mu_lifespan, 1-mu_lifespan), replace = T))
  
  ## Mutate learning speed ## 
  new_population$learning_speed[mutation_learning] <- round(new_population$learning_speed[mutation_learning] + rnorm(sum(mutation_learning), mean = 0, sd = 0.1),2)
  new_population$learning_speed <- pmin(pmax(new_population$learning_speed, min_learningspeed), max_learningspeed) # Limit learning_speed within the bounds        
  
  ## Mutate allocation ##
  new_population$allocation[mutation_allocation] <- round(new_population$allocation[mutation_allocation] + rnorm(sum(mutation_allocation), mean = 0, sd = 0.1),2)
  new_population$allocation <- pmin(pmax(new_population$allocation, min_allocation), max_allocation)
  
  ## Mutate generalisation ##
  new_population$generalisation[mutation_generalisation] <- round(new_population$generalisation[mutation_generalisation] + rnorm(sum(mutation_generalisation), mean = 0, sd = 0.1),2)
  new_population$generalisation <- pmin(pmax(new_population$generalisation, min_generalisation), max_generalisation)
  
  ## Mutate brain size ##
  new_population$brain_size[mutation_brain] <- round(new_population$brain_size[mutation_brain] + rnorm(sum(mutation_brain), mean = 0, sd = 0.1),2)
  new_population$brain_size <- pmin(pmax(new_population$brain_size, min_brainsize), max_brainsize) 
  
  ## Mutate vertical transmission ## 
  new_population$vertical_transmission[mutation_vertical] <- round(new_population$vertical_transmission[mutation_vertical] + rnorm(sum(mutation_vertical), mean = 0, sd = 0.1),2)
  new_population$vertical_transmission <- pmin(pmax(new_population$vertical_transmission, min_verticaltransmission), max_verticaltransmission) 
  
  ## Mutate Lifespan ## 
  new_population$Lifespan[mutation_lifespan] <- round(new_population$Lifespan[mutation_lifespan] + rnorm(sum(mutation_lifespan), mean = 0, sd = 3))
  new_population$Lifespan <- pmin(pmax(new_population$Lifespan, min_lifespan), max_lifespan) 
  
  ## Recalculate the A max and B max based on brain size and allocation strategy
  
  new_population$max_A <- new_population$brain_size * new_population$allocation
  new_population$max_B <- new_population$brain_size * (1 - new_population$allocation)
  
  # Finally, recalculate the prob values from the mutations
  new_population <- new_population %>%
    mutate(Baseline = 1 - (4 * learning_speed),
           prob_A = Baseline * allocation,
           prob_B = Baseline * (1 - allocation),
           max_A = ((brain_size * 10) * allocation) + 1,
           max_B = ((brain_size * 10) * (1 - allocation) + 1),
           max_A = 2000, ## NOTE. When we want to allow allocation strategy to also increase learning capacity. Delete these lines.
           max_B = 2000)
  
  # Then, add the bonus from vertical transmission to prob_A and prob_B
  
  new_population$prob_A <- new_population$prob_A + A_bonus
  new_population$prob_B <- new_population$prob_B + B_bonus
  
  current_population <- new_population
  
}

end_time <- Sys.time()
end_time - start_time

