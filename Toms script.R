hack <- TRUE

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
current_population <- data.frame(
  total_payoff = rep(0, N),
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
  max_A = rep(2000, N), ## NOTE. If we want to allow allocation strategy to also increase learning capacity, delete these lines.
  max_B = rep(2000,N)
)

# 2) Start the model for the current generation

for(gen in 1:generations){
  
  A_prob <- current_population$prob_A
  B_prob <- current_population$prob_B
  A_max <- current_population$max_A
  B_max <- current_population$max_B
  generalisation <- current_population$generalisation
  learning_speed <- current_population$learning_speed
  allocation <- current_population$allocation
  L <- current_population$Lifespan
  
  current_population$total_payoff <- rep(0, N)
  
  if (hack) {
    Tasks <- rbinom(N, L, P)
    current_population$total_payoff <- rbinom(N, Tasks, log_to_prob(A_prob)) + rbinom(N, L-Tasks, log_to_prob(B_prob))
    A_prob <- min(A_prob + Tasks*(learning_rate * learning_speed * allocation) + (L-Tasks)*(learning_rate * learning_speed * (1-allocation) * generalisation), A_max)
    B_prob <- min(B_prob + (L-Tasks)*(learning_rate * learning_speed * allocation * generalisation) + Tasks*(learning_rate * learning_speed * (1-allocation)), B_max)
  } else {
    for(agent in 1:N) { # For each agent in the population.
      Tasks <- rbinom(L[agent], 1, P)
      for(task in Tasks){
        current_population$total_payoff[agent] <- rbinom(1, 1, prob = (task*log_to_prob(A_prob[agent])) + (1-task)*log_to_prob(B_prob[agent]))
        
        A_prob[agent] <- min(A_prob[agent] + task*(learning_rate * learning_speed[agent] * allocation[agent]) + (1-task)*(learning_rate * learning_speed[agent] * (1-allocation[agent]) * generalisation[agent]), A_max)
        B_prob[agent] <- min(B_prob[agent] + (1-task)*(learning_rate * learning_speed[agent] * allocation[agent] * generalisation[agent]) + task*(learning_rate * learning_speed[agent] * (1-allocation[agent])), B_max)
      }
    }
  }
  
  # 4) Apply the fitness penalties
  # NOTE - Need to include brain size if we allow this to evolve also
  
  current_population$total_payoff <- pmax(current_population$total_payoff - cost*(current_population$learning_speed + current_population$generalisation + current_population$vertical_transmission) * current_population$Lifespan, 0.001)
  
  # Add extra fitness in case of negative fitness values
  
  # selection strength won't do anything here
  # current_population$total_payoff <- (current_population$total_payoff + 10) * selection_strength 

  
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

  new_population$learning_speed <- pmin(pmax(new_population$learning_speed + rnorm(N, mean = 0, sd = 0.1), min_learningspeed), max_learningspeed)
  new_population$allocation <- pmin(pmax(new_population$allocation + rnorm(N, mean = 0, sd = 0.1), min_allocation), max_allocation)
  new_population$generalisation <- pmin(pmax(new_population$generalisation + rnorm(N, mean = 0, sd = 0.1), min_generalisation), max_generalisation)
  new_population$brain_size <- pmin(pmax(new_population$brain_size + rnorm(N, mean = 0, sd = 0.1), min_brainsize), max_brainsize)
  new_population$vertical_transmission <- pmin(pmax(new_population$vertical_transmission + rnorm(N, mean = 0, sd = 0.1), min_verticaltransmission), max_verticaltransmission)
  new_population$Lifespan <- pmin(pmax(new_population$Lifespan + rnorm(N, mean = 0, sd = 0.1), min_lifespan), max_lifespan)
  
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

