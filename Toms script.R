rm(list=ls())
hack <- FALSE

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
generations = 200
N = 500
P = 0.5
cost = 0.02
Baseline = -1

learning_rate = 5

min_lifespan = 20
max_lifespan = 20

min_learningspeed = 0
max_learningspeed = 0

min_allocation = 0
max_allocation = 1

min_generalisation = 0
max_generalisation = 0

min_verticaltransmission = 0
max_verticaltransmission = 0

vertical_transmission_rate = 0.5

output <- data.frame(
  generation = c(1:generations),
  allocation = vector(length=generations),
  learning_speed = vector(length=generations),
  generalisation = vector(length=generations)
)

start_time <- Sys.time()
print(paste0("Starting simulations where P = ", P))


for(gen in 1:generations){
  
  # 1) Sample values for a fresh generation

  if (gen==1) {
    fitness = rep(0.001, N)
    total_payoff = rep(0, N)
    learning_speed = rep(min_learningspeed, N)
    allocation = rep(min_allocation, N)
    generalisation = rep(min_generalisation, N)
    vertical_transmission = rep(min_verticaltransmission, N)
    Lifespan = rep(min_lifespan, N)
  } else {
    parents <- sample(c(1:N), replace=T, prob=fitness)
    
    fitness = rep(0.001, N)
    total_payoff = rep(0, N)
    learning_speed <- pmin(pmax(learning_speed[parents] + rnorm(N, mean = 0, sd = 0.1), min_learningspeed), max_learningspeed)
    allocation <- pmin(pmax(allocation[parents] + rnorm(N, mean = 0, sd = 0.1), min_allocation), max_allocation)
    generalisation <- pmin(pmax(generalisation[parents] + rnorm(N, mean = 0, sd = 0.1), min_generalisation), max_generalisation)
    vertical_transmission <- pmin(pmax(vertical_transmission[parents] + rnorm(N, mean = 0, sd = 0.1), min_verticaltransmission), max_verticaltransmission)
    Lifespan <- round(pmin(pmax(Lifespan[parents] + rnorm(N, mean = 0, sd = 0.1), min_lifespan), max_lifespan))
  }
  prob_A = Baseline + 4*(allocation-0.5) - (4*learning_speed)
  prob_B = Baseline + 4*(0.5-allocation) - (4*learning_speed)
  
  
  # 2) Solve tasks
  if (hack) {
    Tasks <- rbinom(N, Lifespan, P)
    total_payoff <- rbinom(N, Tasks, log_to_prob(prob_A)) + rbinom(N, Lifespan-Tasks, log_to_prob(prob_B))
    prob_A <- prob_A + Tasks*learning_rate*learning_speed*allocation + (Lifespan-Tasks)*(learning_rate * learning_speed * (1-allocation) * generalisation)
    prob_B <- prob_B + (Lifespan-Tasks)*(learning_rate * learning_speed * (1-allocation)) + Tasks*(learning_rate * learning_speed * (1-allocation)*generalisation)
  } else {
    for(agent in 1:N) {
      Tasks <- rbinom(Lifespan[agent], 1, P)
      for(task in Tasks){
        total_payoff[agent] <- rbinom(1, 1, prob = (task*log_to_prob(prob_A[agent])) + (1-task)*log_to_prob(prob_B[agent]))
        
        prob_A[agent] <- prob_A[agent] + task*(learning_rate * learning_speed[agent] * allocation[agent]) + (1-task)*(learning_rate * learning_speed[agent] * (1-allocation[agent]) * generalisation[agent])
        prob_B[agent] <- prob_B[agent] + (1-task)*(learning_rate * learning_speed[agent] * (1-allocation[agent])) + task*(learning_rate * learning_speed[agent] * allocation[agent]*generalisation[agent])
      }
    }
  }
  
  # 3) Calculate fitness
  total_payoff <- pmax(total_payoff - cost*(learning_speed + generalisation + vertical_transmission) * Lifespan, 0.001)
  fitness <- total_payoff/sum(total_payoff)
  
  # 4) Save data
  output$allocation[gen] <- mean(allocation)
  output$learning_speed[gen] <- mean(learning_speed)
  output$generalisation[gen] <- mean(generalisation)
  
}

end_time <- Sys.time()
end_time - start_time

plot(output$allocation~output$generation, type="b", ylim=c(0, 1))
lines(output$learning_speed~output$generation, type="b", col="maroon")
lines(output$generalisation~output$generation, type="b", col="cornflowerblue")

