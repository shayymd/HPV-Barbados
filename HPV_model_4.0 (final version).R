#' Title: Age and gender stratified HPV SIS transmission model
#' Author: Shay Stabler-Morris & Kevin van Zandvoort
#' Date: 22/8/2021
#' Description: A Susceptible-Infectious-Susceptible transmission model with 3 age groups and 2 genders that is implemented
#'  in R using deSolve. Model is deterministic and uses differential equations (runs in continuous time).

#' ----- Step 0 - load required packages
library("tidyverse")
library("deSolve") 
library("optimx")
library("ggpubr")


options(dplyr.summarise.inform = FALSE)

#' ----- Step 1 - specify demographics
#' Total Population = 287371
popsize = 239206 # population size of Barbados minus the 0-14 year olds (UN 2020 estimates)

#' Assume an average life expectancy of 79 years (World Bank 2021)
life_exp = 79

#' lowest two age groups spend, on average, 10 years in each age group. Oldest age group spends life_exp-10-10-15 years in that age group. 
age_duration = c(10, 10, life_exp-35)
N = age_duration/sum(age_duration) * popsize # should make N[1] and N[2] the same size.

#' Assume a rectangular age distribution, where deaths only occur in the oldest age group.
#' age_out is equal to 1 over the total duration of their age group.
age_out = 1/(age_duration*365)

#' To keep a stable population size, births equal the deaths occurring in the old age group
#' age_in is created to be an empty numeric vector with the same size as the age_out vector specified above (3 age groups)
age_in = numeric(length(age_out))
age_in[2:3] = age_out[1:2]
age_in[1] = age_out[3]

#' ----- Step 2 - specify disease and transmission parameters
#' Around 90% of infections clear within 2 years
#' The cumulative incidence 'risk' of clearing infection is 90% after two years
risk = 0.9
time = 365*2 #days
rate = -log(1-risk)/time 
dur_infectiousness = 1/rate #days, irrespective of age

# proportion of high risk in each age group
prop_high_risk = c(0.339, 0.252, 0.128)

# The average number of partnerships per year in the low- and high-activity group
partnerships_year = readRDS("./partnerships_year.RDS")

# Proportion of partnerships that are made by contactors of age j with contactees of age i

contacts_age = matrix(c(
  0.822, 0.089, 0.089,
  0.089, 0.911, 0,
  0.089, 0, 0.911), ncol=3, dimnames=list(contactee_age = c("15-24", "25-34", "35+"),
                                  contactor_age = c("15-24", "25-34", "35+")))

#contacts_age = matrix(c(
 # 1, 0, 0,
  #0, 1, 0,
  #0, 0, 1), ncol=3, dimnames=list(contactee_age = c("15-24", "25-34", "35+"),
   #                                       contactor_age = c("15-24", "25-34", "35+")))


#' A function that takes the proportion of partnerships offered by the high-activity group, and calculates the proportion of
#' contacts made by the low/high-activity group for a given Q-value
getMixingForQ = function(prop_high, Q = 0){
  if(Q > 1) stop("max Q is 1")
  if(Q < 0) warning("results for negative Q may be invalid")
  
  low = 1-prop_high
  high = prop_high
  
  x = -(Q-1)
  
  y = matrix(c(
    1-high*x, high*x,
    low*x, 1-low*x
  ), ncol=2, dimnames=list(contactee = c("low", "high"), contactor = c("low", "high")))
  
  if(any(y < 0)) warning("some values are negative")
  
  return(y)
}


#' Set a value for Q (see Moodle MD08 7.3)
Q = 0
  
# Average number of annual partnerships in each group
partnerships_by_group = partnerships_year[c(which(partnerships_year$partner_year_cat == "Low"),
                                              which(partnerships_year$partner_year_cat == "High")),] %>%
  .[["partner_year_mean"]]
names(partnerships_by_group) = c("low", "high")
  
# convert to average number of partnerships per day
partnerships_by_group = partnerships_by_group/365
  
#calculate the proportion of partnerships offered by low/high group in age group 1
prop_partners_1 = partnerships_by_group * c(1 - prop_high_risk[1], prop_high_risk[1])/sum(
  partnerships_by_group * c(1 - prop_high_risk[1], prop_high_risk[1]))
  
#calculate proportion of partnerships offered by low/high group in age group 2
prop_partners_2 = partnerships_by_group * c(1 - prop_high_risk[2], prop_high_risk[2])/sum(
  partnerships_by_group * c(1 - prop_high_risk[2], prop_high_risk[2]))
  
#calculate proportion of partnerships offered by low/high group in age group 3
prop_partners_3 = partnerships_by_group * c(1 - prop_high_risk[3], prop_high_risk[3])/sum(
  partnerships_by_group * c(1 - prop_high_risk[3], prop_high_risk[3]))

#check values for the provided Q value  
getMixingForQ(prop_partners_1[2], Q)
getMixingForQ(prop_partners_2[2], Q)
getMixingForQ(prop_partners_3[2], Q)
  
#' Calculate beta-values (beta_12_lh is the daily number of partners made by people in the low-activity group in age
#'  group 1 with people in the high-activity group in age-group 2 )

# contactor: 15-24 year old low
beta_11_ll = contacts_age[1, 1] * partnerships_by_group["low"] * getMixingForQ(prop_partners_1[2], Q)["low", "low"]
beta_11_lh = contacts_age[1, 1] * partnerships_by_group["low"] * getMixingForQ(prop_partners_1[2], Q)["high", "low"]
beta_12_ll = contacts_age[2, 1] * partnerships_by_group["low"] * getMixingForQ(prop_partners_2[2], Q)["low", "low"]
beta_12_lh = contacts_age[2, 1] * partnerships_by_group["low"] * getMixingForQ(prop_partners_2[2], Q)["high", "low"]
beta_13_ll = contacts_age[3, 1] * partnerships_by_group["low"] * getMixingForQ(prop_partners_3[2], Q)["low", "low"]
beta_13_lh = contacts_age[3, 1] * partnerships_by_group["low"] * getMixingForQ(prop_partners_3[2], Q)["high", "low"]
# contactor: 15-24 year old high
beta_11_hl = contacts_age[1, 1] * partnerships_by_group["high"] * getMixingForQ(prop_partners_1[2], Q)["low", "high"]
beta_11_hh = contacts_age[1, 1] * partnerships_by_group["high"] * getMixingForQ(prop_partners_1[2], Q)["high", "high"]
beta_12_hl = contacts_age[2, 1] * partnerships_by_group["high"] * getMixingForQ(prop_partners_2[2], Q)["low", "high"]
beta_12_hh = contacts_age[2, 1] * partnerships_by_group["high"] * getMixingForQ(prop_partners_2[2], Q)["high", "high"]
beta_13_hl = contacts_age[3, 1] * partnerships_by_group["high"] * getMixingForQ(prop_partners_3[2], Q)["low", "high"]
beta_13_hh = contacts_age[3, 1] * partnerships_by_group["high"] * getMixingForQ(prop_partners_3[2], Q)["high", "high"]

# contactor: 25-34 year old low
beta_21_ll = contacts_age[1, 2] * partnerships_by_group["low"] * getMixingForQ(prop_partners_1[2], Q)["low", "low"]
beta_21_lh = contacts_age[1, 2] * partnerships_by_group["low"] * getMixingForQ(prop_partners_1[2], Q)["high", "low"]
beta_22_ll = contacts_age[2, 2] * partnerships_by_group["low"] * getMixingForQ(prop_partners_2[2], Q)["low", "low"]
beta_22_lh = contacts_age[2, 2] * partnerships_by_group["low"] * getMixingForQ(prop_partners_2[2], Q)["high", "low"]
beta_23_ll = contacts_age[3, 2] * partnerships_by_group["low"] * getMixingForQ(prop_partners_3[2], Q)["low", "low"]
beta_23_lh = contacts_age[3, 2] * partnerships_by_group["low"] * getMixingForQ(prop_partners_3[2], Q)["high", "low"]
# contactor: 25-34 year old high
beta_21_hl = contacts_age[1, 2] * partnerships_by_group["high"] * getMixingForQ(prop_partners_1[2], Q)["low", "high"]
beta_21_hh = contacts_age[1, 2] * partnerships_by_group["high"] * getMixingForQ(prop_partners_1[2], Q)["high", "high"]
beta_22_hl = contacts_age[2, 2] * partnerships_by_group["high"] * getMixingForQ(prop_partners_2[2], Q)["low", "high"]
beta_22_hh = contacts_age[2, 2] * partnerships_by_group["high"] * getMixingForQ(prop_partners_2[2], Q)["high", "high"]
beta_23_hl = contacts_age[3, 2] * partnerships_by_group["high"] * getMixingForQ(prop_partners_3[2], Q)["low", "high"]
beta_23_hh = contacts_age[3, 2] * partnerships_by_group["high"] * getMixingForQ(prop_partners_3[2], Q)["high", "high"]

# contactor: 35+ year old low
beta_31_ll = contacts_age[1, 3] * partnerships_by_group["low"] * getMixingForQ(prop_partners_1[2], Q)["low", "low"]
beta_31_lh = contacts_age[1, 3] * partnerships_by_group["low"] * getMixingForQ(prop_partners_1[2], Q)["high", "low"]
beta_32_ll = contacts_age[2, 3] * partnerships_by_group["low"] * getMixingForQ(prop_partners_2[2], Q)["low", "low"]
beta_32_lh = contacts_age[2, 3] * partnerships_by_group["low"] * getMixingForQ(prop_partners_2[2], Q)["high", "low"]
beta_33_ll = contacts_age[3, 3] * partnerships_by_group["low"] * getMixingForQ(prop_partners_3[2], Q)["low", "low"]
beta_33_lh = contacts_age[3, 3] * partnerships_by_group["low"] * getMixingForQ(prop_partners_3[2], Q)["high", "low"]
# contactor: 35+ year old high
beta_31_hl = contacts_age[1, 3] * partnerships_by_group["high"] * getMixingForQ(prop_partners_1[2], Q)["low", "high"]
beta_31_hh = contacts_age[1, 3] * partnerships_by_group["high"] * getMixingForQ(prop_partners_1[2], Q)["high", "high"]
beta_32_hl = contacts_age[2, 3] * partnerships_by_group["high"] * getMixingForQ(prop_partners_2[2], Q)["low", "high"]
beta_32_hh = contacts_age[2, 3] * partnerships_by_group["high"] * getMixingForQ(prop_partners_2[2], Q)["high", "high"]
beta_33_hl = contacts_age[3, 3] * partnerships_by_group["high"] * getMixingForQ(prop_partners_3[2], Q)["low", "high"]
beta_33_hh = contacts_age[3, 3] * partnerships_by_group["high"] * getMixingForQ(prop_partners_3[2], Q)["high", "high"]
  
contact_matrix = matrix(c(
  0, 0, 0, 0, 0, 0, beta_11_ll,  beta_11_lh,     beta_12_ll,     beta_12_lh,     beta_13_ll,     beta_13_lh,
  0, 0, 0, 0, 0, 0, beta_11_hl,  beta_11_hh,     beta_12_hl,     beta_12_hh,     beta_13_hl,     beta_13_hh,
  0, 0, 0, 0, 0, 0, beta_21_ll,  beta_21_lh,     beta_22_ll,     beta_22_lh,     beta_23_ll,     beta_23_lh,
  0, 0, 0, 0, 0, 0, beta_21_hl,  beta_21_hh,     beta_22_hl,     beta_22_hh,     beta_23_hl,     beta_23_hh,
  0, 0, 0, 0, 0, 0, beta_31_ll,  beta_31_lh,     beta_32_ll,     beta_32_lh,     beta_33_ll,     beta_33_lh,
  0, 0, 0, 0, 0, 0, beta_31_hl,  beta_31_hh,     beta_32_hl,     beta_32_hh,     beta_33_hl,     beta_33_hh,
  
  beta_11_ll,   beta_11_lh,     beta_12_ll,     beta_12_lh,     beta_13_ll,     beta_13_lh, 0, 0, 0, 0, 0, 0,
  beta_11_hl,   beta_11_hh,     beta_12_hl,     beta_12_hh,     beta_13_hl,     beta_13_hh, 0, 0, 0, 0, 0, 0,
  beta_21_ll,   beta_21_lh,     beta_22_ll,     beta_22_lh,     beta_23_ll,     beta_23_lh, 0, 0, 0, 0, 0, 0,
  beta_21_hl,   beta_21_hh,     beta_22_hl,     beta_22_hh,     beta_23_hl,     beta_23_hh, 0, 0, 0, 0, 0, 0,
  beta_31_ll,   beta_31_lh,     beta_32_ll,     beta_32_lh,     beta_33_ll,     beta_33_lh, 0, 0, 0, 0, 0, 0,
  beta_31_hl,   beta_31_hh,     beta_32_hl,     beta_32_hh,     beta_33_hl,     beta_33_hh, 0, 0, 0, 0, 0, 0),
  nrow=2*3*2, ncol=2*3*2,
  dimnames = list(contactee = expand.grid(risk = c("low-risk", "high-risk"), age = c("15-24", "25-34", "35+"), gender = c("females", "males")) %>%
                    as_tibble() %>% mutate(name = paste(risk, gender, age)) %>% .[["name"]],
                  contactor = expand.grid(risk = c("low-risk", "high-risk"), age = c("15-24", "25-34", "35+"), gender = c("females", "males")) %>%
                    as_tibble() %>% mutate(name = paste(risk, gender, age)) %>% .[["name"]]))
  
#' ----- Step 3 - prepare the model
#' When fitting the model, you don't want to introduce vaccination (model the prevalence without vaccination)
parameters = list(
  "age_in" = age_in,
  "age_out" = age_out,
  "contact_matrix" = contact_matrix, 
  "gamma" = 1/dur_infectiousness, # rate at which infectious become susceptible again
  "delta" = 1/(365*14), # rate at which vaccinated become susceptible again -> 14 years, not 20 years
  "vacc_cov_f" = 0, # vaccine coverage in girls; this parameter will change during exploration
  "vacc_cov_m" = 0, # vaccine coverage in boys; to be 90% of girl coverage
  "vacc_eff" = 1 # vaccine efficacy; range from 89-100%
)
  
#' Initial state for the compartments in our model (i.e., the total number of people in each compartment
#'  at the start of the simulation)
state = c(
  "Sf1l" = N[1]*0.5*(1-prop_high_risk[1]) - 50,
  "Sf1h" = N[1]*0.5*(prop_high_risk[1]) - 50,
  "Sf2l" = N[2]*0.5*(1-prop_high_risk[2]) - 50,
  "Sf2h" = N[2]*0.5*(prop_high_risk[2]) - 50,
  "Sf3l" = N[3]*0.5*(1-prop_high_risk[3]) - 50,
  "Sf3h" = N[3]*0.5*(prop_high_risk[3]) - 50,
  "Sm1l" = N[1]*0.5*(1-prop_high_risk[1]) - 50,
  "Sm1h" = N[1]*0.5*(prop_high_risk[1]) - 50,
  "Sm2l" = N[2]*0.5*(1-prop_high_risk[2]) - 50,
  "Sm2h" = N[2]*0.5*(prop_high_risk[2]) - 50,
  "Sm3l" = N[3]*0.5*(1-prop_high_risk[3]) - 50,
  "Sm3h" = N[3]*0.5*(prop_high_risk[3]) - 50,
  "If1l" = 50,
  "If1h" = 50,
  "If2l" = 50,
  "If2h" = 50,
  "If3l" = 50,
  "If3h" = 50,
  "Im1l" = 50,
  "Im1h" = 50,
  "Im2l" = 50,
  "Im2h" = 50,
  "Im3l" = 50,
  "Im3h" = 50,
  "Vf1l" = 0,
  "Vf1h" = 0,
  "Vf2l" = 0,
  "Vf2h" = 0,
  "Vf3l" = 0,
  "Vf3h" = 0,
  "Vm1l" = 0,
  "Vm1h" = 0,
  "Vm2l" = 0,
  "Vm2h" = 0,
  "Vm3l" = 0,
  "Vm3h" = 0,
  "Cf1l" = 0,
  "Cf1h" = 0,
  "Cf2l" = 0,
  "Cf2h" = 0,
  "Cf3l" = 0,
  "Cf3h" = 0,
  "Cm1l" = 0,
  "Cm1h" = 0,
  "Cm2l" = 0,
  "Cm2h" = 0,
  "Cm3l" = 0,
  "Cm3h" = 0)
  
#' We now specify our model (New model using linear algebra to speed up calculations)
#' - this is an R function that takes 3 arguments:
#'   - t: the current day
#'   - state: a vector with the current state of each compartment
#'   - pars: the list with parameter values
    
#' parameters that will be used in the model
Ngroups = 2 * 3 * 2 #number of strata in each compartment
prop_high_risk_remain = prop_high_risk[2:3]/prop_high_risk[1:2] #proportion that will remain in the high-risk group for each age
  
SIS2 = function(t, state, pars){
  #' parameters passed to the model
  age_in = rep(rep(pars$age_in, each=2), 2)
  age_out = rep(rep(pars$age_out, each=2), 2)
  contact_matrix = pars[["contact_matrix"]]
  gamma = pars[["gamma"]]
  delta = pars[["delta"]]
  vacc_cov_f = pars[["vacc_cov_f"]]
  vacc_cov_m = pars[["vacc_cov_m"]]
  vacc_eff = pars[["vacc_eff"]]
  beta_l = pars[["beta_l"]]
  beta_h = pars[["beta_h"]]
  
  #' We also read the state value of all compartments as vectors
  S = state[1:Ngroups + Ngroups*0]
  I = state[1:Ngroups + Ngroups*1]
  V = state[1:Ngroups + Ngroups*2]
  C = state[1:Ngroups + Ngroups*3]
  
  N = S+I+V#+C
  #' N in oldest age group
  N2 = sum(N[c(5, 6, 11, 12)])
  
  #' create a vector of beta values (beta_l will be multiplied with infected people in the low-activity group)
  #'  beta_h will be multiplied with infected people in the high-activity group
  beta = rep(c(beta_l, beta_h), 6)
  
  #' Consider whether you need to use frequency- or density dependent transmission
  
  #' We sum I and C, multiply this by beta, perform a matrix multiplication with the contact matrix, and take the
  #'  product of that vector with the vector of susceptibles
  new_infections = (((I/N)*beta) %*% contact_matrix) * S
  
  #all age groups
  dS = (
    -new_infections
    #gamma is multiplied by 0.9 for the first 6 I-values (females), and 1 for the next 6 (males)
    +I*gamma #* c(rep(0.9, 6), rep(1, 6))
    +V*delta
    -age_out*S
    +age_in*c(0, 0, S[1:4], 0, 0, S[7:10])*rep(c(0, 0, 1, prop_high_risk_remain[1], 1, prop_high_risk_remain[2]), 2))
  
  #high risk ageing into low risk
  dS[c(3,5,9,11)] = dS[c(3,5,9,11)] + age_in[c(3,5,9,11)]*(1-rep(c(prop_high_risk_remain[1], prop_high_risk_remain[2]), 2))*S[c(2,4,8,10)]
  
  #newborns getting (not) vaccinated
  dS[c(1,2,7,8)] = dS[c(1,2,7,8)] + age_in[c(1,2,7,8)]*N2*0.5*
    c(1-prop_high_risk[1], prop_high_risk[1], 1-prop_high_risk[1], prop_high_risk[1])*
    c(1-vacc_cov_f*vacc_eff, (1-vacc_cov_f*vacc_eff), 1-vacc_cov_m*vacc_eff, (1-vacc_cov_m*vacc_eff))
  
  dI = (
    #transmissions (new infections)
    +new_infections
    #transitions (recovery of infection)
    -I*gamma 
    #demographic changes (births and ageing)
    -age_out*I
    +age_in*c(0, 0, I[1:4], 0, 0, I[7:10])*rep(c(0, 0, 1, prop_high_risk_remain[1], 1, prop_high_risk_remain[2]), 2))
  
  #high risk ageing into low risk
  dI[c(3,5,9,11)] = dI[c(3,5,9,11)] + age_in[c(3,5,9,11)]*(1-rep(c(prop_high_risk_remain[1], prop_high_risk_remain[2]), 2))*I[c(2,4,8,10)]
  
  dV = (
    #transitions (recovery of infection)
    -V*delta 
    #demographic changes (births and ageing)
    -age_out*V
    +age_in*c(0, 0, V[1:4], 0, 0, V[7:10])*rep(c(0, 0, 1, prop_high_risk_remain[1], 1, prop_high_risk_remain[2]), 2))
  
  #high risk ageing into low risk
  dV[c(3,5,9,11)] = dV[c(3,5,9,11)] + age_in[c(3,5,9,11)]*(1-rep(c(prop_high_risk_remain[1], prop_high_risk_remain[2]), 2))*V[c(2,4,8,10)]
  
  #newborns getting vaccinated
  dV[c(1,2,7,8)] = dV[c(1,2,7,8)] + age_in[c(1,2,7,8)]*N2*0.5*
    c(1-prop_high_risk[1], prop_high_risk[1], 1-prop_high_risk[1], prop_high_risk[1])*
    c(vacc_cov_f*vacc_eff, vacc_cov_f*vacc_eff, vacc_cov_m*vacc_eff, vacc_cov_m*vacc_eff)
  
  dC = ( # incidence of new infections
    #transitions (recovery of infection)
    #gamma is multiplied by 0.9 for the first 6 I-values (females), and 1 for the next 6 (males)
    #+I*gamma*c(rep(0.1, 6), rep(0, 6))
    #demographic changes (births and ageing)
    +new_infections 
    - C
    #-age_out*C
    #+age_in*c(0, 0, C[1:4], 0, 0, C[7:10])*rep(c(0, 0, 1, prop_high_risk_remain[1], 1, prop_high_risk_remain[2]), 2)
    )
  
  #high risk ageing into low risk
  #dC[c(3,5,9,11)] = dC[c(3,5,9,11)] + age_in[c(3,5,9,11)]*(1-rep(c(prop_high_risk_remain[1], prop_high_risk_remain[2]), 2))*C[c(2,4,8,10)]
  
  #' We return the calculated changes to each compartment as a list
  return(list(c(dS, dI, dV, dC)))
}

#' ----- Step 5 - run the model

#' We can run the model using the lsoda function from the deSolve package
#' - lsoda will use the correct ODE solver, but the technical details are not important
#' - we need to provide it with the initial state, times where we require output, model function, and parameter values
#model_output = deSolve::lsoda(state, times, SIS, parameters)

#' Read the prevalence data to be used to fit the model
prevalence_data = readRDS("./prevalence.RDS")

#' For how long to run the model?
inityears = 40
maxyears = 200

#' This function will be used to run the model for a set of parameter values, and calculate the -log-likelihood
modelLL = function(fitparams=c("beta_l" = log(0.00001), "beta_h" = log(0.00002)), model_parameters = parameters){
  
  #' Number of years the model is initially ran for
  years = inityears
  
  #' We exponentiate the values to ensure they are always positive
  beta_l = exp(fitparams["beta_l"])
  beta_h = exp(fitparams["beta_h"])
  
  parameters[["beta_l"]] = beta_l
  parameters[["beta_h"]] = beta_h
  
  #' We set a maxvar value of 50. maxvar will be the highest value of the variance in the susceptibles
  #'  over the last 365 days modelled.
  #' We only use the model results when the variability is less than or equal to 0.5
  #'  This will ensure that the model has reached an equilibrium
  #'  If no equilibrium is reached, we rerun the model, 10 years longer, until an equilibrium is reached
  maxvar = 50
  while(maxvar > 0.5){
    if(years > inityears)
      message(paste0("maxvar > 0.5; rerunning model; years: ", years))
    
    model_output2 = deSolve::lsoda(state, c(1:(365*years)), SIS2, parameters)
    model_output2 = tibble::as_tibble(model_output2)
    model_output2 = model_output2 %>%
      pivot_longer(cols = c("Sf1l", "Sf1h", "Sf2l", "Sf2h", "Sf3l", "Sf3h", "Sm1l", "Sm1h", "Sm2l", "Sm2h", "Sm3l", "Sm3h", 
                            "If1l", "If1h", "If2l", "If2h", "If3l", "If3h", "Im1l", 'Im1h', 'Im2l', "Im2h", "Im3l", "Im3h",
                            "Vf1l", "Vf1h", "Vf2l", "Vf2h", "Vf3l", "Vf3h", "Vm1l", "Vm1h", "Vm2l", "Vm2h", "Vm3l", "Vm3h", 
                            "Cf1l", "Cf1h", "Cf2l", "Cf2h", "Cf3l", "Cf3h", "Cm1l", "Cm1h", "Cm2l", "Cm2h", "Cm3l", "Cm3h"), names_to = "model_compartment") %>%
      # create a new column with only the compartment (without the age group)
      mutate(compartment = gsub("\\d+", "", model_compartment)) %>%
      # create a new column with only the age group (without the compartment)
      mutate(age_group = gsub("\\D+", "", model_compartment) %>% as.numeric()) %>%
      mutate(state = substr(compartment, 1, 1), gender = substr(compartment, 2, 2),
             activity = substr(compartment, 3, 3))
    
    #' Calculate the maximum variance in the susceptible females (maximum in age- and activity strata)
    maxvar = model_output2 %>% filter(state == "S" & gender == "f" & time %in% c(max(time):(max(time)-365))) %>% group_by(age_group, state, activity) %>% summarise(var = var(value)) %>% .[["var"]] %>% max
    
    years = years + 10
    if(years > maxyears){
      #Return a high -LL value if takes long to reach equilibrium
      message("years has exceeded maxyears, returning early")
      return(1e6)
    }
  }
  
  #' Overwrite the inityears variable in the global environment (outside the scope of this function)
  assign("inityears", years-10, envir = .GlobalEnv)
  
  #Calculate the log likelihood of the model parameters given the data
  log_likelihood = model_output2 %>%
    filter(time == max(time), gender == "f") %>%
    group_by(age_group, activity) %>%
    mutate(N = sum(value)) %>%
    #' We combine both I and C when calculating prevalence, is this what you want?
    #filter(state %in% c("I", "C")) %>%
    filter(state %in% c("I")) %>%
    summarise(value = sum(value), N=mean(N)) %>%
    mutate(modelled_prev = value/N) %>%
    mutate(age_group = recode(age_group,
                              "1" = "15-24",
                              "2" = "25-34",
                              "3" = "35+")) %>%
    mutate(partner_year = recode(activity,
                                 "l" = "Low",
                                 "h" = "High")) %>%
    select(age_group, partner_year, modelled_prev) %>%
    left_join(prevalence_data %>% select(age_group, partner_year, n, total), by=c("age_group", "partner_year")) %>%
    rowwise %>%
    #the LL is calculated here
    mutate(ll = dbinom(total, n, modelled_prev, log=TRUE)) %>%
    .[["ll"]] %>% sum()
  
  #If log_likelihood cannot be calculated, set to a very low number (parameters don't reproduce the data)
  if(is.na(log_likelihood)) log_likelihood = -1e6
  
  #' Return a message to keep track of the iterations of the model
  message(sprintf("log_likelihood: %s; beta_l: %s; beta_h: %s; years: %s", round(-log_likelihood, 5), round(beta_l, 6), round(beta_h, 6), years-10))
  
  return(-log_likelihood)
}

#' fit the transmission model
#' use optimx to find parameter values using the Nelder-Mead algorithm
#' It minimizes the function modelLL (maximizes the log likelihood, same as minimizing the -ll)
#model_fit = optimx::optimx(
 # par = c("beta_l"=log(0.0001), "beta_h"=log(0.0002)),
  #fn = modelLL,
  #method = c("Nelder-Mead"),
  #control = list(maximize=F, maxit = 100)) #max 100 iterations

# saveRDS(model_fit, "model_fit.RDS")

model_fit = readRDS("model_fit.RDS")


########################################
############### Analysis ###############
########################################

#Run the model with the 'optimal' parameter values
parameters[["beta_l"]] = 0.000154 #exp(model_fit[1, "beta_l"])
parameters[["beta_h"]] = 0.861074 #exp(model_fit[1, "beta_h"])


#' We specify for how many days we want to run our model and produce outputs
#' - here, every day for 100 years
times = seq(1, 100 * 365, 1)

model_output2 = deSolve::lsoda(state, times, SIS2, parameters)

#' Get model values in equilibrium
state_equilibrium = model_output2[nrow(model_output2),-1]

### START ANALYSIS HERE ####


#' Run with vaccination from equilibrium
parameters[["vacc_cov_f"]] = 0.50
parameters[["vacc_cov_m"]] = 0.45
parameters[["vacc_eff"]] = 0.89


model_output2 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)

model_output2 = tibble::as_tibble(model_output2)
model_output2 = model_output2 %>%
  pivot_longer(cols = c("Sf1l", "Sf1h", "Sf2l", "Sf2h", "Sf3l", "Sf3h", "Sm1l", "Sm1h", "Sm2l", "Sm2h", "Sm3l", "Sm3h", 
                        "If1l", "If1h", "If2l", "If2h", "If3l", "If3h", "Im1l", 'Im1h', 'Im2l', "Im2h", "Im3l", "Im3h",
                        "Vf1l", "Vf1h", "Vf2l", "Vf2h", "Vf3l", "Vf3h", "Vm1l", "Vm1h", "Vm2l", "Vm2h", "Vm3l", "Vm3h", 
                        "Cf1l", "Cf1h", "Cf2l", "Cf2h", "Cf3l", "Cf3h", "Cm1l", "Cm1h", "Cm2l", "Cm2h", "Cm3l", "Cm3h"), names_to = "model_compartment") %>%
  # create a new column with only the compartment (without the age group)
  mutate(compartment = gsub("\\d+", "", model_compartment)) %>%
  # create a new column with only the age group (without the compartment)
  mutate(age_group = gsub("\\D+", "", model_compartment) %>% as.numeric()) %>%
  mutate(state = substr(compartment, 1, 1), gender = substr(compartment, 2, 2),
         activity = substr(compartment, 3, 3)) 

#date of first timestep
start_date = as.Date("2016-01-01")

#convert model time to corresponding dates
model_output2 = model_output2 %>%
  mutate(date = start_date + time - 1)

#filter model output based on dates
model_output2 %>% # 2030
filter(state == "I", date %in% c(as.Date("2030-01-01"),as.Date("2030-12-31")), gender == "f")

model_output2 %>% # 2100
filter(state == "I", date %in% c(as.Date("2100-01-01"),as.Date("2100-12-31")), gender == "f")







# temporarily used to get vaccine proportions (total and by-age)

HPV_2030 = model_output2 %>% # for 2030
  filter(state == "I", date %in% c(as.Date("2021-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

HPV_2030_prev = (HPV_2030/sum(N)*0.5) %>% print



############################################################

# All strategies at 0.89% efficacy
#V_1 = model_output2 # 0.30
#V_2 = model_output2 # 0.40
#V_3 = model_output2 # 0.50
#V_4 = model_output2 # 0.60
#V_5 = model_output2 # 0.70
#V_6 = model_output2 # 0.80
#V_7 = model_output2 # 0.90
#V_8 = model_output2 # 1
#saveRDS(V_1, "V_1.RDS")
#saveRDS(V_2, "V_2.RDS")
#saveRDS(V_3, "V_3.RDS")
#saveRDS(V_4, "V_4.RDS")
#saveRDS(V_5, "V_5.RDS")
#saveRDS(V_6, "V_6.RDS")
#saveRDS(V_7, "V_7.RDS")
#saveRDS(V_8, "V_8.RDS")

#############################
###### Model Output #########
#############################

# Mid-point prevalence by age and HPV incidence by age (relative to strata size)
########
# total numbers - by age
age1_2030 = model_output2 %>% # for 2030
  filter(age_group =="1", state == "I", date %in% c(as.Date("2030-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

age1_2100 = model_output2 %>% # for 2100
  filter(age_group =="1", state == "I", date %in% c(as.Date("2100-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

age2_2030 = model_output2 %>% # for 2030
  filter(age_group =="2", state == "I", date %in% c(as.Date("2030-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

age2_2100 = model_output2 %>% # for 2100
  filter(age_group =="2", state == "I", date %in% c(as.Date("2100-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

age3_2030 = model_output2 %>% # for 2030
  filter(age_group =="3", state == "I", date %in% c(as.Date("2030-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

age3_2100 = model_output2 %>% # for 2100
  filter(age_group =="3", state == "I", date %in% c(as.Date("2100-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

# Total numbers by age-group
print(c(age1_2030, age2_2030, age3_2030)) %>% sum # probably not interested in these
print(c(age1_2100, age2_2100, age3_2100)) %>% sum


# prevalence per 100,000 by age-group -> First table in results
p_100000_30 = c(age1_2030, age2_2030, age3_2030) 
print((p_100000_30)/(N*0.5)*100000)

#print(age1_2030/(N[1]*0.5))*100000
#print(age2_2030/(N[2]*0.5))*100000
#print(age3_2030/(N[3]*0.5))*100000

print(sum(age1_2030, age2_2030, age3_2030)/(sum(N)*0.5)*100000)


p_100000_100 = c(age1_2100, age2_2100, age3_2100)
print((p_100000_100)/(N*0.5)*100000)
#print(age1_2100/(N[1]*0.5))*100000
#print(age2_2100/(N[2]*0.5))*100000
#print(age3_2100/(N[3]*0.5))*100000
print(sum(age1_2100, age2_2100, age3_2100)/(sum(N)*0.5)*100000)

# Cumulative incidence per 100,000 by age-group

# incidence - 15-24
age1_incidence_bytime = model_output2 %>% 
  filter(age_group == "1", state == "C", gender == "f") 

age1_cases_2030 = age1_incidence_bytime %>% # for 2030
  filter(date %in% c(as.Date("2030-01-01"):as.Date("2030-12-31"))) %>% 
  .[["value"]] %>% sum

age1_cases_2030_per100000 = (age1_cases_2030/(N[1]*0.50))*100000

age1_cases_2100 = age1_incidence_bytime %>% # for 2100
  filter(date %in% c(as.Date("2100-01-01"):as.Date("2100-12-31"))) %>% 
  .[["value"]] %>% sum

age1_cases_2100_per100000 = (age1_cases_2100/(N[1]*0.50))*100000

# incidence - 25-34
age2_incidence_bytime = model_output2 %>% 
  filter(age_group == "2", state == "C", gender == "f") 

age2_cases_2030 = age2_incidence_bytime %>% # for 2030
  filter(date %in% c(as.Date("2030-01-01"):as.Date("2030-12-31"))) %>% 
  .[["value"]] %>% sum

age2_cases_2030_per100000 = (age2_cases_2030/(N[2]*0.50))*100000

age2_cases_2100 = age2_incidence_bytime %>% # for 2100
  filter(date %in% c(as.Date("2100-01-01"):as.Date("2100-12-31"))) %>% 
  .[["value"]] %>% sum

age2_cases_2100_per100000 = (age2_cases_2100/(N[2]*0.50))*100000

# incidence - 35+
age3_incidence_bytime = model_output2 %>% 
  filter(age_group == "3", state == "C", gender == "f") 

age3_cases_2030 = age3_incidence_bytime %>% # for 2030
  filter(date %in% c(as.Date("2030-01-01"):as.Date("2030-12-31"))) %>% 
  .[["value"]] %>% sum

age3_cases_2030_per100000 = (age3_cases_2030/(N[3]*0.50))*100000

age3_cases_2100 = age3_incidence_bytime %>% # for 2100
  filter(date %in% c(as.Date("2100-01-01"):as.Date("2100-12-31"))) %>% 
  .[["value"]] %>% sum

age3_cases_2100_per100000 = (age3_cases_2100/(N[3]*0.50))*100000

print(c(age1_cases_2030_per100000, age1_cases_2100_per100000))
print(c(age2_cases_2030_per100000, age2_cases_2100_per100000))
print(c(age3_cases_2030_per100000, age3_cases_2100_per100000))


##########################################################


# Combined values - incidence -> per 100000 (C compartment) for all ages -> main table in results

prop_cancer = (0.172+0.136)/2
cancer_delay = 1/(-log(1 - 0.99)/10) 

# Midpoint point-prevalence per 100,000 by year

HPV_2030 = model_output2 %>% # for 2030
  filter(state == "I", date %in% c(as.Date("2030-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

HPV_2030_prev = (HPV_2030/119603)*100000

HPV_2100 = model_output2 %>% # for 2100
  filter(state == "I", date %in% c(as.Date("2100-07-01")), gender == "f") %>% 
  .[["value"]] %>% sum

HPV_2100_prev = (HPV_2100/119603)*100000

# Cumulative incidence per 100,000 by year
# HPV incidence

HPV_incidence_bytime = model_output2 %>% 
  filter(state == "C", gender == "f") 

HPV_cases_2030 = HPV_incidence_bytime %>% # for 2030
  filter(date %in% c(as.Date("2030-01-01"):as.Date("2030-12-31"))) %>% 
  .[["value"]] %>% sum

HPV_cases_2030_per100000 = (HPV_cases_2030/119603)*100000

HPV_cases_2100 = HPV_incidence_bytime %>% # for 2100
  filter(date %in% c(as.Date("2100-01-01"):as.Date("2100-12-31"))) %>% 
  .[["value"]] %>% sum

HPV_cases_2100_per100000 = (HPV_cases_2100/119603)*100000

# CANCER

cancer_incidence_bytime = model_output2 %>% 
  filter(state == "C", gender == "f") %>% 
  mutate(cancer_incidence = (value*prop_cancer), cancer_time = date + (cancer_delay*365))

model_output2 %>% # 2030
filter(state == "C", date %in% c(as.Date("2030-01-01"),as.Date("2030-12-31")), gender == "f") %>% 
  mutate(cancer_incidence = (value*prop_cancer), cancer_time = date + (cancer_delay*365))

cancer_cases_2030 = cancer_incidence_bytime %>%
  filter(cancer_time >= as.Date("2030-01-01") & cancer_time < as.Date("2031-01-01")) %>% 
.[["cancer_incidence"]] %>% sum

cancer_cases_2030_per100000 = (cancer_cases_2030/119603)*100000

model_output2 %>% # 2100
  filter(state == "C", date %in% c(as.Date("2100-01-01"),as.Date("2100-12-31")), gender == "f") %>% 
  mutate(cancer_incidence = (value*prop_cancer), cancer_time = date + (cancer_delay*365))

cancer_cases_2100 = cancer_incidence_bytime %>%
  filter(cancer_time >= as.Date("2100-01-01") & cancer_time < as.Date("2101-01-01")) %>% 
  .[["cancer_incidence"]] %>% sum

cancer_cases_2100_per100000 = (cancer_cases_2100/119603)*100000

#print(c(HPV_2030_prev, HPV_2100_prev, 
        #HPV_cases_2030_per100000, HPV_cases_2100_per100000, 
        #cancer_cases_2030_per100000, cancer_cases_2100_per100000))

#print(c(HPV_2030_prev, HPV_2100_prev))
#print(c(HPV_cases_2030_per100000, HPV_cases_2100_per100000))
#print(c(cancer_cases_2030_per100000, cancer_cases_2100_per100000))
 

print(c(HPV_2030_prev))
print(c(HPV_2100_prev))
print(c(HPV_cases_2030_per100000))
print(c(HPV_cases_2100_per100000))
print(c(cancer_cases_2030_per100000))
print(c(cancer_cases_2100_per100000))


#print(c(HPV_2030_prev, HPV_cases_2030_per100000,cancer_cases_2030_per100000))
      
#print(c(HPV_2100_prev, HPV_cases_2100_per100000,cancer_cases_2100_per100000))

#############################

#####################################
########## Checks and plots #########
#####################################

#' One useful sanity check to see if there are any errors in your model is to see if population size remains constant
model_output2 %>%
  filter(state != "C") %>% 
  group_by(time) %>% #First summarize your data by summing the value column in each age group
  summarise(N = sum(value)) %>%
  tail() #The population size should remain equal to popsize i.e. 239206

#to check equilibrium estimates
#########

prev = model_output2 %>%
  filter(time == max(time), gender == "f") %>%
  group_by(age_group, activity)%>%
  mutate(N = sum(value)) %>%
  filter(state %in% c("I")) %>%
  summarise(value = sum(value), N=mean(N)) %>%
  #' Rename modelled_prev to modelled
  mutate(modelled = value/N) %>%
  mutate(age_group = recode(age_group,
                            "1" = "15-24",
                            "2" = "25-34",
                            "3" = "35+")) %>%
  mutate(partner_year = recode(activity,
                               "l" = "Low",
                               "h" = "High"))%>%
  select(age_group, partner_year, modelled) %>%
  left_join(prevalence_data, by=c("age_group", "partner_year"))

#########

#' Show the fit of the model to the data (red dots are modelled prevalence, black dots are data with 95% CI)  
model_output2 %>%
  filter(time == max(time), gender == "f") %>%
  group_by(age_group, activity) %>%
  mutate(N = sum(value)) %>%
  filter(state %in% c("I")) %>%
  summarise(value = sum(value), N=mean(N)) %>%
  #' Rename modelled_prev to modelled
  mutate(modelled = value/N) %>%
  mutate(age_group = recode(age_group,
                            "1" = "15-24",
                            "2" = "25-34",
                            "3" = "35+")) %>%
  mutate(partner_year = recode(activity,
                               "l" = "Low",
                               "h" = "High")) %>%
  select(age_group, partner_year, modelled) %>%
  left_join(prevalence_data, by=c("age_group", "partner_year")) %>%
  #' Rename est to observed
  mutate(observed = est) %>%
  #' Pivot longer based on observed and modelled
  pivot_longer(cols = c("observed", "modelled"), names_to = "variable") %>%
  ggplot(aes(x=age_group))+
  #' Specify aesthetics based on the variable (observed or modelled), this will create a legend
  #' - ymin and ymax only exist for the observed values, not for the modelled values
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5)+
  #' Specify aesthetics based on the variable (observed or modelled), this will create a legend
  geom_point(aes(y = value, shape = variable, colour = variable, fill = variable))+
  facet_grid(.~partner_year)+
  labs(x="Age-group",
       y="HPV-16/18 prevalence",
       colours = fill)+
  #' Manually map shapes to values in 'variable' (used for the shape aesthetic)
  scale_shape_manual(values = c(observed = 22, modelled = 23))+
  #' Manually map colours to values in 'variable' (used for the colour aesthetic)
  scale_colour_manual(values = c(observed = "#000000", modelled = "#ff0000"))+
  #' Manually map colours to values in 'variable' (used for the fill aesthetic)
  scale_fill_manual(values = c(observed = "#000000", modelled = "#ff0000"))+
  theme_bw()


# Proportions HPV prevalence in women over time, stratified by age and activity

#####################
# Version 2 (3 plots)

inf_colors = c("Low-activity" = "blue", "High-activity" = "red")

ggplot()+
  geom_line(data = model_output2 %>% 
              filter(state == "I", gender == "f", activity =="l", age_group =="1") %>% 
              mutate(age_group = recode(age_group, "1" = "15-24 years", "2" = "25-34 years", "3" = "35+ years")), 
            mapping = aes(x=(time/365)+2016, y=value/(N[1]*0.5*(1-prop_high_risk[1])), colour = "Low-activity"), size=0.5)+
  geom_line(data = model_output2 %>% 
              filter(state == "I", gender == "f", activity =="h", age_group =="1") %>% 
              mutate(age_group = recode(age_group, "1" = "15-24 years", "2" = "25-34 years","3" = "35+ years")), 
            mapping = aes(x=(time/365)+2016, y=value/(N[1]*0.5*(prop_high_risk[1])), colour = "High-activity"), size=0.5)+
  geom_line(data = model_output2 %>% 
              filter(state == "I", gender == "f", activity =="l", age_group =="2") %>% 
              mutate(age_group = recode(age_group, "1" = "15-24 years", "2" = "25-34 years", "3" = "35+ years")), 
            mapping = aes(x=(time/365)+2016, y=value/(N[2]*0.5*(1-prop_high_risk[2])), colour = "Low-activity"), size=0.5)+
  geom_line(data = model_output2 %>% 
              filter(state == "I", gender == "f", activity =="h", age_group =="2") %>% 
              mutate(age_group = recode(age_group, "1" = "15-24 years",  "2" = "25-34 years", "3" = "35+ years")), 
            mapping = aes(x=(time/365)+2016, y=value/(N[2]*0.5*(prop_high_risk[2])), colour = "High-activity"), size=0.5) +
  geom_line(data = model_output2 %>% 
              filter(state == "I", gender == "f", activity =="l", age_group =="3") %>% 
              mutate(age_group = recode(age_group, "1" = "15-24 years", "2" = "25-34 years", "3" = "35+ years")), 
            mapping = aes(x=(time/365)+2016, y=value/(N[3]*0.5*(1-prop_high_risk[3])), colour = "Low-activity"), size=0.5)+
  geom_line(data = model_output2 %>% 
              filter(state == "I", gender == "f", activity =="h", age_group =="3") %>% 
              mutate(age_group = recode(age_group, "1" = "15-24 years", "2" = "25-34 years", "3" = "35+ years")), 
            mapping = aes(x=(time/365)+2016, y=value/(N[3]*0.5*(prop_high_risk[3])), colour = "High-activity"), size=0.5) +
  labs(x = "Year",
       y = "Proportion of Effectively Vaccinated Women",
       colour = "Sexual-activity") +
  facet_grid(.~age_group)+
  scale_color_manual(values = inf_colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  scale_y_continuous(labels=scales::percent)+
  theme_bw()+ 
  theme(text = element_text(size = 12))+
  theme(axis.text = element_text(size = 8))

#####################
# Mid-point prevalence by age and activity -> used with figure 4 (check 2016 [2016-01-01], 2030 [2030-07-01], 2100[2100-07-01].)
prop_prev = model_output2 %>% # 15-24 years
  filter(age_group =="1", state == "I", date %in% c(as.Date("2100-07-01")), gender == "f") %>%print() # change date for analysis
# low activity
print(c(prop_prev$value[1]/(N[1]*0.5*(1-prop_high_risk[1]))))
# high activity
print(c(prop_prev$value[2]/(N[1]*0.5*(prop_high_risk[1]))))


prop_prev = model_output2 %>%  # 25-34 years
  filter(age_group =="2", state == "I", date %in% c(as.Date("2100-07-01")), gender == "f") %>%print() # change date for analysis
# low activity
print(c(prop_prev$value[1]/(N[2]*0.5*(1-prop_high_risk[2]))))
# high activity
print(c(prop_prev$value[2]/(N[2]*0.5*(prop_high_risk[2]))))


prop_prev = model_output2 %>%  # 35+ years
  filter(age_group =="3", state == "I", date %in% c(as.Date("2100-07-01")), gender == "f") %>%print() # change date for analysis
# low activity
print(c(prop_prev$value[1]/(N[3]*0.5*(1-prop_high_risk[3]))))
# high activity
print(c(prop_prev$value[2]/(N[3]*0.5*(prop_high_risk[3]))))
#####################################

############################################## 
######### Sensitivity Analysis ###############
############################################## 

#' Create a (named) list with parameters options
scenarios = list(
  "Current" = c(
    "vacc_cov_f" = 0.30, "vacc_cov_m" = 0.27, "vacc_eff" = 0.89),
  "40% girls, 36% boys" = c(
    "vacc_cov_f" = 0.40, "vacc_cov_m" = 0.36, "vacc_eff" = 0.89),
  "50% girls, 45% boys" = c(
    "vacc_cov_f" = 0.50, "vacc_cov_m" = 0.45, "vacc_eff" = 0.89),
  "60% girls, 54% boys" = c(
    "vacc_cov_f" = 0.60, "vacc_cov_m" = 0.54, "vacc_eff" = 0.89),
  "70% girls, 63% boys" = c(
    "vacc_cov_f" = 0.70, "vacc_cov_m" = 0.63, "vacc_eff" = 0.89),
  "80% girls, 72% boys" = c(
    "vacc_cov_f" = 0.80, "vacc_cov_m" = 0.72, "vacc_eff" = 0.89),
  "90% girls, 81% boys" = c(
    "vacc_cov_f" = 0.90, "vacc_cov_m" = 0.81, "vacc_eff" = 0.89),
  "100% girls, 90% boys" = c(
    "vacc_cov_f" = 0.100, "vacc_cov_m" = 0.90, "vacc_eff" = 0.89),
  "30% girls, no boys" = c(
    "vacc_cov_f" = 0.30, "vacc_cov_m" = 0.0, "vacc_eff" = 0.89))

#' Create an empty list for the output of each analysis
out = list()

#' Loop through each value of the list
for(i in 1:length(scenarios)){
  #' Copy main parameter list for this analysis
  parameters_sens = parameters
  
  #' Overwrite with values for this analysis
  parameters_sens[["vacc_cov_f"]] = scenarios[[i]]["vacc_cov_f"]
  parameters_sens[["vacc_cov_m"]] = scenarios[[i]]["vacc_cov_m"]
  parameters_sens[["vacc_eff"]] = scenarios[[i]]["vacc_eff"]
  
  #' Run the model with the variables for this analysis
  model_sens = deSolve::lsoda(state_equilibrium, times, SIS2, parameters_sens)
  
  #' Convert data.frame to a tibble
  model_sens = tibble::as_tibble(model_sens)
  
  model_sens = model_sens %>%
    pivot_longer(cols = c("Sf1l", "Sf1h", "Sf2l", "Sf2h", "Sf3l", "Sf3h", "Sm1l", "Sm1h", "Sm2l", "Sm2h", "Sm3l", "Sm3h",
                          "If1l", "If1h", "If2l", "If2h", "If3l", "If3h", "Im1l", 'Im1h', 'Im2l', "Im2h", "Im3l", "Im3h",
                          "Vf1l", "Vf1h", "Vf2l", "Vf2h", "Vf3l", "Vf3h", "Vm1l", "Vm1h", "Vm2l", "Vm2h", "Vm3l", "Vm3h",
                          "Cf1l", "Cf1h", "Cf2l", "Cf2h", "Cf3l", "Cf3h", "Cm1l", "Cm1h", "Cm2l", "Cm2h", "Cm3l", "Cm3h"),
                 names_to = "model_compartment") %>%
    #this changes the numbers from data-type 'deSolve' to data-type 'double'
    mutate(value = as.numeric(value), time = as.numeric(time))
  
  #' Add column with scenario name
  model_sens = model_sens %>% mutate(scen = names(scenarios)[i])
  
  #' Append to output list
  out[[length(out)+1]] = model_sens
}

#' Combine output in a single tibble
out = out %>% dplyr::bind_rows()

#' Process output
out = out %>%
  # create a new column with only the compartment (without the age group)
  mutate(compartment = gsub("\\d+", "", model_compartment)) %>%
  # create a new column with only the age group (without the compartment)
  mutate(age_group = gsub("\\D+", "", model_compartment) %>% as.numeric()) %>%
  mutate(state = substr(compartment, 1, 1), gender = substr(compartment, 2, 2), activity = substr(compartment, 3, 3))

#' Convert timesteps to date
start_date = as.Date("2016-01-01")
out %>% mutate(date = start_date + time - 1)

#' combine ages and risk groups if you wish
out_minimised = out %>% 
  mutate(date = start_date + time - 1) %>% 
  group_by(scen, date, state, gender) %>%
  summarise(value = sum(value))

out2 = out %>% 
  group_by(scen, time, state, gender) %>%
  summarise(value = sum(value))

saveRDS(out, "out.RDS")
saveRDS(out_minimised, "out_minimised.RDS")


vac_colors = c("Current" = "blue", "40% girls, 36% boys" = "red", "50% girls, 45% boys" = "purple",
           "60% girls, 54% boys" = "brown", "70% girls, 63% boys" = "orange", "80% girls, 72% boys" = "gold",
           "90% girls, 81% boys" = "green", "100% girls, 90% boys" = "pink")

ggplot()+ 
  geom_line(data = out_minimised %>% 
              filter(scen == "30% girls, 27% boys", state == "I", gender == "f"), 
            mapping = aes(x=date, y=value/(sum(N)*0.5)*100000))+
  labs(x = "Year",
       y = "Proportion of infected women") +
  theme_bw()




#################################
## Extra Datasets manipulation ##
#################################

parameters[["vacc_cov_f"]] = 0.30
parameters[["vacc_cov_m"]] = 0.27
parameters[["vacc_eff"]] = 0.89
model_curr_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_curr_0.89 = tibble::as_tibble(model_curr_0.89)
model_curr_0.89 = model_curr_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)
  
saveRDS(model_curr_0.89, "model_curr_0.89.RDS")

parameters[["vacc_cov_f"]] = 0.30
parameters[["vacc_cov_m"]] = 0.27
parameters[["vacc_eff"]] = 1
model_curr_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_curr_1 = tibble::as_tibble(model_curr_1)
model_curr_1 = model_curr_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_curr_1, "model_curr_1.RDS")

parameters[["vacc_cov_f"]] = 0.40
parameters[["vacc_cov_m"]] = 0.36
parameters[["vacc_eff"]] = 0.89
model_40_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_40_0.89 = tibble::as_tibble(model_40_0.89)
model_40_0.89 = model_40_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_40_0.89, "model_40_0.89.RDS")

parameters[["vacc_cov_f"]] = 0.40
parameters[["vacc_cov_m"]] = 0.36
parameters[["vacc_eff"]] = 1
model_40_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_40_1 = tibble::as_tibble(model_40_1)
model_40_1 = model_40_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_40_1, "model_40_1.RDS")

parameters[["vacc_cov_f"]] = 0.50
parameters[["vacc_cov_m"]] = 0.45
parameters[["vacc_eff"]] = 0.89
model_50_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_50_0.89 = tibble::as_tibble(model_50_0.89)
model_50_0.89 = model_50_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_50_0.89, "model_50_0.89.RDS")

parameters[["vacc_cov_f"]] = 0.50
parameters[["vacc_cov_m"]] = 0.45
parameters[["vacc_eff"]] = 1
model_50_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_50_1 = tibble::as_tibble(model_50_1)
model_50_1 = model_50_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_50_1, "model_50_1.RDS")

parameters[["vacc_cov_f"]] = 0.60
parameters[["vacc_cov_m"]] = 0.54
parameters[["vacc_eff"]] = 0.89
model_60_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_60_0.89 = tibble::as_tibble(model_60_0.89)
model_60_0.89 = model_60_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_60_0.89, "model_60_0.89.RDS")

parameters[["vacc_cov_f"]] = 0.60
parameters[["vacc_cov_m"]] = 0.54
parameters[["vacc_eff"]] = 1
model_60_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_60_1 = tibble::as_tibble(model_60_1)
model_60_1 = model_60_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_60_1, "model_60_1.RDS")

parameters[["vacc_cov_f"]] = 0.70
parameters[["vacc_cov_m"]] = 0.63
parameters[["vacc_eff"]] = 0.89
model_70_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_70_0.89 = tibble::as_tibble(model_70_0.89)
model_70_0.89 = model_70_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_70_0.89, "model_70_0.89.RDS")

parameters[["vacc_cov_f"]] = 0.70
parameters[["vacc_cov_m"]] = 0.63
parameters[["vacc_eff"]] = 1
model_70_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_70_1 = tibble::as_tibble(model_70_1)
model_70_1 = model_70_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_70_1, "model_70_1.RDS")

parameters[["vacc_cov_f"]] = 0.80
parameters[["vacc_cov_m"]] = 0.72
parameters[["vacc_eff"]] = 0.89
model_80_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_80_0.89 = tibble::as_tibble(model_80_0.89)
model_80_0.89 = model_80_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_80_0.89, "model_80_0.89.RDS")

parameters[["vacc_cov_f"]] = 0.80
parameters[["vacc_cov_m"]] = 0.72
parameters[["vacc_eff"]] = 1
model_80_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_80_1 = tibble::as_tibble(model_80_1)
model_80_1 = model_80_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_80_1, "model_80_1.RDS")

parameters[["vacc_cov_f"]] = 0.90
parameters[["vacc_cov_m"]] = 0.81
parameters[["vacc_eff"]] = 0.89
model_90_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_90_0.89 = tibble::as_tibble(model_90_0.89)
model_90_0.89 = model_90_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_90_0.89, "model_90_0.89.RDS")

parameters[["vacc_cov_f"]] = 0.90
parameters[["vacc_cov_m"]] = 0.81
parameters[["vacc_eff"]] = 1
model_90_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_90_1 = tibble::as_tibble(model_90_1)
model_90_1 = model_90_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_90_1, "model_90_1.RDS")

parameters[["vacc_cov_f"]] = 1
parameters[["vacc_cov_m"]] = 0.90
parameters[["vacc_eff"]] = 0.89
model_100_0.89 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_100_0.89 = tibble::as_tibble(model_100_0.89)
model_100_0.89 = model_100_0.89 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_100_0.89, "model_100_0.89.RDS")

parameters[["vacc_cov_f"]] = 1
parameters[["vacc_cov_m"]] = 0.90
parameters[["vacc_eff"]] = 1
model_100_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_100_1 = tibble::as_tibble(model_100_1)
model_100_1 = model_100_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)
 
saveRDS(model_100_1, "model_100_1.RDS")

parameters[["vacc_cov_f"]] = 0
parameters[["vacc_cov_m"]] = 0
parameters[["vacc_eff"]] = 0.89
model_100_1 = deSolve::lsoda(state_equilibrium, times, SIS2, parameters)
model_100_1 = tibble::as_tibble(model_100_1)
model_100_1 = model_100_1 %>% 
  mutate(Sf = Sf1l+Sf1h+Sf2l+Sf2h+Sf3l+Sf3h) %>% 
  mutate(Sm = Sm1l+Sm1h+Sm2l+Sm2h+Sm3l+Sm3h) %>% 
  mutate(If = If1l+If1h+If2l+If2h+If3l+If3h) %>% 
  mutate(Im = Im1l+Im1h+Im2l+Im2h+Im3l+Im3h) %>% 
  mutate(Vf = Vf1l+Vf1h+Vf2l+Vf2h+Vf3l+Vf3h) %>% 
  mutate(Vm = Vm1l+Vm1h+Vm2l+Vm2h+Vm3l+Vm3h) %>% 
  mutate(Cf = Cf1l+Cf1h+Cf2l+Cf2h+Cf3l+Cf3h) %>% 
  mutate(Vf1 = Vf1l+Vf1h) %>% 
  mutate(Vf2 = Vf2l+Vf2h) %>% 
  mutate(Vf3 = Vf3l+Vf3h)

saveRDS(model_100_1, "model_0.RDS")
#################################

############################
########## Plots ###########
############################

model_curr_0.89 = readRDS("model_curr_0.89.RDS")
model_curr_1 = readRDS("model_curr_1.RDS")
model_40_0.89 = readRDS("model_40_0.89.RDS")
model_40_1 = readRDS("model_40_1.RDS")
model_50_0.89 = readRDS("model_50_0.89.RDS")
model_50_1 = readRDS("model_50_1.RDS")
model_60_0.89 = readRDS("model_60_0.89.RDS")
model_60_1 = readRDS("model_60_1.RDS")
model_70_0.89 = readRDS("model_70_0.89.RDS")
model_70_1 = readRDS("model_70_1.RDS")
model_80_0.89 = readRDS("model_80_0.89.RDS")
model_80_1 = readRDS("model_80_1.RDS")
model_90_0.89 = readRDS("model_90_0.89.RDS")
model_90_1 = readRDS("model_90_1.RDS")
model_100_0.89 = readRDS("model_100_0.89.RDS")
model_100_1 = readRDS("model_100_1.RDS")

colors = c("Current" = "blue", "40% girls, 36% boys" = "red", "50% girls, 45% boys" = "purple",
           "60% girls, 54% boys" = "brown", "70% girls, 63% boys" = "orange", "80% girls, 72% boys" = "gold",
           "90% girls, 81% boys" = "green", "100% girls, 90% boys" = "pink")

# Number of susceptibles after 0.89% efficacy
ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "Current"), size=1.5) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "40% girls, 36% boys"), size=1.5) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "50% girls, 45% boys"), size=1.5) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "60% girls, 54% boys"), size=1.5) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "70% girls, 63% boys"), size=1.5) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "80% girls, 72% boys"), size=1.5) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "90% girls, 81% boys"), size=1.5) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=Sf, colour = "100% girls, 90% boys"), size=1.5) +
  labs(x = "Year",
       y = "Number of Susceptible Females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  theme_bw()

# Number of susceptibles after 100% efficacy
ggplot()+
  geom_line(data = model_curr_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "Current"), size=1.5) +
  geom_line(data = model_40_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "40% girls, 36% boys"), size=1.5) +
  geom_line(data = model_50_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "50% girls, 45% boys"), size=1.5) +
  geom_line(data = model_60_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "60% girls, 54% boys"), size=1.5) +
  geom_line(data = model_70_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "70% girls, 63% boys"), size=1.5) +
  geom_line(data = model_80_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "80% girls, 72% boys"), size=1.5) +
  geom_line(data = model_90_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "90% girls, 81% boys"), size=1.5) +
  geom_line(data = model_100_1, mapping = aes(x=(time/365)+2016, y=Sf, colour = "100% girls, 90% boys"), size=1.5) +
  labs(x = "Year",
       y = "Number of Susceptible Females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  theme_bw()


# Number of infectious after 0.89% efficacy
ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "Current"), size=1.5) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "40% girls, 36% boys"), size=1.5) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "50% girls, 45% boys"), size=1.5) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "60% girls, 54% boys"), size=1.5) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "70% girls, 63% boys"), size=1.5) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "80% girls, 72% boys"), size=1.5) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "90% girls, 81% boys"), size=1.5) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=If, colour = "100% girls, 90% boys"), size=1.5) +
  labs(x = "Year",
       y = "Number of Infectious Females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  theme_bw()

#proportion of infectious after 0.89% efficacy -> I currently use this one

colors = c("Current" = "blue", "40% girls, 36% boys" = "red", "50% girls, 45% boys" = "purple",
           "60% girls, 54% boys" = "brown", "70% girls, 63% boys" = "orange", "80% girls, 72% boys" = "gold",
           "90% girls, 81% boys" = "green", "100% girls, 90% boys" = "pink")

ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "Current"), size=1) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "40% girls, 36% boys"), size=1, alpha=0.7) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "50% girls, 45% boys"), size=1, alpha=0.7) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "60% girls, 54% boys"), size=1, alpha=0.7) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "70% girls, 63% boys"), size=1, alpha=0.7) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "80% girls, 72% boys"), size=1, alpha=0.7) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "90% girls, 81% boys"), size=1, alpha=0.7) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "100% girls, 90% boys"), size=1, alpha=0.7) +
  labs(x = "Year",
       y = "Proportion of Infectious Women",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw()+
  theme(text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 8))

#proportion of infectious in sensitivity analysis compared -> i currently use this one
colors = c("Current" = "blue", "40% girls, 36% boys" = "red", "50% girls, 45% boys" = "purple",
           "60% girls, 54% boys" = "brown", "70% girls, 63% boys" = "orange")

lines = c("89% Efficacy" = 1, "100% Efficacy" = 9)

ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "Current", linetype = "89% Efficacy"), size=1) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "40% girls, 36% boys", linetype = "89% Efficacy"), size=1, alpha = 0.5) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "50% girls, 45% boys", linetype = "89% Efficacy"), size=1, alpha = 0.5) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "70% girls, 63% boys", linetype = "89% Efficacy"), size=1, alpha = 0.5) + 
  geom_line(data = model_curr_1, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "Current", linetype = "100% Efficacy"), size=0.6) +
  geom_line(data = model_40_1, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "40% girls, 36% boys", linetype = "100% Efficacy"), size=0.6) +
  geom_line(data = model_50_1, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "50% girls, 45% boys", linetype = "100% Efficacy"), size=0.6) +
  geom_line(data = model_70_1, mapping = aes(x=(time/365)+2016, y=If/119603, colour = "70% girls, 63% boys", linetype = "100% Efficacy"), size=0.6) + 
  labs(x = "Year",
       y = "Proportion of Infectious Women",
       colour = "Vaccine Strategy",
       linetype = "Sensitivity Analysis") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = lines)+
  coord_cartesian(xlim = c(2020, 2100)) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw()+
  theme(text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 8))

# Number of Vaccinated after 89% efficacy -> currently use this one
ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "Current"), size=1.5) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "40% girls, 36% boys"), size=1.5) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "50% girls, 45% boys"), size=1.5) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "60% girls, 54% boys"), size=1.5) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "70% girls, 63% boys"), size=1.5) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "80% girls, 72% boys"), size=1.5) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "90% girls, 81% boys"), size=1.5) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=Vf/(sum(N)*0.5), colour = "100% girls, 90% boys"), size=1.5) +
  labs(x = "Year",
       y = "Proportion of Vaccinated Women aged >15 years",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw()+
  theme(text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 8))

############################################################

# Number of Vaccinated after 0.89% efficacy -> i use this one
colors = c("Current" = "blue", "40% girls, 36% boys" = "red", "50% girls, 45% boys" = "purple",
           "60% girls, 54% boys" = "brown", "70% girls, 63% boys" = "orange", "80% girls, 72% boys" = "gold",
           "90% girls, 81% boys" = "green", "100% girls, 90% boys" = "pink")
lines = c("15-24 years" = 6, "25-34 years" = 3, "35+ years" = 1)

v1 = ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "Current"), size=1) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "40% girls, 36% boys"), size=1) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "50% girls, 45% boys"), size=1) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "60% girls, 54% boys"), size=1) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "70% girls, 63% boys"), size=1) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "80% girls, 72% boys"), size=1) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "90% girls, 81% boys"), size=1) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=Vf1/(N[1]*0.5), colour = "100% girls, 90% boys"), size=1) +
  labs(x = "Year",
       y = "Proportion of Effectively Vaccinated Females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  #scale_linetype_manual(values = lines)+
  coord_cartesian(xlim = c(2020, 2100), ylim = c(0.0, 0.5)) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw()+
  theme(text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 8))


v2 = ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "Current"), size=1) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "40% girls, 36% boys"), size=1) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "50% girls, 45% boys"), size=1) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "60% girls, 54% boys"), size=1) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "70% girls, 63% boys"), size=1) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "80% girls, 72% boys"), size=1) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "90% girls, 81% boys"), size=1) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=Vf2/(N[2]*0.5), colour = "100% girls, 90% boys"), size=1) +
  labs(x = "Year",
       y = "Proportion of Effectively Vaccinated Females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  #scale_linetype_manual(values = lines)+
  coord_cartesian(xlim = c(2020, 2100), ylim = c(0.0, 0.5)) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw()+
  theme(text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 8))

v3 = ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "Current"), size=1) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "40% girls, 36% boys"), size=1) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "50% girls, 45% boys"), size=1) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "60% girls, 54% boys"), size=1) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "70% girls, 63% boys"), size=1) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "80% girls, 72% boys"), size=1) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "90% girls, 81% boys"), size=1) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=Vf3/(N[3]*0.5), colour = "100% girls, 90% boys"), size=1) +
  labs(x = "Year",
       y = "Proportion of Effectively Vaccinated Females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  #scale_linetype_manual(values = lines)+
  coord_cartesian(xlim = c(2020, 2100), ylim = c(0.0, 0.5)) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw()+
  theme(text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 8))



figure = ggarrange(plotlist = v1, v2, v3,
          labels(c("15-24 years", "25-34 years", "35+ years")),
          ncol = 3, crow = 1)

cowplot::plot_grid(plotlist = v1, v2, v3,
                   labels(c("15-24 years", "25-34 years", "35+ years")),
                   ncol = 3, crow = 1)

# Number of Vaccinated after 100% efficacy
ggplot()+
  geom_line(data = model_curr_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "Current"), size=1.5) +
  geom_line(data = model_40_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "40% girls, 36% boys"), size=1.5) +
  geom_line(data = model_50_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "50% girls, 45% boys"), size=1.5) +
  geom_line(data = model_60_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "60% girls, 54% boys"), size=1.5) +
  geom_line(data = model_70_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "70% girls, 63% boys"), size=1.5) +
  geom_line(data = model_80_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "80% girls, 72% boys"), size=1.5) +
  geom_line(data = model_90_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "90% girls, 81% boys"), size=1.5) +
  geom_line(data = model_100_1, mapping = aes(x=(time/365)+2016, y=Vf, colour = "100% girls, 90% boys"), size=1.5) +
  labs(x = "Year",
       y = "Number of Vaccinated Females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  theme_bw()

# Number of Cancer after 89%% efficacy
ggplot()+
  geom_line(data = model_curr_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "Current"), size=1.5) +
  geom_line(data = model_40_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "40% girls, 36% boys"), size=1.5) +
  geom_line(data = model_50_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "50% girls, 45% boys"), size=1.5) +
  geom_line(data = model_60_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "60% girls, 54% boys"), size=1.5) +
  geom_line(data = model_70_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "70% girls, 63% boys"), size=1.5) +
  geom_line(data = model_80_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "80% girls, 72% boys"), size=1.5) +
  geom_line(data = model_90_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "90% girls, 81% boys"), size=1.5) +
  geom_line(data = model_100_0.89, mapping = aes(x=(time/365)+2016, y=Cf, colour = "100% girls, 90% boys"), size=1.5) +
  labs(x = "Year",
       y = "Number of Cancerous females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  theme_bw()


# Number of Cancer after 100% efficacy
ggplot()+
  geom_line(data = model_curr_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "Current"), size=1.5) +
  geom_line(data = model_40_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "40% girls, 36% boys"), size=1.5) +
  geom_line(data = model_50_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "50% girls, 45% boys"), size=1.5) +
  geom_line(data = model_60_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "60% girls, 54% boys"), size=1.5) +
  geom_line(data = model_70_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "70% girls, 63% boys"), size=1.5) +
  geom_line(data = model_80_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "80% girls, 72% boys"), size=1.5) +
  geom_line(data = model_90_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "90% girls, 81% boys"), size=1.5) +
  geom_line(data = model_100_1, mapping = aes(x=(time/365)+2016, y=Cf, colour = "100% girls, 90% boys"), size=1.5) +
  labs(x = "Year",
       y = "Number of Cancerous females",
       colour = "Vaccine Strategy") +
  scale_color_manual(values = colors) +
  coord_cartesian(xlim = c(2020, 2100)) +
  theme_bw()


# read in all strategy version at 0.89 efficacy
V_1 = readRDS("V_1.RDS")
V_2 = readRDS("V_2.RDS")
V_3 = readRDS("V_3.RDS")
V_4 = readRDS("V_4.RDS")
V_5 = readRDS("V_5.RDS")
V_6 = readRDS("V_6.RDS")
V_7 = readRDS("V_7.RDS")
V_8 = readRDS("V_8.RDS")
