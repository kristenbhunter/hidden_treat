# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : simulation.r
#
# Programmer Names   : Kristen Hunter, kristenbhunter@gmail.com
#                      Katy McKeough
#                      Nicole Pashley
#
# Last Updated       : May 2020
#
# Purpose            : Runs the numerical simulations in "hidden treatment" paper
#
# Input              : Takes in the base directory where github repo is located
# Output             : Saves out plots in a subdirectory called "plots"
#
# References         :
#
# Platform           : R
# Version            : v3.6.3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(plyr)
library(reshape2)
library(wesanderson)

### change to the directory of the github repo on your own computer
base.dir = '/Users/khunter/Dropbox/hidden_treat/'

source(paste(base.dir, 'simulation_functions.R', sep = ''))
plot.dir = paste(base.dir, 'plots/', sep = '')

set.seed(040218)

########################################
# simulation parameters
########################################

# experiment size
r = 25
# simulation reps
sim.reps = 1000
# effects
theta.a = 2
theta.b = 2
theta.ab = 2
real.theta = list(theta.a = theta.a, theta.b = theta.b, theta.ab = theta.ab)

#### generate rho matrix
rho.00.10 = rho.00.01 = rho.00.11 =
  rho.10.01 = rho.10.11 = rho.01.11 = 0.4

rho.matrix = diag(4)
rho.matrix[1,2] = rho.matrix[2,1] = rho.00.10
rho.matrix[1,3] = rho.matrix[3,1] = rho.00.01
rho.matrix[1,4] = rho.matrix[4,1] = rho.00.11
rho.matrix[2,3] = rho.matrix[3,2] = rho.10.01
rho.matrix[2,4] = rho.matrix[4,2] = rho.10.11
rho.matrix[3,4] = rho.matrix[4,3] = rho.01.11

########################################
# run sim for a particular piB
########################################

run.sim.piB = function(sim.reps, N, piB,
                       Y.list,
                       szz.constants,
                       real.theta = NULL,
                       title = NULL,
                       diagnostic.plot = FALSE,
                       diagnostic.piB = c(0.05, 0.1, 0.2))
{
  # store results
  output = NULL

  # loop through simulation reps
  for(j in 1:sim.reps)
  {
    # treatment assignment
    WA = sample(c(rep(1, N/2), rep(0, N/2)))

    # for case 1, we'll take any assignment,
    # but for case 2 we reject assignments with N < 2
    WB.1 = rbinom(N, 1, piB)
    WB.2 = rep(0, N)

    while( (sum(  WA    & WB.2) < 2) |
           (sum(  WA    & (1-WB.2)) < 2) |
           (sum( (1-WA) & WB.2) < 2) |
           (sum( (1-WA) & (1-WB.2) ) < 2) )
    {
      WB.2 = rbinom(N, 1, piB)
    }

    # number of units for each treatment assignment

    # case I
    N1dot = sum(WA)
    N0dot = sum(1 - WA)

    # case II
    N00.2 = sum( (1-WA) * (1-WB.2) )
    N01.2 = sum( (1-WA) * (WB.2)   )
    N10.2 = sum( (WA)   * (1-WB.2) )
    N11.2 = sum( (WA)   * (WB.2)   )
    N.2.list = list('N00' = N00.2, 'N01' = N01.2, 'N10' = N10.2, 'N11' = N11.2)

    # estimated pib and nz for case 2
    piB.hat = (N01.2 + N11.2)/N
    nz00.2 = (1 - piB.hat) * N0dot
    nz01.2 = piB.hat       * N0dot
    nz10.2 = (1 - piB.hat) * N1dot
    nz11.2 = piB.hat       * N1dot

    # calculate observed Y
    Yobs.1 =
      Y.list[['Y00']] * (1 - WA) * (1 - WB.1) +
      Y.list[['Y10']] * (WA)     * (1 - WB.1) +
      Y.list[['Y01']] * (1 - WA) * (WB.1) +
      Y.list[['Y11']] * (WA)     * (WB.1)

    Yobs.2 =
      Y.list[['Y00']] * (1 - WA) * (1 - WB.2) +
      Y.list[['Y10']] * (WA)     * (1 - WB.2) +
      Y.list[['Y01']] * (1 - WA) * (WB.2) +
      Y.list[['Y11']] * (WA)     * (WB.2)

    # true thetas
    theta.a = 1/2*(Y.list[['Y11']] - Y.list[['Y01']] +
                   Y.list[['Y10']] - Y.list[['Y00']])
    theta.a.fs = mean(theta.a)

    # Case I theta estimator
    y1.obs.1 = sum( Yobs.1 * (WA)   )/N1dot
    y1.obs.0 = sum( Yobs.1 * (1-WA) )/N0dot
    theta.a.1 = y1.obs.1 - y1.obs.0

    # Case II theta estimator
    y2.obs.1 = (sum(Yobs.2 * (WA)   * (WB.2)   )/N11.2 +
                  sum(Yobs.2 * (WA)   * (1-WB.2) )/N10.2)/2
    y2.obs.0 = (sum(Yobs.2 * (1-WA) * (WB.2)   )/N01.2 +
                  sum(Yobs.2 * (1-WA) * (1-WB.2) )/N00.2)/2
    theta.a.2 = y2.obs.1 - y2.obs.0

    # Case I expected value
    exp.theta.1 =
      (1 - piB) * mean(Y.list[['Y10']] - Y.list[['Y00']]) +
      (piB)     * mean(Y.list[['Y11']] - Y.list[['Y01']])

    # Case I variance
    var.hat.theta.1 =
      calc.s2z(Yobs.1, WA, za = 1)/N1dot +
      calc.s2z(Yobs.1, WA, za = 0)/N0dot

    # Case II expected value
    exp.theta.2 = theta.a.fs

    # Case II variances
    var.hat.1.theta.2 = (1/4)*(
      calc.s2z(Yobs.2, WA, WB.2, za = 1, zb = 1)/N11.2 +
      calc.s2z(Yobs.2, WA, WB.2, za = 1, zb = 0)/N10.2 +
      calc.s2z(Yobs.2, WA, WB.2, za = 0, zb = 1)/N01.2 +
      calc.s2z(Yobs.2, WA, WB.2, za = 0, zb = 0)/N00.2
    )

    # Case II variance version 2 (found in supplementary material)
    var.hat.2.theta.2 = (1/4)*(
      calc.expnzinv(za = 1, zb = 1, WA, piB.hat)*calc.s2z(Yobs.2, WA, WB.2, za = 1, zb = 1) +
      calc.expnzinv(za = 1, zb = 0, WA, piB.hat)*calc.s2z(Yobs.2, WA, WB.2, za = 1, zb = 0) +
      calc.expnzinv(za = 0, zb = 1, WA, piB.hat)*calc.s2z(Yobs.2, WA, WB.2, za = 0, zb = 1) +
      calc.expnzinv(za = 0, zb = 0, WA, piB.hat)*calc.s2z(Yobs.2, WA, WB.2, za = 0, zb = 0)
    )

    # true variance calculations
    var.theta.1 = calc.var.theta.1(Y.list, N1dot, N0dot, piB, szz.constants)
    var.theta.2 = calc.var.theta.2(Y.list, N.2.list, piB)

    # save out the estimates, variances, and true value
    output = rbind(output, data.frame(
      sim.id = rep(j, 2),
      type = c('Case I','Case II'),
      estimate = c(theta.a.1, theta.a.2),
      exp.estimate = c(exp.theta.1, exp.theta.2),
      variance = c(var.hat.theta.1, var.hat.1.theta.2),
      true.variance = c(var.theta.1, var.theta.2),
      true.value = rep(theta.a.fs, 2)
    ))

  } # end sim.reps loop

  # save out measures of performance
  output$piB = piB
  output$ci.lower = output$estimate - qnorm(0.975) * sqrt(output$variance)
  output$ci.upper = output$estimate + qnorm(0.975) * sqrt(output$variance)
  output$cover =
    (output$true.value >= output$ci.lower) &
    (output$true.value <= output$ci.upper)
  output$width = output$ci.upper - output$ci.lower

  # summaries of performance
  output.summary = ddply(output, c('type'), summarise,
                         piB = unique(piB),
                         cover = mean(cover, na.rm = TRUE),
                         width = mean(width, na.rm = TRUE),
                         rel.bias = mean((exp.estimate - true.value)/true.value, na.rm = TRUE),
                         var.rel.bias = mean((variance - true.variance)/true.variance, na.rm = TRUE),
                         mse = mean((exp.estimate - true.value)^2 + true.variance, na.rm = TRUE)
  )

  return(output.summary)
}

########################################
# main simulation function
########################################

sim.function = function(r, y.model,
                        real.theta = NULL,
                        rho.matrix = NULL,
                        title = NULL,
                        diagnostic.plot = FALSE,
                        diagnostic.piB = c(0.05, 0.1, 0.2)
)
{
  # number of units
  N = 4 * r
  # for computation speed: calculates values used
  # in later equations that are constant given a certain N
  szz.constants = fancyS2zz.constants.fun(N)

  # sequence of probabilities to try
  piB = seq(0.05, 0.95, 0.05)

  # generate potential outcomes
  Y.all = generate.y(N, y.model = y.model,
                     real.theta = real.theta,
                     rho.matrix = rho.matrix)
  Y.list = Y.all[['Y']]
  real.theta = Y.all[['real.theta']]

  # check:
  # finite sample values of main and interaction effects in generated data
  # should match the inputs
  theta.a.fs = 1/2*(mean(Y.list[['Y11']]) - mean(Y.list[['Y01']]) +
                      mean(Y.list[['Y10']]) - mean(Y.list[['Y00']]))
  theta.b.fs = 1/2*(mean(Y.list[['Y11']]) - mean(Y.list[['Y10']]) +
                      mean(Y.list[['Y01']]) - mean(Y.list[['Y00']]))
  theta.ab.fs = 1/2*(mean(Y.list[['Y11']]) - mean(Y.list[['Y01']]) -
                       mean(Y.list[['Y10']]) + mean(Y.list[['Y00']]))
  if(round(real.theta[['theta.a']], 2)  != round(theta.a.fs, 2) |
     round(real.theta[['theta.b']], 2)  != round(theta.b.fs, 2) |
     round(real.theta[['theta.ab']], 2) != round(theta.ab.fs, 2))
  {
    print(theta.a.fs)
    print(theta.b.fs)
    print(theta.ab.fs)
    stop('Error in potential outcome assignment')
  }

  # save out results
  results = NULL

  # loop through piB
  print(paste('piB iteration out of', length(piB)))
  for(i in 1:length(piB))
  {
    print(i)
    results = rbind(results, run.sim.piB(sim.reps, N, piB[i],
                                         Y.list = Y.list,
                                         real.theta = real.theta,
                                         szz.constants,
                                         title = title,
                                         diagnostic.plot,
                                         diagnostic.piB))
  }
  return(results)
}

########################################
# strict additivity, zero interaction
########################################

title = 'Strict additivity with no interaction'
real.theta[['theta.ab']] = 0
results.strict.no = sim.function(r, y.model = 'strict',
                                 real.theta = real.theta, title = title)
results.strict.no$model = title

########################################
# strict additivity, interaction
########################################

title = 'Strict additivity with interaction'
real.theta[['theta.ab']] = theta.ab
results.strict.int = sim.function(r, y.model = 'strict',
                                  real.theta = real.theta, title = title)
results.strict.int$model = title

########################################
# positively correlated outcomes
########################################

title = 'Positively correlated outcomes'
results.corr = sim.function(r, y.model = 'correlated', real.theta = real.theta,
                            rho.matrix = rho.matrix, title = title)
results.corr$model = title

########################################
# plot all results on one plot
########################################

all.results = rbind(results.strict.no, results.strict.int, results.corr)
results.plot = plot.performance(all.results)

plot.file = paste(plot.dir, 'performance.png', sep = '')
png(plot.file, width = 1200, height = 800)
results.plot
dev.off()

########################################
# detailed plots for supplement
########################################

png.width = 1000
png.height = 600

png(paste(plot.dir, 'performance_strict_nointeraction.png', sep = ''), width = png.width, height = png.height)
plot.detailed.performance(results.strict.no)
dev.off()

png(paste(plot.dir, 'performance_strict_interaction.png', sep = ''), width = png.width, height = png.height)
plot.detailed.performance(results.strict.int)
dev.off()

png(paste(plot.dir, 'performance_correlated.png', sep = ''), width = png.width, height = png.height)
plot.detailed.performance(results.corr)
dev.off()

