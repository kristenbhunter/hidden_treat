# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : simulation_functions.r
#
# Programmer Names   : Kristen Hunter, kristenbhunter@gmail.com
#                      Katy McKeough
#                      Nicole Pashley
#
# Last Updated       : May 2020
#
# Purpose            : Supplementary functions for simulations
#
# Input              : None
# Output             : None
#
# Platform           : R
# Version            : v3.6.3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


########################################
# generate Y(bz)
# Section 2.1
########################################

Yz = function(Y.list, za, zb)
{
  Yz =
    (1 - za) * (1 - zb) * Y.list$Y00 +
    za  * (1 - zb) * Y.list$Y10 +
    (1 - za) *      zb  * Y.list$Y01 +
    za  *      zb  * Y.list$Y11
  return(Yz)
}

########################################
# Calculate covariance of W indicators
# Lemma C.1
########################################

cov.W = function(za1, zb1, za2, zb2, N1dot, N0dot, piB, same.i = TRUE)
{

  # useful expressions
  N = N1dot + N0dot
  EWB1 = piB * zb1 + (1 - piB) * (1 - zb1)
  EWB2 = piB * zb2 + (1 - piB) * (1 - zb2)
  cov.WA.diff = (N1dot * N0dot)/(N^2 * (N - 1))
  cov.WA.same = (N1dot * N0dot)/(N^2)

  if(!same.i & za1 != za2)
    cov = EWB1 * EWB2 * cov.WA.diff
  else if (!same.i & za1 == za2)
    cov = -EWB1 * EWB2 * cov.WA.diff
  else if(same.i & (za1 != za2 | zb1 != zb2))
    cov = -EWB1 * EWB2 * cov.WA.same
  else if(same.i & za1 == za2 & zb1 == zb2)
  {
    Nzadot = za1 * N1dot + (1 - za1) * N0dot
    cov = (EWB1 * Nzadot * (N - EWB1 * Nzadot))/N^2
  }

  return(cov)
}

########################################
# Calculate simplified covariance of potential outcomes
# Eqn 11
########################################

calc.S2zz = function(za1, zb1, za2, zb2, Y.list, N)
{
  Yz1 = Yz(Y.list, za1, zb1)
  Yz2 = Yz(Y.list, za2, zb2)

  S2zz = sum( (Yz1 - mean(Yz1)) * (Yz2 - mean(Yz2)) )/(N-1)

  return(S2zz)
}

########################################
# Calculate s^2(z)
# Section 3.1, equation before Equation 17
########################################

calc.s2z = function(Yobs, WA, WB = NULL, za, zb = NULL)
{
  # accommodate either S^2(z) or S^2(z_A)
  if(is.null(WB))
  {
    Wz = (WA == za)
  } else
  {
    # Wza
    Wz = (WA == za) & (WB == zb)
  }

  s2z = var(Yobs[Wz])

  return(s2z)
}

########################################
# Calculate Exp(1/Nz)
# Lemma A.1
########################################

calc.expnzinv = function(za, zb, WA, piB.hat)
{
  Nzadot = ifelse(za == 1, sum(WA), sum(1 - WA))
  exp.WB = ifelse(zb == 1, piB.hat, 1 - piB.hat)

  constant = 1/( 1 - piB.hat^Nzadot - (1 - piB.hat)^Nzadot)
  calc.term = function(n){
    expr = (1/n) * choose(Nzadot, n) * exp.WB^n * (1- exp.WB)^(Nzadot - n)
    return(expr)
  }
  n.vec = seq(1, Nzadot - 1)
  sum.terms = sapply(n.vec, calc.term)
  return(constant * sum(sum.terms))
}


########################################
# Calculate covariance of potential outcomes
# Section D.1
########################################

calc.fancy.S2zz = function(za1, zb1, za2, zb2, Y.list, N1dot, N0dot, piB, szz.constants)
{
  Yz1 = Yz(Y.list, za1, zb1)
  Yz2 = Yz(Y.list, za2, zb2)

  cov.samei = cov.W(za1, zb1, za2, zb2, N1dot, N0dot, piB, same.i = TRUE)
  cov.diffi = cov.W(za1, zb1, za2, zb2, N1dot, N0dot, piB, same.i = FALSE)

  szz.val1 = sum(cov.samei * Yz1 * Yz2)
  szz.val2 = sum(cov.diffi * Yz1[szz.constants[['combos']][,1]] * Yz2[szz.constants[['combos']][,2]])

  szz = szz.val1 + szz.val2

  if(za1 == za2 & zb1 == zb2)
    return(szz)
  else
    return(-szz)

}

########################################
# calculate constants for fancy S2zz above
#
# This is a time-saving function.
# It moves out parts of the function that don't need to be repeated
# so we can move it out of for loops and speed things up
########################################

fancyS2zz.constants.fun = function(N)
{
  combinations = unique(t(combn(rep(seq(1,N),2),2)))
  combinations = combinations[combinations[,1] != combinations[,2],]
  return(list('combos' = combinations))
}

########################################
# Calculate true variance of theta.1
# Eqn 14
# alternative version at end of section D.1
########################################

calc.var.theta.1 = function(Y.list, N1dot, N0dot, piB, szz.constants)
{
  var.1 =
    calc.fancy.S2zz(1, 1, 1, 1, Y.list, N1dot, N0dot, piB, szz.constants)/N1dot^2 +
    calc.fancy.S2zz(1, 0, 1, 0, Y.list, N1dot, N0dot, piB, szz.constants)/N1dot^2 +
    calc.fancy.S2zz(0, 1, 0, 1, Y.list, N1dot, N0dot, piB, szz.constants)/N0dot^2 +
    calc.fancy.S2zz(0, 0, 0, 0, Y.list, N1dot, N0dot, piB, szz.constants)/N0dot^2

  var.2 = 2*(
    calc.fancy.S2zz(1, 1, 1, 0, Y.list, N1dot, N0dot, piB, szz.constants)/N1dot^2 +
    calc.fancy.S2zz(0, 1, 0, 0, Y.list, N1dot, N0dot, piB, szz.constants)/N0dot^2)

  var.3 = 2*(
    calc.fancy.S2zz(1, 1, 0, 1, Y.list, N1dot, N0dot, piB, szz.constants)/(N1dot*N0dot) +
    calc.fancy.S2zz(1, 1, 0, 0, Y.list, N1dot, N0dot, piB, szz.constants)/(N1dot*N0dot) +
    calc.fancy.S2zz(1, 0, 0, 1, Y.list, N1dot, N0dot, piB, szz.constants)/(N1dot*N0dot) +
    calc.fancy.S2zz(1, 0, 0, 0, Y.list, N1dot, N0dot, piB, szz.constants)/(N1dot*N0dot))

  return(var.1 - var.2 + var.3)
}

########################################
# Calculate true variance of theta.2
# Eqn 18
########################################

calc.var.theta.2 = function(Y.list, N.list, piB)
{
  # N.list = N.2.list
  N = N.list$N00 + N.list$N01 + N.list$N10 + N.list$N11

  var.1 =
    calc.S2zz(1, 1, 1, 1, Y.list, N) * (N - N.list$N11)/N.list$N11 +
    calc.S2zz(1, 0, 1, 0, Y.list, N) * (N - N.list$N10)/N.list$N10 +
    calc.S2zz(0, 1, 0, 1, Y.list, N) * (N - N.list$N01)/N.list$N01 +
    calc.S2zz(0, 0, 0, 0, Y.list, N) * (N - N.list$N00)/N.list$N00

  var.2 = 2*(
    calc.S2zz(1, 1, 1, 0, Y.list, N) +
    calc.S2zz(0, 1, 0, 0, Y.list, N))

  var.3 = 2*(
    calc.S2zz(1, 1, 0, 1, Y.list, N) +
    calc.S2zz(1, 1, 0, 0, Y.list, N) +
    calc.S2zz(1, 0, 0, 1, Y.list, N) +
    calc.S2zz(1, 0, 0, 0, Y.list, N))

  return( (1/(4 * N)) * (var.1 - var.2 + var.3))
}

########################################
# generate potential outcomes
# under different simulation models
########################################

generate.y = function(N,
                      y.model,
                      real.theta = NULL,
                      sigma = 1,
                      rho.matrix = NULL)
{
  if(y.model == 'strict')
  {
    theta.0 = rnorm(N, mean = 0, sd = sigma)
    Y00 = theta.0 +
      real.theta[['theta.ab']]
    Y10 = theta.0 +
      real.theta[['theta.a']]
    Y01 = theta.0 +
      real.theta[['theta.b']]
    Y11 = theta.0 +
      real.theta[['theta.a']] +
      real.theta[['theta.b']] +
      real.theta[['theta.ab']]
  } else
  {
    means = c(real.theta[['theta.ab']],
              real.theta[['theta.a']],
              real.theta[['theta.b']],
              real.theta[['theta.a']] +
              real.theta[['theta.b']] +
              real.theta[['theta.ab']])
    Y = rmvnorm(N, mean = means, sigma = sigma * rho.matrix)
    Y00 = Y[,1]
    Y10 = Y[,2]
    Y01 = Y[,3]
    Y11 = Y[,4]

    real.theta = list()
    real.theta[['theta.a']] = mean(Y10 + Y11 - Y01 - Y00)/2
    real.theta[['theta.b']] = mean(Y01 + Y11 - Y10 - Y00)/2
    real.theta[['theta.ab']] = mean(Y11 - Y10 - Y01 + Y00)/2
  }

  return(list(Y = list(Y00 = Y00, Y10 = Y10, Y01 = Y01, Y11 = Y11),
              real.theta = real.theta))
}

########################################
# plot simulation performance
# Figure 1
########################################

plot.performance = function(all.results)
{
  # pick plot colors
  plot.colors = c(wes_palette('Darjeeling1')[2], wes_palette('Darjeeling1')[4])

  # rearrange plot data to easily put into ggplot facets
  plot.data = all.results[,c('type', 'piB', 'cover', 'width', 'model')]
  colnames(plot.data) = c('Case', 'piB', 'Coverage', 'Interval Width', 'Model')
  plot.data = melt(plot.data, c('Case', 'piB', 'Model'))

  # plot 95% coverage on coverage plots
  coverage.lines = data.frame(
    Model = rep(unique(plot.data$Model), 2),
    variable = c(rep('Coverage', 3), rep('Interval Width', 3)),
    nominal = c(rep(0.95, 3), rep(NA, 3)))

  # change names to be in correct order for the figure
  plot.data$Model = factor(
    plot.data$Model,
    levels = c('Strict additivity with no interaction', 'Strict additivity with interaction', 'Positively correlated outcomes'))

  # make names more legible
  model.names = c(
    'Strict additivity with no interaction' = 'Strict additivity\nwith no interaction',
    'Strict additivity with interaction' = 'Strict additivity\nwith interaction',
    'Positively correlated outcomes' = 'Positively correlated\noutcomes',
    'Coverage' = 'Coverage', 'Interval Width' = 'Interval Width'
  )

  # final plot
  results.plot =
    ggplot(plot.data, aes(x = piB, y = value, color = Case)) +
    facet_grid(variable ~ Model, scales = 'free_y', labeller = as_labeller(model.names)) +
    xlab(expression(pi[B])) +
    geom_line(size = 1.5) +
    ylab('') +
    geom_hline(data = coverage.lines, aes(yintercept = nominal)) +
    scale_colour_manual(name = "Case", values = plot.colors) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    theme(text = element_text(size = 30)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 30), axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      legend.title = element_text(size = 20), legend.text = element_text(size = 20)
    )

  return(results.plot)
}

########################################
# plot detailed simulation performance
# supplementary materials graphs
# Figures 1-3 in supplement
########################################

#############
# generate shared legend
#############

g.legend = function(a.gplot) {
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

plot.detailed.performance = function(results, title)
{
  plot.colors = c(wes_palette('Darjeeling1')[2], wes_palette('Darjeeling1')[4])

  plots = list()

  # coverage
  plots[[1]] = ggplot(results, aes(x = piB, y = cover, color = factor(type))) +
    geom_path(size = 1.5) +
    geom_hline(yintercept = 0.95) +
    ggtitle('Coverage') +
    xlab('') +
    ylab('') +
    scale_colour_manual(name = "Case", values = plot.colors) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ylim(0, 1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 18), legend.text = element_text(size = 15),
      plot.title = element_text(size = 25)
    )

  # interval width
  plots[[2]] = ggplot(results, aes(x = piB, y = width, color = factor(type))) +
    geom_path(size = 1.5) +
    ggtitle('Average interval width') +
    xlab('') +
    ylab('') +
    scale_colour_manual(name = "Case", values = plot.colors) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ylim(0, 2) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 18), legend.text = element_text(size = 15),
      plot.title = element_text(size = 25)
    )

  # average relative bias
  plots[[3]] = ggplot(results, aes(x = piB, y = rel.bias, color = factor(type))) +
    geom_path(size = 1.5) +
    ggtitle('Average relative bias') +
    xlab('') +
    ylab('') +
    scale_colour_manual(name = "Case", values = plot.colors) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    geom_hline(yintercept = 0) +
    ylim(-1,1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 18), legend.text = element_text(size = 15),
      plot.title = element_text(size = 25)
    )

  # MSE
  plots[[4]] = ggplot(results, aes(x = piB, y = mse, color = factor(type))) +
    geom_path(size = 1.5) +
    ggtitle('MSE') +
    xlab('') +
    ylab('') +
    scale_colour_manual(name = "Case", values = plot.colors) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ylim(0,4) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 18), legend.text = element_text(size = 15),
      plot.title = element_text(size = 25)
    )

  # Bias of variance estimate
  plots[[5]] = ggplot(results, aes(x = piB, y = var.rel.bias, color = factor(type))) +
    geom_path(size = 1.5) +
    ggtitle('Average relative bias\nof variance estimator') +
    xlab('') +
    ylab('') +
    scale_colour_manual(name = "Case", values = plot.colors) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ylim(-1,1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 18), legend.text = element_text(size = 15),
      plot.title = element_text(size = 25)
    )

  shared.legend = g.legend(plots[[1]])

  plot.object = grid.arrange(
    arrangeGrob(
      plots[[1]] + theme(legend.position = 'none'),
      plots[[2]] + theme(legend.position = 'none'),
      plots[[3]] + theme(legend.position = 'none'),
      plots[[4]] + theme(legend.position = 'none'),
      plots[[5]] + theme(legend.position = 'none'),
      shared.legend,
      ncol = 3
    ),
    bottom = textGrob(expression(pi[B]), gp = gpar(fontsize = 25))
  )

  return(plot.object)

}
