#' sinh arcshinh Probability Density Function
#' @param x: vector of quantiles.
#' @param delta: shape parameter.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return density of the sinh-arcsinh distribution with shape parameter delta
#' @export
dsas <- function (x, delta, log = FALSE)
{
  logPDF <- dnorm(sinh(delta * asinh(x)), log = T) + log(delta) -
    0.5 * log(1 + x^2) + log(cosh(delta * asinh(x)))
  ifelse(log, return(logPDF), return(exp(logPDF)))
}

#' twopiece sinh arcshinh Probability Density Function
#' @param x: vector of quantiles.
#' @param delta: shape parameter.
#' @param par1: scale parameter 1
#' @param par2: scale parameter 2 or skewness parameter
#' @param param: parameterisation. Defaults to "tp". Information about the "eps" and "isf" parameterizations can be found in the References.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return density of the twopiece sinh-arcsinh distribution
#' @export
dtpsas <- function (x, mu, par1, par2, delta, param = "tp", log = FALSE)
{
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0 & delta > 0, logPDF <- log(2) +
             ifelse(x < mu, dsas((x - mu)/par1, delta, log = T),
                    dsas((x - mu)/par2, delta, log = T)) - log(par1 +
                                                                 par2), logPDF <- "invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1 & delta > 0, logPDF <- ifelse(x <
                                                                      mu, dsas((x - mu)/(sigma * (1 + gamma)), delta, log = T),
                                                                    dsas((x - mu)/(sigma * (1 - gamma)), delta, log = T)) -
             log(sigma), logPDF <- "invalid arguments: par1 or/and delta is/are no positive or/and abs(par2) is no less that 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0 & delta > 0, logPDF <- log(2) +
             ifelse(x < mu, dsas((x - mu)/(sigma * gamma), delta,
                                 log = T), dsas((x - mu)/(sigma/gamma), delta,
                                                log = T)) - log(sigma * (gamma + 1/gamma)), logPDF <- "invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization isf")
  }
  ifelse(is.numeric(logPDF), ifelse(log, return(logPDF), return(exp(logPDF))),
         logPDF)
}


#' sinh-arcsinh Cumulative Distribution Function
#' @param x: vector of quantiles.
#' @param delta: shape parameter.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @return cumulative distribution of the twopiece sinh-arcsinh distribution
#' @export
psas <- function (x, delta, log.p = FALSE)
{
  logCDF <- pnorm(sinh(delta * asinh(x)), log.p = T)
  ifelse(log.p, return(logCDF), return(exp(logCDF)))
}

#' twopiece sinh-arcsinh Cumulative Distribution Function
#' @param x: vector of quantiles.
#' @param delta: shape parameter.
#' @param par1: scale parameter 1
#' @param par2: scale parameter 2 or skewness parameter
#' @param param: parameterisation. Defaults to "tp". Information about the "eps" and "isf" parameterizations can be found in the References.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @return cumulative distribution of the twopiece sinh-arcsinh distribution with shape parameter delta
#' @export
ptpsas <- function (x, mu, par1, par2, delta, param = "tp", log.p = FALSE)
{
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0, CDF <- ifelse(x < mu, 2 *
                                                par1 * psas((x - mu)/par1, delta, log.p = F)/(par1 +
                                                                                                par2), (par1 + par2 * (2 * psas((x - mu)/par2, delta,
                                                                                                                                log.p = F) - 1))/(par1 + par2)), CDF <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1, CDF <- ifelse(x <
                                                       mu, (1 + gamma) * psas((x - mu)/(sigma * (1 + gamma)),
                                                                              delta, log.p = F), gamma + (1 - gamma) * psas((x -
                                                                                                                               mu)/(sigma * (1 - gamma)), delta, log.p = F)), CDF <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0, CDF <- ifelse(x < mu, 2 *
                                                  gamma^2 * psas((x - mu)/(sigma * gamma), delta, log.p = F)/(1 +
                                                                                                                gamma^2), (gamma^2 - 1 + 2 * psas((x - mu)/(sigma/gamma),
                                                                                                                                                  delta, log.p = F))/(1 + gamma^2)), CDF <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  ifelse(is.numeric(CDF), ifelse(log.p, return(log(CDF)), return(CDF)),
         CDF)
}


#' sinh-arcsinh Quantile Function
#' @param p: vector of probabilities.
#' @param delta: shape parameter.
#' @return quantile function of the sinh-arcsinh distribution with shape parameter delta
#' @export
qsas <- function (p, delta)
{
  Q <- sinh(asinh(qnorm(p))/delta)
  return(Q)
}

#' twopiece sinh-arcsinh Quantile Function
#' @param p: vector of probabilities
#' @param delta: shape parameter.
#' @param par1: scale parameter 1
#' @param par2: scale parameter 2 or skewness parameter
#' @param param: parameterisation. Defaults to "tp". Information about the "eps" and "isf" parameterizations can be found in the References.
#' @return quantile function of the twopiece sinh-arcsinh distribution
#' @export
qtpsas <- function (p, mu, par1, par2, delta, param = "tp")
{
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0, Q <- ifelse(p < par1/(par1 +
                                                        par2), mu + par1 * qsas(0.5 * p * (par1 + par2)/par1,
                                                                                delta), mu + par2 * qsas(0.5 * ((par1 + par2) * (1 +
                                                                                                                                   p) - 2 * par1)/par2, delta)), Q <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1, Q <- ifelse(p < 0.5 *
                                                     (1 + gamma), mu + sigma * (1 + gamma) * qsas(p/(1 +
                                                                                                       gamma), delta), mu + sigma * (1 - gamma) * qsas((p -
                                                                                                                                                          gamma)/(1 - gamma), delta)), Q <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0, Q <- ifelse(p < gamma^2/(1 +
                                                             gamma^2), mu + sigma * gamma * qsas(0.5 * p * (1 +
                                                                                                              gamma^2)/gamma^2, delta), mu + sigma * qsas(0.5 *
                                                                                                                                                            (p * (1 + gamma^2) + 1 - gamma^2), delta)/gamma),
           Q <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(Q)
}


# sinh-arcsinh Random Number Generation
#' @param n: number of observations.
#' @param delta: shape parameter.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return random number generation for the sinh-arcsinh distribution with shape parameter delta
#' @export
rsas <- function (n, delta, log = FALSE)
{
  sample <- sinh(asinh(rnorm(n, 0, 1))/delta)
  return(sample)
}

# twopiece sinh-arcsinh Random Number Generation
#' @param n: number of observations.
#' @param delta: shape parameter.
#' @param par1: scale parameter 1
#' @param par2: scale parameter 2 or skewness parameter
#' @param param: parameterisation. Defaults to "tp". Information about the "eps" and "isf" parameterisations can be found in the References.
#' @return random number generation for the twopiece sinh-arcsinh distribution
#' @export
rtpsas <- function (n, mu, par1, par2, delta, param = "tp")
{
  if (param == "tp") {
    ifelse(par1 > 0 & par2 > 0 & delta > 0, sample <- ifelse(runif(n) <
                                                               par1/(par1 + par2), mu - par1 * abs(rsas(n, delta)),
                                                             mu + par2 * abs(rsas(n, delta))), sample <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & abs(gamma) < 1 & delta > 0, sample <- ifelse(runif(n) <
                                                                      0.5 * (1 + gamma), mu - sigma * (1 + gamma) * abs(rsas(n,
                                                                                                                             delta)), mu + sigma * (1 - gamma) * abs(rsas(n, delta))),
           sample <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma = par1
    gamma = par2
    ifelse(sigma > 0 & gamma > 0 & delta > 0, sample <- ifelse(runif(n) <
                                                                 gamma^2/(1 + gamma^2), mu - sigma * gamma * abs(rsas(n,
                                                                                                                      delta)), mu + sigma * abs(rsas(n, delta))/gamma),
           sample <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(sample)
}
