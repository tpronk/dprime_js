// Quantile of of X given a Gaussian distribution with mean = 0 and SD = 1
// Adapted from https://github.com/errcw/gaussian/blob/master/lib/gaussian.js
var ppf = function(x) {
  return -1 * Math.sqrt(2) * ierfc(2 * x);
};

var erfc = function(x) {
  var z = Math.abs(x);
  var t = 1 / (1 + z / 2);
  var r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 +
          t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 +
          t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
          t * (-0.82215223 + t * 0.17087277)))))))))
  return x >= 0 ? r : 2 - r;
};

var ierfc = function(x) {
  if (x >= 2) { return -100; }
  if (x <= 0) { return 100; }

  var xx = (x < 1) ? x : 2 - x;
  var t = Math.sqrt(-2 * Math.log(xx / 2));

  var r = -0.70711 * ((2.30753 + t * 0.27061) /
          (1 + t * (0.99229 + t * 0.04481)) - t);

  for (var j = 0; j < 2; j++) {
    var err = erfc(r) - xx;
    r += err / (1.12837916709551257 * Math.exp(-(r * r)) - r * err);
  }

  return (x < 1) ? r : -r;
};

// D-prime using loglinear approach to correct for extreme values. Adapted from:
// https://cran.r-project.org/web/packages/splithalfr/vignettes/gng_dprime.html
// n_hit  - Number of hits
// n_miss - Number of Misses
// n_cr   - Number of correct rejections
// n_fa   - Number of false alarms
dprime = function(n_hit, n_miss, n_cr, n_fa) {
  var p_hit = (n_hit + 0.5) / ((n_hit + 0.5) + n_miss + 1);
  var p_fa = (n_fa + 0.5) / ((n_fa + 0.5) + n_cr + 1);
  return (ppf(p_hit) - ppf(p_fa))
}

/*
# Cases below tested against this R-function
dprime <- function(n_hit, n_miss, n_cr, n_fa) {
  p_hit <- (n_hit + 0.5) / ((n_hit + 0.5) + n_miss + 1)
  p_fa <- (n_fa + 0.5) / ((n_fa + 0.5) + n_cr + 1)  
  return (qnorm(p_hit) - qnorm(p_fa))
}
dprime(10, 10, 10, 10);
# R: 0, JS: 0
dprime(20, 15, 5, 3);
# R: 0.4911764, JS: 0.4911764
dprime(100, 0, 0, 0);
# R: 2.762656, JS: 2.762656
*/