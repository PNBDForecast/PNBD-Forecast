#include <Rcpp.h>
#include <time.h>

using namespace Rcpp;

/* STEP 1: SLICE SAMPLER DEFINITION */

NumericVector slice_sampler(double (*logfn)(NumericVector, NumericVector),
                            NumericVector params,
                            NumericVector x0,
                            int steps = 10,
                            double w = 1,
                            double lower = -INFINITY,
                            double upper = INFINITY) {
  
  double u, r0, r1, logy, logz, logys;
  NumericVector x, xs, L, R;
  
  x = clone(x0);
  L = clone(x0);
  R = clone(x0);
  logy = logfn(x, params);
  
  for (int i = 0; i < steps; i++) {
    
    for (int j = 0; j < x0.size(); j++) {
      // draw uniformly from [0, y]
      logz = logy - rexp(1)[0];
      
      // expand search range
      u = runif(1)[0] * w;
      L[j] = x[j] - u;
      R[j] = x[j] + (w-u);
      while ( L[j] > lower && logfn(L, params) > logz )
        L[j] = L[j] - w;
      while ( R[j] < upper && logfn(R, params) > logz )
        R[j] = R[j] + w;
      
      // sample until draw is within valid range
      r0 = std::max(L[j], lower);
      r1 = std::min(R[j], upper);
      
      xs = clone(x);
      int cnt = 0;
      do {
        cnt++;
        xs[j] = runif(1, r0, r1)[0];
        logys = logfn(xs, params);
        if ( logys > logz )
          break;
        if ( xs[j] < x[j] )
          r0 = xs[j];
        else
          r1 = xs[j];
      } while (cnt<1e4);
      if (cnt==1e4) ::Rf_error("slice_sample_cpp loop did not finish");
      
      x = clone(xs);
      logy = logys;
    }
  }
  
  return x;
}

/* STEP 2: DEFINING THE POSTERIOR DISTRIBUTIONS:
I. PNBD-model.
1. posterior for gamma distributed heterogeneity parameters -> post_gamma_parameters,
2. posterior for lambda -> post_lambda
3. posterior for mu -> post_mu
*/

// I.1 PNBD: posterior for gamma distributed heterogeneity parameters -> post_gamma_parameters


double post_gamma_parameters(NumericVector log_data, NumericVector params) {
  double shape = exp(log_data[0]);
  double rate = exp(log_data[1]);
  double len_x = params[0];
  double sum_x = params[1];
  double sum_log_x = params[2];
  double hyper1 = params[3];
  double hyper2 = params[4];
  double hyper3 = params[5];
  double hyper4 = params[6];
  return len_x * (shape * log(rate) - lgamma(shape)) + (shape-1) * sum_log_x - rate * sum_x +
    (hyper1 - 1) * log(shape) - (shape * hyper2) +
    (hyper3 - 1) * log(rate) - (rate * hyper4);  //gemeinsame Posterior von r/alpha bzw. s/beta, Ma/Liu Seite (2f) Formel (5)/(6)
}


// [[Rcpp::export]]
NumericVector slice_sampler_gamma_parameters(
    NumericVector data,
    NumericVector init,
    NumericVector hyper,
    double steps = 20,
    double w = 1)
{
  NumericVector params = NumericVector::create(data.size(), sum(data), sum(log(data)), hyper[0], hyper[1], hyper[2], hyper[3]);
  return exp(slice_sampler(post_gamma_parameters, params, log(init), steps, w, -INFINITY, INFINITY));
}

// I.2, I.3 PNBD: posteriors for lambda and mu

double post_lambda_pnbd(NumericVector data, NumericVector params) {
  double lambda_ = data[0];
  double x       = params[0];
  double tx      = params[1];
  double Tcal    = params[2];
  //  double lambda = params[3];
  double mu      = params[4];
  double r       = params[5];
  double alpha   = params[6];
  //  double s       = params[7];
  //  double beta    = params[8];
  if ( log(mu+lambda_) - log(mu) < 1e-10 ) {
    return -INFINITY; // avoid numeric underflow
  } else {
    return (r-1+x) * log(lambda_) - (lambda_*alpha) - log(lambda_+mu) +
      log(mu*exp(-tx*(lambda_+mu))+lambda_*exp(-Tcal*(lambda_+mu)));  // Ma/Liu (2007) Seite 2 Formel (3)
  }
}

double post_mu_pnbd(NumericVector data, NumericVector params) {
  double mu_    = data[0];
  // double x      = params[0];
  double tx     = params[1];
  double Tcal   = params[2];
  double lambda = params[3];
  //  double mu     = params[4];
  //  double r      = params[5];
  //  double alpha  = params[6];
  double s      = params[7];
  double beta   = params[8];
  if ( log(lambda+mu_) - log(lambda) < 1e-10 ) {
    return -INFINITY; // avoid numeric underflow
  } else {
    return (s-1) * log(mu_) - (mu_*beta) - log(lambda+mu_) +
      log(mu_*exp(-tx*(lambda+mu_))+lambda*exp(-Tcal*(lambda+mu_)));   // Ma/Liu (2007) Seite 2 Formel (4)
  }
}


// [[Rcpp::export]]
NumericVector slice_sample_pnbd(String what,
                                NumericVector x, NumericVector tx, NumericVector Tcal,
                                NumericVector lambda, NumericVector mu,
                                double r, double alpha, double s, double beta) {
  int N = x.size();
  NumericVector out(N);
  for (int i=0; i<N; i++) {
    // Rcpp::Rcout << i << " x:" << x[i] << " tx:" << tx[i] << " Tcal:" << Tcal[i] << " lambda:" << lambda[i] << " mu:" << mu[i] << " r:" << r << " alpha:" << alpha << " s:" << s << " beta:" << beta << " - " << std::endl;
    NumericVector params = NumericVector::create(x[i], tx[i], Tcal[i], lambda[i], mu[i], r, alpha, s, beta);
    if (what == "lambda") {
      out[i] = slice_sampler(post_lambda_pnbd, params, NumericVector::create(lambda[i]), 3, 3 * sqrt(r) / alpha, 1e-5, 1e+5)[0];
    } else if (what == "mu") {
      out[i] = slice_sampler(post_mu_pnbd, params, NumericVector::create(mu[i]), 6, 3 * sqrt(s) / beta, 1e-5, 1e+5)[0];
    }
  }
  return out;
}

double simpson38(double (*fn)(double), double a, double b) {
  // http://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson.27s_3.2F8_rule_.28for_n_intervals.29
  double n = 12.0;
  double integral = (3.0/8.0) * ((b-a)/n) *
    (fn(a) +
    3 * fn(a+(1/n)*(b-a)) +
    3 * fn(a+(2/n)*(b-a)) +
    2 * fn(a+(3/n)*(b-a)) +
    3 * fn(a+(4/n)*(b-a)) +
    3 * fn(a+(5/n)*(b-a)) +
    2 * fn(a+(6/n)*(b-a)) +
    3 * fn(a+(7/n)*(b-a)) +
    3 * fn(a+(8/n)*(b-a)) +
    2 * fn(a+(9/n)*(b-a)) +
    3 * fn(a+(10/n)*(b-a)) +
    3 * fn(a+(11/n)*(b-a)) +
    fn(b));
  return(integral);
}



