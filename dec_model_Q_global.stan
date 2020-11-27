data {
  int<lower=0> N; // N is the number of lines in the inut table. More than one per site at the moment, maybe can be more efficient than this
  real<lower=0> max_sd;
  vector[N] time_gone;
  vector[N] mass_left;
  vector[N] u0;
  vector[N] tmax; //tmax is taken deterministically from the diameter at death, but an error is added later
  vector[N] scaled_d; // scaled wood density to correct tmax
}

parameters {
    //real<lower=0> tmax;
    real<lower=0> beta;
    real<lower=0> eta_11;
    real<lower=0> q0;
    real<lower=0> e0;
    real<lower=0> fc;
    real<lower=0> delay;
    real<lower=0> tmax_error;
    real<lower=0> u0_error;

  }

transformed parameters {
    real<lower=1, upper=1.6> zeta; // declare zeta
    zeta = (1-e0)/(beta*eta_11*e0);
 }


model {
    //real zeta; // declare zeta

    //zeta = (1-e0)/(beta*eta_11*e0);
    beta        ~ normal(7,7*0.1);
    eta_11      ~ normal(0.36,0.36*0.1);
    q0          ~ normal(1.101848,0.1275926);
    e0          ~ normal(0.30,0.25);
    fc          ~ normal(0.5,0.5*0.1);
    delay       ~ uniform(0,5);
    tmax_error  ~ normal(1,1*0.15);
    u0_error    ~ normal(1,1*0.5);

    for (i in 1:N){ // loop for each data point (material remamining and time coordinates)
      vector[N] alpha; // alpha chages every N
      vector[N] time_dec;

      alpha[i] = fc*beta*eta_11*u0[i]*u0_error*q0^beta; // alpha is recalculated every N because of varying u0 in different sites
      time_dec[i]=time_gone[i]-delay;
      if(time_dec[i]<=0)
        time_dec[i]=0;

      // decomposition model running
      if(time_dec[i]<(tmax[i]*scaled_d[i]*tmax_error)){
      mass_left[i] ~ normal(((2/(tmax[i]*scaled_d[i]*tmax_error))*(1/(alpha[i]*(1-zeta)))*((1+alpha[i]*time_dec[i])^(1-zeta)-
                           (1-(time_dec[i]/(tmax[i]*scaled_d[i]*tmax_error))))+
                           ((2/(tmax[i]*scaled_d[i]*tmax_error)^2)*(1/(alpha[i]^2*(1-zeta)*(2-zeta)))*(1-(1+alpha[i]*time_dec[i])^(2-zeta)))+
                           (1-(time_dec[i]/(tmax[i]*scaled_d[i]*tmax_error)))^2),
                           // after this the deviation of the mass left (error)
                           max_sd); //the sigma is the standard deviation in the data
      } else {
      mass_left[i]~ normal((2/(tmax[i]*scaled_d[i]*tmax_error))*(1/(alpha[i]*(1-zeta)))*(1+alpha[i]*time_dec[i])^(1-zeta)+
                           ((2/((tmax[i]*scaled_d[i]*tmax_error)^2))*(1/(alpha[i]^2*(1-zeta)*(2-zeta)))*(((1+alpha[i]*(time_dec[i]-(tmax[i]*scaled_d[i]*tmax_error)))^(2-zeta))-((1+alpha[i]*time_dec[i])^(2-zeta)))),
                           // after this the deviation of the mass left (error)
                           max_sd); //the sigma is the standard deviation in the data

      }
    }
}
