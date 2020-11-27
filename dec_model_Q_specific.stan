data {
  int<lower=0> N; // N is the number of lines in the inut table. More than one per site at the moment, maybe can be more efficient than this
  real<lower=0> max_sd; 
  vector[N] time_gone;
  vector[N] mass_left;
  vector[N] u0;
  vector[N] tmax; //tmax is taken deterministically from the diameter at death, but an error is added later
  vector[N] scaled_d; // scaled wood density to correct tmax
  vector[N] decay_tree_class; // type of the decaying material
  vector[N] tree_class; // wood type (for different q0)
  vector[N] delay_class; // delay class (snag or log have different delay)
}

parameters {

    real<lower=5.5, upper=8.5> beta;
    real<lower=0.36-0.1, upper=0.36+0.1> eta_11;
    real<lower=0.5, upper=1.5> q0_b;
    real<lower=0.5, upper=1.5> q0_s;
    real<lower=0.5, upper=1.5> q0_p;
    real<lower=0.5, upper=1.5> q0_t;
    real<lower=0.1, upper=0.6> e0;
    real<lower=0.4, upper=0.6> fc;
    real<lower=0, upper=5> delayS; //delay for snags
    real<lower=0, upper=5> delayL; // delay for logs
    

    // starting to declare the local parameters here
    real<lower=0.5, upper=1.5> tmax_error1;
    real<lower=0.5, upper=1.5> tmax_error2;
    real<lower=0.5, upper=1.5> tmax_error3;
    real<lower=0.5, upper=1.5> tmax_error4;
    real<lower=0.5, upper=1.5> tmax_error5;
    real<lower=0.5, upper=1.5> tmax_error6;
    real<lower=0.5, upper=1.5> tmax_error7;
    real<lower=0.5, upper=1.5> tmax_error8;
    real<lower=0.5, upper=1.5> tmax_error9;
    real<lower=0.5, upper=1.5> tmax_error10;

    real<lower=0.5, upper=1.5> u0_error1;
    real<lower=0.5, upper=1.5> u0_error2;
    real<lower=0.5, upper=1.5> u0_error3;
    real<lower=0.5, upper=1.5> u0_error4;
    real<lower=0.5, upper=1.5> u0_error5;
    real<lower=0.5, upper=1.5> u0_error6;
    real<lower=0.5, upper=1.5> u0_error7;
    real<lower=0.5, upper=1.5> u0_error8;
    real<lower=0.5, upper=1.5> u0_error9;
    real<lower=0.5, upper=1.5> u0_error10;
    


  }
  
transformed parameters { 
    real<lower=1, upper=1.6> zeta; // declare zeta
    zeta = (1-e0)/(beta*eta_11*e0);
 } 


model {
    real tmax_error; // declare tmax error to be overwritten
    real q0; // declare q0 error to be overwritten
    
    beta        ~ normal(7,7*0.1);
    eta_11      ~ normal(0.36,0.36*0.1);
    q0_b        ~ normal(1.101848,0.1275926);
    q0_s        ~ normal(1.101848,0.1275926);
    q0_p        ~ normal(1.101848,0.1275926);
    q0_t        ~ normal(1.101848,0.1275926);
    e0          ~ normal(0.55,0.25);
    fc          ~ normal(0.5,0.5*0.1);
    delayS      ~ uniform(0,5);
    delayL      ~ uniform(0,5);

    // starting to declare the local parameters here
    u0_error1   ~ normal(1,1*0.5);
    u0_error2   ~ normal(1,1*0.5);
    u0_error3   ~ normal(1,1*0.5);
    u0_error4   ~ normal(1,1*0.5);
    u0_error5   ~ normal(1,1*0.5);
    u0_error6   ~ normal(1,1*0.5);
    u0_error7   ~ normal(1,1*0.5);
    u0_error8   ~ normal(1,1*0.5);
    u0_error9   ~ normal(1,1*0.5);
    u0_error10  ~ normal(1,1*0.5);
    
    tmax_error1   ~ normal(1,1*0.5);
    tmax_error2   ~ normal(1,1*0.5);
    tmax_error3   ~ normal(1,1*0.5);
    tmax_error4   ~ normal(1,1*0.5);
    tmax_error5   ~ normal(1,1*0.5);
    tmax_error6   ~ normal(1,1*0.5);
    tmax_error7   ~ normal(1,1*0.5);
    tmax_error8   ~ normal(1,1*0.5);
    tmax_error9   ~ normal(1,1*0.5);
    tmax_error10  ~ normal(1,1*0.5);

    
    for (i in 1:N){ // loop for each data point (material remamining and time coordinates)
      vector[N] alpha; // alpha chages every N
      vector[N] time_dec;
      vector[N] SS;

  // if statements for q0
            if(tree_class[i]==3){ 
              q0  = q0_b;} // if wood is birch
              else if(tree_class[i]==2){
              q0  = q0_s;} // if wood is spruce
              else if(tree_class[i]==1){
              q0  = q0_p;} // if wood is pine
              else if(tree_class[i]==4){
              q0  = q0_t;} // if the data is from Tarasov et al.
            

  // if statements for u0 error
            if(decay_tree_class[i]==3){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error1*q0^beta; // alpha is recalculated every N because of varying u0 in different sites
              tmax_error  = tmax_error1;}
            else if(decay_tree_class[i]==2){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error2*q0^beta;
              tmax_error  = tmax_error2;}
            else if(decay_tree_class[i]==1){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error3*q0^beta;
              tmax_error  = tmax_error3;}
            else if(decay_tree_class[i]==5){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error4*q0^beta;
              tmax_error  = tmax_error4;}
            else if(decay_tree_class[i]==7){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error5*q0^beta;
              tmax_error  = tmax_error5;}
            else if(decay_tree_class[i]==6){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error6*q0^beta;
              tmax_error  = tmax_error6;}
            else if(decay_tree_class[i]==9){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error7*q0^beta;
              tmax_error  = tmax_error7;}
            else if(decay_tree_class[i]==10){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error8*q0^beta;
              tmax_error  = tmax_error8;}
            else if(decay_tree_class[i]==11){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error9*q0^beta;
              tmax_error  = tmax_error9;}
            else if(decay_tree_class[i]==16){
              alpha[i] = fc*beta*eta_11*u0[i]*u0_error10*q0^beta;
              tmax_error  = tmax_error10;}
          
          
          //adding a new condition for the maximum steady state
           SS[i]=((1/alpha[i])*(1/(zeta-1)))+((tmax[i]*scaled_d[i]*tmax_error)/3);

   // if statements for delay
            if(delay_class[i]==1){
              time_dec[i]=time_gone[i]-delayS;} // when the wood is a snag of any kind
            else if(delay_class[i]==0){
              time_dec[i]=time_gone[i]-delayL;} // when the wood is a log of any kind
  
        // resetting the time to zero in case delay makes it negative
        if(time_dec[i]<=0)
        time_dec[i]=0;




      // decomposition model running
      if(time_dec[i]<(tmax[i]*scaled_d[i]*tmax_error)){
      mass_left[i] ~ normal(((2/(tmax[i]*scaled_d[i]*tmax_error))*(1/(alpha[i]*(1-zeta)))*((1+alpha[i]*time_dec[i])^(1-zeta)-
                          (1-(time_dec[i]/(tmax[i]*scaled_d[i]*tmax_error))))+
                          ((2/(tmax[i]*scaled_d[i]*tmax_error)^2)*(1/(alpha[i]^2*(1-zeta)*(2-zeta)))*(1+alpha[i]*time_dec[i])^(2-zeta))+
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
