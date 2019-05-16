# age-structure SEIR model with apc-based force of infection
# Apr 12, 2019

#loading packages####
packages = c("tidyverse", "pomp2", "reshape2", "plyr", "ggthemes",
             "grid", "gridExtra")

sapply(packages, library, character.only = TRUE)
rm(packages)

# read.csv("sample_covar.csv") -> sample_covar
read.csv("sample_data.csv") -> sample_data
matrix(scan('foi_matrix.txt'),ncol=5,byrow=TRUE) -> foi_mat
matrix(scan('pops.txt'),ncol=5,byrow=TRUE) -> pops

Csnippet("
  //double *X = &X1;
  double *Clocal = &C1;
  int i, j;
  int n = agegroups;
  int time;
  //int m = duration;
  double foi;
  double *FOIlocal = get_userdata_double(\"foi_mat\");
  double *Slocal = &X1;
  double *Elocal = Slocal + n;
  double *Ilocal = Elocal + n;
  double *Rlocal = Ilocal + n;
  // age K, year J
  #define FOI(J,K) FOIlocal[((J)-1)*n+(K)]
  #define S(K) Slocal[(K)]
  #define E(K) Elocal[(K)]
  #define I(K) Ilocal[(K)]
  #define R(K) Rlocal[(K)]
  #define C(K) Clocal[(K)]
  // transitions between compartments
  double dSE[n], dEI[n], dIE[n], dIR[n], dRE[n];
  // natural deaths
  double dSD[n], dED[n], dID[n], dRD[n];
  //transitions between age groups/recruiments and loss
  double dSS[n], dEE[n], dII[n], dRR[n];
  double *L = &lambda1;
  //Rprintf(\"%lg %d %lg %lg %lg %lg %lg\\n\",t, time, L[0], L[1], L[2], L[3], L[4]);

  // epidemic process
  for (i=0; i<n; i++) {
    // transitions between compartments
    time = (int) (floor(t/12) + 1);     // FIX ME
    foi = FOI(time,i);
    //Rprintf(\"%lg %d %lg\\n\",t, time, foi);
    dSE[i] = rbinom(S(i), foi*dt);
    dEI[i] = rbinom(E(i), alpha*dt);
    dIE[i] = rbinom(I(i), phi*dt);
    dIR[i] = rbinom(I(i), gamma*dt);
    dRE[i] = rbinom(R(i), L[i]*dt);
    // natural deaths
    dSD[i] = rbinom(S(i), mu*dt);
    dED[i] = rbinom(E(i), mu*dt);
    dID[i] = rbinom(I(i), mu*dt);
    dRD[i] = rbinom(R(i), mu*dt);
    // recruiments and loss
    dSS[0] = rbinom(N0, p*dt) - rbinom(S(0), p*dt);
    dEE[0] = -rbinom(E(0), p*dt);
    dII[0] = -rbinom(I(0), p*dt);
    dRR[0] = -rbinom(R(0), p*dt);
    if (i>0) {
      dSS[i] = rbinom(S(i-1), p*dt) - rbinom(S(i), p*dt);
      dEE[i] = rbinom(E(i-1), p*dt) - rbinom(E(i), p*dt);
      dII[i] = rbinom(I(i-1), p*dt) - rbinom(I(i), p*dt);
      dRR[i] = rbinom(R(i-1), p*dt) - rbinom(R(i), p*dt);
      //Rprintf(\"%lg %d %lg %lg %lg %lg\\n\",t, time, S(i-1), S(i), R(i-1), R(i));
    }
    S(i) += dSS[i] - dSE[i] - dSD[i];
    E(i) += dEE[i] + dSE[i] + dIE[i] + dRE[i] - dEI[i] - dED[i];
    I(i) += dII[i] + dEI[i] - dIE[i] - dIR[i] - dID[i];
    R(i) += dRR[i] + dIR[i] - dRE[i] - dRD[i];
    C(i) += dIR[i];    // this is C
    //Rprintf(\"%lg %d %lg %lg %lg %lg %lg\\n\",t, time, S(i), E(i), I(i), R(i), C(i));
  }
") -> rproc

# pops: 539936.6, 647457.1, 526054.5, 302406.8, 102332.7

Csnippet("
  double *X = &X1;
  const double *P = get_userdata_double(\"pops\");
  int n = agegroups;
  const double *ratio = &S1_0;   // 20, S, E, I, R

  for (int i=0; i<n; i++) {                   
    X[i] = nearbyint(ratio[i*n+1]*P[i]);                      //S(i)
    X[n+i] = nearbyint(ratio[i*n+2]*P[i]);                    //E(i)
    X[2*n+i] = nearbyint(ratio[i*n+3]*P[i]);                  //I(i)
    X[3*n+i] = nearbyint(ratio[i*n+4]*P[i]);                  //R(i)
  }
") -> rinit

agegroups <- 5
# stateX <- c(sprintf("X%d", 1:(4*agegroups)))

Csnippet("
  double *case = &C1;
  double *d = &data1;
  int n = agegroups;
  lik = 1;
  for (i=0; i<n; i++) {
    lik = lik * dnbinom_mu(d[i],theta,rho*case[i],give_log);
  }
  lik = (give_log) ? log(lik) : lik;
") -> dmeas

Csnippet("
  //double *case = &C1;
  //double *d = &data1;
  //int n = agegroups;
  //for (i=0; i<n; i++) {
  //  d[i] = rnbinom_mu(theta,rho*case[i]);
  //}
  data1 = rnbinom_mu(theta,rho*C1);
  data2 = rnbinom_mu(theta,rho*C2);
  data3 = rnbinom_mu(theta,rho*C3);
  data4 = rnbinom_mu(theta,rho*C4);
  data5 = rnbinom_mu(theta,rho*C5);
") -> rmeas


params <- c(0.0068/12,1/3,0.001/12, 1/2, 
            rep(0.001,5),
            # 0.001,
            0.95, 15, 1/15/12, 5, 453677,
            rep(c(0.6,0.3,0.003,1-0.9-0.003), times=5)
            )
# rw <- rep(0.02,length(params))
paramnames <- c("mu", "alpha", "phi", "gamma",
                # "lambda",
                sprintf("lambda%d", 1:agegroups),
                "rho", "theta", "p", "agegroups", "N0",
                paste0(rep(c("S", "E", "I", "R"),5), rep(1:5, each=4),"_0")
                )
names(params) <- paramnames


sample_data %>%
  # filter(time < 5) %>%
  pfilter(times="time", t0=0, Np=1000,
    # monthly params
    params=params,
    rinit=rinit,
    rprocess=euler(rproc, delta.t=1/10),
    dmeasure=dmeas,
    rmeasure=rmeas,
    accumvars=c(sprintf("C%d", 1:agegroups)),
    statenames=c(sprintf("X%d", 1:(4*agegroups)),sprintf("C%d", 1:agegroups)),
    paramnames=paramnames,
    foi_mat = foi_mat,
    pops =  c(611955, 695555, 546655, 261520, 88495) # pops in 2006
) -> pf

# pf %>% as.data.frame() %>%
#   gather(var,val,-time) %>%  
#   filter(!is.na(val)) %>%
#   ggplot(aes(x=time,y=val))+
#   geom_line()+
#   facet_wrap(~var,scales="free_y")

ests <- c("alpha", "phi", "gamma", 
          # "lambda",
          sprintf("lambda%d", 1:agegroups),
          paste0(rep(c("S", "E", "I", "R"),5), rep(1:5, each=4),"_0"),
          "rho", "theta"
          )

bake(file="mf_test.rds",seed=186593237,
     {
       pf %>%
         mif2(Nmif=50,rw.sd=rw.sd(
           alpha=0.02, phi=0.02, gamma=0.02, 
           # lambda=0.02,
           lambda1=0.02, lambda2=0.02, lambda3=0.02, lambda4=0.02, lambda5=0.02,
           S1_0=ivp(0.02), E1_0=ivp(0.02), I1_0=ivp(0.02), R1_0=ivp(0.02), 
           S2_0=ivp(0.02), E2_0=ivp(0.02), I2_0=ivp(0.02), R2_0=ivp(0.02), 
           S3_0=ivp(0.02), E3_0=ivp(0.02), I3_0=ivp(0.02), R3_0=ivp(0.02), 
           S4_0=ivp(0.02), E4_0=ivp(0.02), I4_0=ivp(0.02), R4_0=ivp(0.02), 
           S5_0=ivp(0.02), E5_0=ivp(0.02), I5_0=ivp(0.02), R5_0=ivp(0.02), 
           rho=0.02, theta=2
           ),
           partrans=parameter_trans(
             logit=c("phi",
                     sprintf("lambda%d", 1:agegroups),
                     "rho"),
             log=c("theta", "gamma", "alpha"),
             toEst=Csnippet("to_log_barycentric(&T_S1_0,&S1_0,4);
                             to_log_barycentric(&T_S2_0,&S2_0,4);
                             to_log_barycentric(&T_S3_0,&S3_0,4);
                             to_log_barycentric(&T_S4_0,&S4_0,4);
                             to_log_barycentric(&T_S5_0,&S5_0,4);"),
             fromEst=Csnippet("from_log_barycentric(&S1_0,&T_S1_0,4);
                               from_log_barycentric(&S2_0,&T_S2_0,4);
                               from_log_barycentric(&S3_0,&T_S3_0,4);
                               from_log_barycentric(&S4_0,&T_S4_0,4);
                               from_log_barycentric(&S5_0,&T_S5_0,4);")
           ),
           paramnames=ests,
           cooling.type="geometric",cooling.fraction.50=0.9) -> mf1
     }) -> mf1

mf1 %>% 
  traces() %>% 
  # as.data.frame() %>%
  # arrange(-loglik) %>% head()
  melt() %>%
  ggplot(aes(x=iteration,y=value))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")

