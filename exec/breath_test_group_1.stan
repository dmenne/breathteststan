// Group fit of Breath Test Curves to Exponential beta with partial 
// hierarchical model

data{
  int<lower=0> n; // Number of data
  int<lower=0> n_pat; // Number of patients
  int<lower=0> n_group; // Number groups
  real<lower=0> dose;
  real<lower=2> student_t_df; // When < 10, use student, else normal
  int<lower=0> pat_i[n];
  int<lower=0> group_i[n];
  vector<lower=0>[n] minute;
  vector<lower=-30>[n] pdr;
}

parameters{
  vector[n_pat] m_pat_raw;
  vector[n_group] m_group;
  real<lower=0> sigma_m_pat;
  real<lower=0> mu_m;
  
  vector[n_pat] k_pat_raw;
  vector[n_group] k_group;
  real<lower=0> sigma_k_pat;
  real<lower=0> mu_k;
  
  vector[n_pat] beta_pat_raw;
  vector[n_group] beta_group;
  real<lower=0> sigma_beta_pat;
  real<lower=1> mu_beta;
  
  real <lower=0> sigma;
}

transformed parameters {
  vector[n] pdr1;
  vector[n_pat] m_pat;
  vector[n_pat] k_pat;
  vector[n_pat] beta_pat;
  real<lower = 0> mn; // Impose constraints
  real<lower = 0>  kn; 
  real<lower = 1>  betan;
  m_pat = sigma_m_pat * m_pat_raw;
  k_pat = sigma_k_pat * k_pat_raw;
  beta_pat = sigma_beta_pat * beta_pat_raw;
  for (i in 1:n){
    int pat;
    int group;
    real exp_ktn;
    pat = pat_i[i];
    group = group_i[i];
    mn  = mu_m + m_pat[pat] + m_group[group];
    kn = mu_k + k_pat[pat] + k_group[group];
    betan = mu_beta + beta_pat[pat] + beta_group[group]; 
    exp_ktn = exp(-kn* minute[i]);
    pdr1[i] = dose*mn*kn*betan*exp_ktn * pow(1 - exp_ktn,(betan -1));
  }
}

model {
  m_pat_raw ~ normal(0, 1);
  sigma_m_pat ~ normal(0, 10);
  m_group ~ normal(0, 10);
  mu_m ~ normal(40,20);
  
  k_pat_raw ~ normal(0, 1);
  sigma_k_pat ~ normal(0, 0.01);
  k_group ~ normal(0, 0.005);
  mu_k ~ normal(0.007, 0.002); // Approximately the default 1/0.65 h
  
  beta_pat_raw ~ normal(0, 1);
  sigma_beta_pat ~ normal(0, 1);
  beta_group ~ normal(0, 0.5);
  mu_beta ~ normal(2, 0.4);
  
  sigma ~ cauchy(0,5);
  if (student_t_df < 10) {
    pdr ~ student_t(student_t_df, pdr1, sigma); 
  } else {
    pdr ~ normal(pdr1, sigma);
  }
}

