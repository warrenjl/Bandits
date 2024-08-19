##############
#Seed
##############
set.seed(6822)

##################################
#Simulation Settings:
#p_setting: "same_p", "varying_p"
#r_setting:  "same_r", "varying_r"
##################################
p_setting<-"varying_p"
r_setting<-"varying_r"

pi_true<-10.00/11.00

###########################################
#Simulation Information
###########################################
n_locs<-10  
n_days<-10000
n_sims<-1000

#p: Republic of the Congo
a<-52.00*0.20
b<-(1000.00 - 52.00)*0.20

ee_times<-c((n_locs*1),  #1 Visit Per Zone
            (n_locs*5))  #5 Visits Per Zone

################################################################################################################################
#Simulating Study
################################################################################################################################
y_se<-
y_eg<-
y_ucb<-
y_rs<-
y_ts<-
y_rts_p<-
y_rts_nb<-
y_c<-
r_se<-
r_eg<-
r_ucb<-
r_rs<-
r_ts<-
r_rts_p<-
r_rts_nb<-
r_c<-matrix(NA,
            nrow = n_sims,
            ncol = n_days)
y_ee<-
r_ee<-list(0)
for(j in 1:length(ee_times)){
   y_ee[[j]]<-
   r_ee[[j]]<-matrix(NA,
                     nrow = n_sims,
                     ncol = n_days)
   }

for(i in 1:n_sims){
  
   ###################################
   #Prior Information
   ###################################
   a_beta_ts<-rep(1.00,
                  times = n_locs)
   b_beta_ts<-rep(1.00,
                  times = n_locs)
   
   a_beta_rts_p<-rep(1.00,
                     times = n_locs)
   b_beta_rts_p<-rep(1.00,
                     times = n_locs)
   c_gamma_rts_p<-rep(1.00,
                      times = n_locs)
   d_gamma_rts_p<-rep(0.0001,
                      times = n_locs)
  
   a_beta_rts_nb<-rep(1.00,
                      times = n_locs)
   b_beta_rts_nb<-rep(1.00,
                      times = n_locs)
   c_beta_rts_nb<-rep(1.00,
                      times = n_locs)
   d_beta_rts_nb<-rep(1.00,
                      times = n_locs)
   e_beta_rts_nb<-rep(0.00,
                      times = n_locs)
  
   s<-100
   w_cat_rts_nb<-matrix((1.00/s),
                        nrow = n_locs,
                        ncol = s)
  
   ############################################################################################################
   #True Parameters
   ############################################################################################################
   if(p_setting == "varying_p"){
     p_true<-rbeta(n = n_locs,
                   shape1 = a,
                   shape2 = b)
     }

   if(p_setting == "same_p"){
     p_true<-rep(0.052,
                 times = n_locs)
     }

   if(r_setting == "varying_r"){
     r_true<-exp(rnorm(n = n_locs,
                       mean = 0.00,
                       sd = 0.75))
     }

   if(r_setting == "same_r"){
     r_true<-rep(1,
                 times = n_locs)
     }
   
   counter<-1
   active_set<-c(1:n_locs)
   g_se<-rep(NA,
             times = n_days)
   epsilon<-0.05
   g_eg<-rep(NA,
              times = n_days)
   g_ucb<-rep(NA,
              times = n_days)
   g_c<-c(1:n_locs)[(pi_true*r_true*p_true/(1.00 - pi_true)) == max(pi_true*r_true*p_true/(1.00 - pi_true))][1]
   g_ee<-rep(NA,
             times = length(ee_times))
   
   for(j in 1:n_days){
     
      ###################################
      #Complete Data on Day j
      ###################################
      m<-rnbinom(n = n_locs,
                 size = r_true,
                 prob = (1.00 - pi_true))
      y<-rbinom(n = n_locs,
                size = m,
                prob = p_true)
      
      ##########################################################################################
      #Successive Elimination:  From "Introduction to Multi-Armed Bandits"
      ##########################################################################################
      if(counter <= length(active_set)){
        
        g_se[j]<-active_set[counter]
        y_se[i,j]<-y[g_se[j]]
        r_se[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                   (pi_true*r_true*p_true/(1.00 - pi_true))[g_se[j]]
        counter<-counter +
                 1
        
        }
      
      if(counter > length(active_set)){
        
        lcb<-
        ucb<-rep(NA,
                 times = length(active_set))
        for(k in 1:length(active_set)){
          
           k_m<-sum(y_se[i, c(1:j)][g_se[1:j] == active_set[k]])/sum(g_se[1:j] == active_set[k])
           ucb[k]<-k_m +
                   sqrt(2*log(j)/sum(g_se[1:j] == active_set[k]))
           lcb[k]<-k_m -
                   sqrt(2*log(j)/sum(g_se[1:j] == active_set[k]))
          
           }
   
        delete<-rep(NA,
                    times = length(active_set))     
        for(k in 1:length(active_set)){
           delete[k]<-max(as.numeric(ucb[k] < lcb))
           }
        active_set<-active_set[delete == 0]
        counter<-1
        
        }
      
      #####################################################################################
      #epsilon Greedy:  From "Algorithms for the multi-armed bandit problem"
      #####################################################################################
      if(j <= n_locs){
        
        g_eg[j]<-j
        y_eg[i,j]<-y[g_eg[j]]
        r_eg[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                   (pi_true*r_true*p_true/(1.00 - pi_true))[g_eg[j]]
        
        }
      
      if(j > n_locs){
        
        dec_val<-rep(NA,
                     times = n_locs)
        for(k in 1:n_locs){
           dec_val[k]<-sum(y_eg[i, c(1:(j-1))][g_eg[1:(j-1)] == k])/sum(g_eg[1:(j-1)] == k)
           }
        
        choice<-rbinom(n = 1,
                       size = 1,
                       prob = epsilon)
        if(choice == 0){
          g_eg[j]<-c(1:n_locs)[dec_val == max(dec_val)][1]
          }
        
        if(choice == 1){
          g_eg[j]<-sample(c(1:n_locs),
                          size = 1)
          }
        
        y_eg[i,j]<-y[g_eg[j]]
        r_eg[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                   (pi_true*r_true*p_true/(1.00 - pi_true))[g_eg[j]]
        
        }
      
      #################################################################################
      #UCB1:  From "Algorithms for the multi-armed bandit problem"
      #################################################################################
      if(j <= n_locs){
        
        g_ucb[j]<-j
        y_ucb[i,j]<-y[g_ucb[j]]
        r_ucb[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                    (pi_true*r_true*p_true/(1.00 - pi_true))[g_ucb[j]]
        
        }
      
      if(j > n_locs){
        
        dec_val<-rep(NA,
                     times = n_locs)
        for(k in 1:n_locs){
          
           k_m<-sum(y_ucb[i, c(1:(j-1))][g_ucb[1:(j-1)] == k])/sum(g_ucb[1:(j-1)] == k)
           dec_val[k]<-k_m +
                       sqrt(2*log(j-1)/sum(g_ucb[1:(j-1)] == k))
                       
           }
        
        g_ucb[j]<-c(1:n_locs)[dec_val == max(dec_val)][1]
        y_ucb[i,j]<-y[g_ucb[j]]
        r_ucb[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                    (pi_true*r_true*p_true/(1.00 - pi_true))[g_ucb[j]]
        
        }
      
      #########################################################
      #Random Sampling
      #########################################################
      g_rs<-sample(c(1:n_locs),
                   size = 1)
      y_rs[i,j]<-y[g_rs]
      r_rs[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                 (pi_true*r_true*p_true/(1.00 - pi_true))[g_rs]
      
      ######################################################################
      #Explore then Exploit
      ######################################################################
      for(k in 1:length(ee_times)){
         
         if(j <= ee_times[k]){
        
           g_ee[k]<-j%%n_locs
           g_ee[k][g_ee[k] == 0]<-n_locs
           y_ee[[k]][i,j]<-y[g_ee[k]]
           r_ee[[k]][i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                           (pi_true*r_true*p_true/(1.00 - pi_true))[g_ee[k]]
           
           }
        
         if(j == ee_times[k]){
             
           ee_tot<-rep(NA,
                       times = n_locs)
           for(l in 1:n_locs){
              ee_tot[l]<-sum(y_ee[[k]][i, seq(l, ee_times[k], n_locs)])
              }
           g_ee[k]<-c(1:n_locs)[ee_tot == max(ee_tot)][1]
          
           }
      
         if(j > ee_times[k]){
           
           y_ee[[k]][i,j]<-y[g_ee[k]]
           r_ee[[k]][i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                           (pi_true*r_true*p_true/(1.00 - pi_true))[g_ee[k]]
           
           }
        
         }
      
      #########################################################
      #Thompson Sampling
      #########################################################
      draw_p<-rbeta(n = n_locs,
                    shape1 = a_beta_ts,
                    shape2 = b_beta_ts)
      g_ts<-c(1:n_locs)[draw_p == max(draw_p)]
      y_ts[i,j]<-y[g_ts]
      r_ts[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                 (pi_true*r_true*p_true/(1.00 - pi_true))[g_ts]
      
      a_beta_ts[g_ts]<-a_beta_ts[g_ts] +
                       y[g_ts]
      b_beta_ts[g_ts]<-b_beta_ts[g_ts] +
                       (m[g_ts] - y[g_ts])
      
      #####################################################################
      #Random Thompson Sampling: Poisson
      #####################################################################
      draw_p<-rbeta(n = n_locs,
                    shape1 = a_beta_rts_p,
                    shape2 = b_beta_rts_p)
      draw_lambda<-rgamma(n = n_locs,
                          shape = c_gamma_rts_p,
                          rate = d_gamma_rts_p)
      g_rts_p<-c(1:n_locs)[(draw_lambda*draw_p) == max(draw_lambda*draw_p)]
      y_rts_p[i,j]<-y[g_rts_p]
      r_rts_p[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                    (pi_true*r_true*p_true/(1.00 - pi_true))[g_rts_p]
      
      a_beta_rts_p[g_rts_p]<-a_beta_rts_p[g_rts_p] +
                             y[g_rts_p]
      b_beta_rts_p[g_rts_p]<-b_beta_rts_p[g_rts_p] +
                             (m[g_rts_p] - y[g_rts_p])
      
      c_gamma_rts_p[g_rts_p]<-c_gamma_rts_p[g_rts_p] +
                              m[g_rts_p]
      d_gamma_rts_p[g_rts_p]<-d_gamma_rts_p[g_rts_p] +
                              1
      
      ##########################################################################################################################
      #Random Thompson Sampling: Negative Binomial
      ##########################################################################################################################
      draw_r<-rep(NA,
                  times = n_locs)
      for(k in 1:n_locs){
         draw_r[k]<-sample(c(1:s),
                           prob = w_cat_rts_nb[k,],
                           size = 1)
         }
      draw_pi<-rbeta(n = n_locs,
                     shape1 = c_beta_rts_nb,
                     shape2 = (d_beta_rts_nb + e_beta_rts_nb*draw_r))
      draw_p<-rbeta(n = n_locs,
                    shape1 = a_beta_rts_nb,
                    shape2 = b_beta_rts_nb)
      g_rts_nb<-c(1:n_locs)[(draw_pi*draw_r*draw_p/(1.00 - draw_pi)) == max((draw_pi*draw_r*draw_p/(1.00 - draw_pi)))]
      y_rts_nb[i,j]<-y[g_rts_nb]
      r_rts_nb[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                     (pi_true*r_true*p_true/(1.00 - pi_true))[g_rts_nb]
      
      log_w<-rep(NA,
                 times = s)
      for(k in 1:s){
         log_w[k]<-log(w_cat_rts_nb[g_rts_nb, k]) +
                   log(choose((m[g_rts_nb] + k - 1), m[g_rts_nb])) + 
                   lbeta((m[g_rts_nb] + c_beta_rts_nb[g_rts_nb]), (d_beta_rts_nb[g_rts_nb] + k*(e_beta_rts_nb[g_rts_nb] + 1))) +
                   -lbeta(c_beta_rts_nb[g_rts_nb], (d_beta_rts_nb[g_rts_nb] + k*e_beta_rts_nb[g_rts_nb]))
         }
      for(k in 1:s){
        w_cat_rts_nb[g_rts_nb, k]<-1.00/sum(exp(log_w - log_w[k]))
        }
      w_cat_rts_nb[g_rts_nb,][is.na(w_cat_rts_nb[g_rts_nb,]) == 1]<-0.00
      
      c_beta_rts_nb[g_rts_nb]<-c_beta_rts_nb[g_rts_nb] +
                               m[g_rts_nb]
      e_beta_rts_nb[g_rts_nb]<-e_beta_rts_nb[g_rts_nb] +
                               1
      
      a_beta_rts_nb[g_rts_nb]<-a_beta_rts_nb[g_rts_nb] +
                               y[g_rts_nb]
      b_beta_rts_nb[g_rts_nb]<-b_beta_rts_nb[g_rts_nb] +
                               (m[g_rts_nb] - y[g_rts_nb])
      
      #######################################################
      #Clairvoyant
      #######################################################
      y_c[i,j]<-y[g_c]
      r_c[i,j]<-max(pi_true*r_true*p_true/(1.00 - pi_true)) -
                (pi_true*r_true*p_true/(1.00 - pi_true))[g_c]
      
      } 
  
   print(round(100*(i/n_sims), 2))
   
   }

#############################################################################
#Summarizing Output
#############################################################################
library(matrixStats)

ci_factor<-qt(0.975,
              df = (n_sims - 1))

se_est<-rowMeans(apply(X = r_se, MARGIN = 1, FUN = cumsum))
se_se<-rowSds(apply(X = r_se, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
se_ci_lower<-se_est - 
             ci_factor*se_se
se_ci_upper<-se_est + 
             ci_factor*se_se

eg_est<-rowMeans(apply(X = r_eg, MARGIN = 1, FUN = cumsum))
eg_se<-rowSds(apply(X = r_eg, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
eg_ci_lower<-eg_est - 
             ci_factor*eg_se
eg_ci_upper<-eg_est + 
             ci_factor*eg_se

ucb_est<-rowMeans(apply(X = r_ucb, MARGIN = 1, FUN = cumsum))
ucb_se<-rowSds(apply(X = r_ucb, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
ucb_ci_lower<-ucb_est - 
              ci_factor*ucb_se
ucb_ci_upper<-ucb_est + 
              ci_factor*ucb_se

rs_est<-rowMeans(apply(X = r_rs, MARGIN = 1, FUN = cumsum))
rs_se<-rowSds(apply(X = r_rs, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
rs_ci_lower<-rs_est - 
             ci_factor*rs_se
rs_ci_upper<-rs_est + 
             ci_factor*rs_se

ts_est<-rowMeans(apply(X = r_ts, MARGIN = 1, FUN = cumsum))
ts_se<-rowSds(apply(X = r_ts, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
ts_ci_lower<-ts_est - 
             ci_factor*ts_se
ts_ci_upper<-ts_est + 
             ci_factor*ts_se

rts_p_est<-rowMeans(apply(X = r_rts_p, MARGIN = 1, FUN = cumsum))
rts_p_se<-rowSds(apply(X = r_rts_p, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
rts_p_ci_lower<-rts_p_est - 
                ci_factor*rts_p_se
rts_p_ci_upper<-rts_p_est + 
                ci_factor*rts_p_se

rts_nb_est<-rowMeans(apply(X = r_rts_nb, MARGIN = 1, FUN = cumsum))
rts_nb_se<-rowSds(apply(X = r_rts_nb, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
rts_nb_ci_lower<-rts_nb_est - 
                 ci_factor*rts_nb_se
rts_nb_ci_upper<-rts_nb_est + 
                 ci_factor*rts_nb_se

c_est<-rowMeans(apply(X = r_c, MARGIN = 1, FUN = cumsum))
c_se<-rowSds(apply(X = r_c, MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
c_ci_lower<-c_est - 
            ci_factor*c_se
c_ci_upper<-c_est + 
            ci_factor*c_se

ee_est<-
ee_ci_lower<-
ee_ci_upper<-list(0)
for(j in 1:length(ee_times)){
   
   ee_est[[j]]<-rowMeans(apply(X = r_ee[[j]], MARGIN = 1, FUN = cumsum))
   ee_se<-rowSds(apply(X = r_ee[[j]], MARGIN = 1, FUN = cumsum))/sqrt(n_sims)
   ee_ci_lower[[j]]<-ee_est[[j]] - 
                     ci_factor*ee_se
   ee_ci_upper[[j]]<-ee_est[[j]] + 
                     ci_factor*ee_se
  
   }




