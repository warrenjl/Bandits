case_study_fun<-function(n_sim,
                         final_y,
                         final_m){

###############################
#Seed
###############################
set.seed(24897)

y_sim_tot<-
y_sim_reg<-matrix(NA,
                  nrow = n_sim,
                  ncol = 10)
for(sim in 1:n_sim){
  
###########################################
#Simulation Information
###########################################
n_locs<-ncol(final_y)
n_days<-nrow(final_y)

ee_times<-c((n_locs*1),  #1 Visit Per Zone
            (n_locs*2))  #2 Visits Per Zone

#############################################################################################################################
#Simulating Study
#############################################################################################################################
y_se<-
y_eg<-
y_ucb<-
y_rs<-
y_ts<-
y_rts_p<-
y_rts_nb<-
y_c<-rep(NA,
         length = n_days)
y_ee<-list(0)
for(j in 1:length(ee_times)){
   y_ee[[j]]<-rep(NA,
                  times = n_days)
   }

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
  
counter<-1
active_set<-c(1:n_locs)
g_se<-rep(NA,
          times = n_days)
epsilon<-0.05
g_eg<-rep(NA,
          times = n_days)
g_ucb<-rep(NA,
           times = n_days)
g_ee<-rep(NA,
          times = length(ee_times))
   
for(i in 1:n_days){
     
   #######################
   #Complete Data on Day i
   #######################
   m<-final_m[i,]
   y<-final_y[i,]
   
   ####################################################################################
   #Clairvoyant
   ####################################################################################
   y_c[i]<-y[c(1:length(unique(data$Town)))[colSums(final_y) == max(colSums(final_y))]]
   
   ######################################################################################
   #Successive Elimination:  From "Introduction to Multi-Armed Bandits"
   ######################################################################################
   if(counter <= length(active_set)){
     
     g_se[i]<-active_set[counter]
     y_se[i]<-y[g_se[i]]
     counter<-counter +
              1
     
     }
   
   if(counter > length(active_set)){
     
     lcb<-
     ucb<-rep(NA,
              times = length(active_set))
     for(j in 1:length(active_set)){
       
       k_m<-sum(y_se[c(1:i)][g_se[1:i] == active_set[j]])/sum(g_se[1:i] == active_set[j])
       ucb[j]<-k_m +
               sqrt(2*log(i)/sum(g_se[1:i] == active_set[j]))
       lcb[j]<-k_m -
               sqrt(2*log(i)/sum(g_se[1:i] == active_set[j]))
       
       }
     
     delete<-rep(NA,
                 times = length(active_set))     
     for(j in 1:length(active_set)){
        delete[j]<-max(as.numeric(ucb[j] < lcb))
        }
     active_set<-active_set[delete == 0]
     counter<-1
     
     }
   
   ##################################################################################
   #epsilon Greedy:  From "Algorithms for the multi-armed bandit problem"
   ##################################################################################
   if(i <= n_locs){
     
     g_eg[i]<-i
     y_eg[i]<-y[g_eg[i]]
     
     }
   
   if(i > n_locs){
     
     dec_val<-rep(NA,
                  times = n_locs)
     for(j in 1:n_locs){
        dec_val[j]<-sum(y_eg[c(1:(i-1))][g_eg[1:(i-1)] == j])/sum(g_eg[1:(i-1)] == j)
        }
     
     choice<-rbinom(n = 1,
                    size = 1,
                    prob = epsilon)
     if(choice == 0){
       g_eg[i]<-c(1:n_locs)[dec_val == max(dec_val)][1]
       }
     
     if(choice == 1){
       g_eg[i]<-sample(c(1:n_locs),
                       size = 1)
       }
     
     y_eg[i]<-y[g_eg[i]]
     
     }
   
   ##############################################################################
   #UCB1:  From "Algorithms for the multi-armed bandit problem"
   ##############################################################################
   if(i <= n_locs){
     
     g_ucb[i]<-i
     y_ucb[i]<-y[g_ucb[i]]
     
     }
   
   if(i > n_locs){
     
     dec_val<-rep(NA,
                  times = n_locs)
     for(j in 1:n_locs){
       
        k_m<-sum(y_ucb[c(1:(i-1))][g_ucb[1:(i-1)] == j])/sum(g_ucb[1:(i-1)] == j)
        dec_val[j]<-k_m +
                    sqrt(2*log(i-1)/sum(g_ucb[1:(i-1)] == j))
       
        }
     
     g_ucb[i]<-c(1:n_locs)[dec_val == max(dec_val)][1]
     y_ucb[i]<-y[g_ucb[i]]
     
     }
   
   #########################
   #Random Sampling
   #########################
   g_rs<-sample(c(1:n_locs),
                size = 1)
   y_rs[i]<-y[g_rs]
      
   ##############################################################
   #Explore then Exploit
   ##############################################################
   for(j in 1:length(ee_times)){
         
      if(i <= ee_times[j]){
        
        g_ee[j]<-i%%n_locs
        g_ee[j][g_ee[j] == 0]<-n_locs
        y_ee[[j]][i]<-y[g_ee[j]]
           
        }
        
      if(i == ee_times[j]){
             
        ee_tot<-rep(NA,
                    times = n_locs)
        for(k in 1:n_locs){
           ee_tot[k]<-sum(y_ee[[j]][seq(k, ee_times[j], n_locs)])
           }
        g_ee[j]<-c(1:n_locs)[ee_tot == max(ee_tot)][1]
          
        }
      
      if(i > ee_times[j]){
        y_ee[[j]][i]<-y[g_ee[j]]
        }
        
      }
      
   ########################################
   #Thompson Sampling
   ########################################
   draw_p<-rbeta(n = n_locs,
                 shape1 = a_beta_ts,
                 shape2 = b_beta_ts)
   g_ts<-c(1:n_locs)[draw_p == max(draw_p)]
   y_ts[i]<-y[g_ts]
      
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
   y_rts_p[i]<-y[g_rts_p]
      
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
   for(j in 1:n_locs){
      draw_r[j]<-sample(c(1:s),
                        prob = w_cat_rts_nb[j,],
                        size = 1)
      }
   draw_pi<-rbeta(n = n_locs,
                  shape1 = c_beta_rts_nb,
                  shape2 = (d_beta_rts_nb + e_beta_rts_nb*draw_r))
   draw_p<-rbeta(n = n_locs,
                 shape1 = a_beta_rts_nb,
                 shape2 = b_beta_rts_nb)
   g_rts_nb<-c(1:n_locs)[(draw_pi*draw_r*draw_p/(1.00 - draw_pi)) == max((draw_pi*draw_r*draw_p/(1.00 - draw_pi)))]
   y_rts_nb[i]<-y[g_rts_nb]
      
   log_w<-rep(NA,
              times = s)
   for(j in 1:s){
      log_w[j]<-log(w_cat_rts_nb[g_rts_nb, j]) +
                log(choose((m[g_rts_nb] + j - 1), m[g_rts_nb])) + 
                lbeta((m[g_rts_nb] + c_beta_rts_nb[g_rts_nb]), (d_beta_rts_nb[g_rts_nb] + j*(e_beta_rts_nb[g_rts_nb] + 1))) +
                -lbeta(c_beta_rts_nb[g_rts_nb], (d_beta_rts_nb[g_rts_nb] + j*e_beta_rts_nb[g_rts_nb]))
      }
   for(j in 1:s){
      w_cat_rts_nb[g_rts_nb, j]<-1.00/sum(exp(log_w - log_w[j]))
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
      
   } 

y_sim_tot[sim, 1]<-sum(y_rts_p)
y_sim_tot[sim, 2]<-sum(y_rts_nb)
y_sim_tot[sim, 3]<-sum(y_ee[[1]])
y_sim_tot[sim, 4]<-sum(y_ee[[2]])
y_sim_tot[sim, 5]<-sum(y_rs)
y_sim_tot[sim, 6]<-sum(y_ts)
y_sim_tot[sim, 7]<-sum(y_ucb)
y_sim_tot[sim, 8]<-sum(y_eg)
y_sim_tot[sim, 9]<-sum(y_se)
y_sim_tot[sim, 10]<-sum(y_c)

y_sim_reg[sim, 1]<-sum(y_c - y_rts_p)
y_sim_reg[sim, 2]<-sum(y_c - y_rts_nb)
y_sim_reg[sim, 3]<-sum(y_c - y_ee[[1]])
y_sim_reg[sim, 4]<-sum(y_c - y_ee[[2]])
y_sim_reg[sim, 5]<-sum(y_c - y_rs)
y_sim_reg[sim, 6]<-sum(y_c - y_ts)
y_sim_reg[sim, 7]<-sum(y_c - y_ucb)
y_sim_reg[sim, 8]<-sum(y_c - y_eg)
y_sim_reg[sim, 9]<-sum(y_c - y_se)
y_sim_reg[sim, 10]<-sum(y_c - y_c)

print(sim/n_sim)

}

colMeans(y_sim_tot)
colMeans(y_sim_reg)

return(list(y_sim_tot,
            y_sim_reg))

}
