# Setup BRMS for Erick

rm(list = ls())

pacman::p_load(devtools, 
               tidyverse, 
               metafor, 
               patchwork, 
               R.rsp, 
               emmeans,
               orchaRd,
               metafor,
               dplyr,
               broom,
               tidyverse,
               here,
               broom,
               ggplot2,
               data.table,
               Matrix, # added
               matrixcalc,# added
               ape, # added
               multcomp,
               rotl,#added
               brms
) 

dat <- read.csv(here("data", "full_dataset_brms.csv")) 

# getting tree
tree_pred <- read.tree(here("Tree", "tree_brms1.tre"))
tree_prey <- read.tree(here("Tree", "tree_brms2.tre"))

# getting branch length and correlation matrix
tree1b <- compute.brlen(tree_pred)
cor1 <-  vcv(tree1b, corr=T) #predator - corresponding col is Predator_1 for RE

tree2b <- compute.brlen(tree_prey)
cor2 <-  vcv(tree2b, corr=T) #prey - corresponding col is Prey_1 for RE

#Non-phylo RE can be Predator2_1 and Prey2_1

#Your three VCV's
#Predator-preyoverlap
#vcv1 <- vcalc(vi_diff, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
             # data = dat)

#Predator single species
vcv1 <- vcalc(Predator_vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
              data = dat)

#Prey single species
vcv2 <- vcalc(Prey_vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
              data = dat)

# Two models below - none of them are perfect but both shows strong correlation between predator and prey

# NOW THIS MODEL IS GOOD

# simple correaliton model assumoing corealiton at residual level

form_prey <- bf(Prey_shift_rev_logit ~ 
                    1 + 
                  (1 |p| Study_ID) + 
                 # (1 | Obs_ID) + # this is included in brms
                  (1 | gr(Predator_1, cov = cor1)) + # phylo effect of Predactor
                  (1 | Predator2_1) + # non-phylo effect
                  (1 | Predator2_1) + # non-phylo effect
                  (1 | gr(Prey_1, cov = cor2)) + # phylo effect of Prey
                  (1 | Prey2_1) + # non-phylo effect of Prey
                  fcor(vcv2)  # a way to get VCV into brms
)


form_pred <- bf(Predator_shift_rev_logit ~ 
                    1 + 
                  (1 |p| Study_ID) + 
                  #(1 | Obs_ID) + # correlation between prey + pred
                  (1 | gr(Predator_1, cov = cor1)) + # phylo effect of Predactor
                  (1 | Predator2_1) + # non-phylo effect
                  (1 | gr(Prey_1, cov = cor2)) + # phylo effect of Prey
                  (1 | Prey2_1) + # non-phylo effect of Prey
                  fcor(vcv1)# a way to get VCV into brms
)                

from <- form_prey + form_pred + set_rescor(FALSE)

# 

prior0 <- get_prior(from,
                    data = dat,
                    data2 = list(vcv1 = vcv1, vcv2 = vcv2, cor1 = cor1, cor2 = cor2),
                    family = gaussian())

#set fcor to have variance of 1
prior0$prior[16] <- "constant(1)"
prior0$prior[29] <- "constant(1)"

#prior0

# # Define priors
# prior <- c(
#   prior(constant(1), class = "sigma")
# ) #set fcor to have variance of 1


out <- brm(from, 
           family = gaussian(),
           data = dat, 
           data2 = list(vcv1 = vcv1, vcv2 = vcv2, cor1 = cor1, cor2 = cor2),
           chains = 2, 
           cores = 2, 
           iter = 4000, 
           warmup = 1000,
           prior = prior0
           #set_rescor100(rescor = FLASE)
)
out


# saving the result as an RDS file

saveRDS(out, here("R", "out.rds"))

# reading RDS

out <- readRDS(here("R", "out.rds"))

out

# DO NOT RUN below


# I like this one better
# alternative model

form_prey <- bf(Prey_shift_rev_logit| se(sqrt(Prey_vi)) ~ 
                   1 + 
                  (1 |p| Study_ID) + # this is where correlation is expected
                  (1 | Obs_ID) + # this is within study level (no correlation expected)
                  (1 | gr(Predator_1, cov = cor1)) + # phylo effect of Predactor
                  (1 | Predator2_1) + # non-phylo effect
                  (1 | Predator2_1) + # non-phylo effect
                  (1 | gr(Prey_1, cov = cor2)) + # phylo effect of Prey
                  (1 | Prey2_1) #+ # non-phylo effect of Prey
                  #fcor(vcv2)  # a way to get VCV into brms
)


form_pred <- bf(Predator_shift_rev_logit | se(sqrt(Predator_vi)) ~ 
                   1 + 
                  (1 |p| Study_ID) + 
                  (1 | Obs_ID) + # correlation between prey + pred
                  (1 | gr(Predator_1, cov = cor1)) + # phylo effect of Predactor
                  (1 | Predator2_1) + # non-phylo effect
                  (1 | gr(Prey_1, cov = cor2)) + # phylo effect of Prey
                  (1 | Prey2_1) #+ # non-phylo effect of Prey
                  #fcor(vcv1)# a way to get VCV into brms
)                

from <- form_prey + form_pred+ set_rescor(TRUE)

out2 <- brm(from, 
           family = gaussian(),
           data = dat, 
           data2 = list(vcv1 = vcv1, vcv2 = vcv2, cor1 = cor1, cor2 = cor2),
           chains = 2, 
           cores = 2, 
           iter = 4000, 
           warmup = 1000
           #set_rescor100(rescor = FLASE)
)
out2

# saving the result as an RDS file

saveRDS(out2, here("R", "out2.rds"))

# reading RDS

out2 <- readRDS(here("R", "out2.rds"))

out2 

##THinking about plotting this regression 

library(tidybayes)

library("functional")
library("scales")
inverse_logit_brks_trans <-
  trans_new("inverse logit",
            transform = plogis,
            inverse = qlogis,
            breaks = Compose(plogis, extended_breaks(), qlogis),
            format = Compose(plogis, format_format()))



bayes <- (model_fit <- dat %>%
    add_predicted_draws(out) %>%  # adding the posterior distribution
    ggplot(aes(x = Prey_shift_rev_logit, y = Predator_shift_rev_logit), transfm ="invlogit") +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                     alpha = 0.5, colour = "black") + 
      # geom_smooth(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
      #             #                 alpha = 0.5, colour = "black") ) + #
    geom_point(data = dat, colour = "darkseagreen4") +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("Predator temporal shift in response to disturbance") +  
    xlab("Prey temporal shift in response to disturbance") +
    theme_bw() +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

conditional_effects(out, effects = "Prey_shift_rev_logit:Predator_shift_rev_logit",
                    conditions = conditions)

bayes <- bayes + scale_y_continuous(trans = inverse_logit_brks_trans)

bayes <- bayes + scale_x_continuous(trans = inverse_logit_brks_trans)

bayes

ggsave("correlation_fig.pdf", width = 4, height = 4)




