
pacman::p_load(tidyverse, 
               here,
               patchwork,
               ggplot2,
               pander,
               brms,
               metafor,
               # reading xls
               readxl,
               metafor)

### an example of a location-scale model
dat <- dat.bangertdrowns2004

# 
dat$ni100 <- dat$ni/100
dat$si <- sqrt(dat$vi)

vcv <- diag(dat$vi)

# id is ID for e and m 
rownames(vcv) <- dat$id
colnames(vcv) <- dat$id


#####################
# random effect model
#####################

# metafor 1
ma1 <- rma(yi, vi, data=dat)
summary(ma1)

# metafor 2
ma2 <- rma.mv(yi = yi, V = vcv, random = ~ 1 | id, data = dat)
summary(ma2)

# brms 1

form3 <- bf(yi | se(si, sigma=TRUE) ~ 1 +  (1|id) # this is e
)

prior3 <- default_prior(form3, 
                        data = dat, 
                        #data2 = list(vcv = vcv),
                        family = gaussian()
)

ma3 <- brm(form3, 
           data = dat, 
           chains = 2, 
           cores = 2, 
           iter = 3000, 
           warmup = 1000,
           prior = prior3
)

print(summary(ma3), digits = 4)

# brms 2

form4 <- bf(yi  ~ 1 +  (1|id)  + fcor(vcv)
)

prior4 <- default_prior(form4, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

ma4 <- brm(form4, 
           data = dat, 
           data2 = list(vcv = vcv),
           chains = 2, 
           cores = 2, 
           iter = 3000, 
           warmup = 2000,
           prior = prior4
)

# brms 3


form5 <- bf(yi  ~ 1 +   (1|gr(id, cov = vcv)) # this is m
            # residual will be es_ID SD
)

prior5 <- default_prior(form5, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior5$prior[3] = "constant(1)"
prior5 

ma5 <- brm(form5, 
           data = dat, 
           data2 = list(vcv = vcv),
           chains = 2, 
           cores = 2, 
           iter = 6000, 
           warmup = 3000,
           prior = prior5
)

print(summary(ma5), digits = 4)


# metafor

mr1 <- rma(yi, vi, mods = ~ ni100, scale = ~ ni100 , data=dat)
summary(mr1)

# brms

form2 <- bf(yi ~ 1  + ni100 + (1|gr(id, cov = vcv)), # this is m
            sigma ~ 1 + ni100
)


prior2 <- default_prior(form2, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)


# fixing the varaince to 1 (meta-analysis)
prior2$prior[5] = "constant(1)"
prior2 
# fit model

mr2 <- brm(form2, 
           data = dat, 
           data2 = list(vcv = vcv),
           chains = 2, 
           cores = 2, 
           iter = 20000, 
           warmup = 10000,
           #backend = "cmdstanr",
           prior = prior2,
           #threads = threading(9),
           control = list(adapt_delta = 0.99, max_treedepth = 15)
)

print(summary(mr2), digits = 4)
