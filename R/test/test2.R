# comparing metafor and brms

devtools::install_github("donaldRwilliams/blsmeta")

# some tests

pacman::p_load(tidyverse, 
               here,
               patchwork,
               ggplot2,
               pander,
               brms,
               metafor,
               # reading xls
               readxl,
               blsmeta,
               metafor)

# patrice's data set (categorical - biological)

dat <- read.csv(here("data", "thermal.csv"))

head(dat)

# sampling error
dat$si <- sqrt(dat$Var_dARR)


# we should use this for real one 
# vcv <- vcalc(vi = Var_dARR, 
#              cluster = population_ID, 
#              obs = es_ID, rho = 0.5, 
#              data = dat)

vcv <- diag(dat$Var_dARR)

# es_ID is ID for e and m 
# while study ID is ID for study
rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID


# use metafor

# random effect model
# note it models tau2 (simga_u^2) but brms models sigma_u 


# meta-analysis

#####################
# random effect model
#####################

# metafor 1
#------------
ma1 <- rma(yi = dARR, 
           vi = Var_dARR, 
           test="t",
           data = dat)

summary(ma1)


# metafor 2
#------------
ma2 <- rma.mv(yi = dARR, 
              V = vcv,
              random = ~ 1 | es_ID,
              test="t",
              data = dat)

summary(ma2)


# brms 1
#------------

form3 <- bf(dARR | se(si) ~ 1 +  (1|es_ID) # this is e
)

prior3 <- default_prior(form3, 
                        data = dat, 
                        #data2 = list(vcv = vcv),
                        family = gaussian()
)

ma3 <- brm(form3, 
            data = dat, 
            #data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 35000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior3,
            #threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 20)
)


# save as rds
saveRDS(ma3, here("Rdata", "ma3.rds"))

# load rds
ma3 <- readRDS(here("Rdata", "ma3.rds"))


print(summary(ma3), digits = 4)

# brms 2
#------------
form4 <- bf(dARR  ~ 1 +  (1|es_ID) # this is e
            + fcor(vcv)
)

prior4 <- default_prior(form4, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)


# we need to fix the variance to 1 = fcor(vcv)
prior4$prior[5] = "constant(1)"
prior4 

ma4 <- brm(form4, 
           data = dat, 
           data2 = list(vcv = vcv),
           chains = 2, 
           cores = 2, 
           iter = 15000, 
           warmup = 5000,
           #backend = "cmdstanr",
           prior = prior4,
           #threads = threading(9),
           #control = list(adapt_delta = 0.95, max_treedepth = 15)
)

# save as rds
saveRDS(ma4, here("Rdata", "ma4.rds"))

# load rds
ma4 <- readRDS(here("Rdata", "ma4.rds"))

print(summary(ma4), digit = 4)


# brm 3 (using Shinichi's trick)
#------------
form5 <- bf(dARR  ~ 1 +   (1|gr(es_ID, cov = vcv)), # this is m
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
           iter = 35000, 
           warmup = 5000,
           prior = prior5
)

print(summary(ma5), digits = 4)

# save as rds
saveRDS(ma4, here("Rdata", "ma5.rds"))

# above shows my trick works - let's move to l-s meta-regressoin

# metafor
#------------
mr1 <- rma(yi = dARR, 
            vi = Var_dARR, 
            mods = ~ habitat,
            scale = ~ habitat,
            test="t",
            data = dat)

summary(mr1)


# brms
#------------
form2 <- bf(dARR
            ~ 1   + habitat + 
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + habitat
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
            iter = 35000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior2,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

# save as rds

saveRDS(mr2, here("Rdata", "mr2.rds"))

# load the rds

mr2 <- readRDS(here("Rdata", "mr2.rds"))

summary(mr2)

# blmata


res0 <- blsmeta(yi = dARR,
               vi = Var_dARR,
               es_id = es_ID,
               study_id = study_ID,
               mods_scale2 = ~ 1,
               data = dat)


###############
# potential mistake - the model which does not run

# brms

form1 <- bf(dARR | se(si) ~ 1 + habitat, # this is e
            sigma ~ 1 + habitat
)

prior1 <- default_prior(form1, 
                        data = dat, 
                        family = gaussian()
)

ma1 <- brm(form1, 
           data = dat, 
           chains = 2, 
           cores = 2, 
           iter = 5000, 
           warmup = 2000,
           prior = prior1,
           control = list(adapt_delta = 0.99, max_treedepth = 20)
)

# this is wrong too - as it mulitple sigma^2 to sampling variance 
#(in meta-analysis, sampling varaince is assumed to be known)

form2 <- bf(dARR | se(si, sigma = TRUE) ~ 1 + habitat, # this is e
            sigma ~ 1 + habitat
)

prior2 <- default_prior(form1, 
                        data = dat, 
                        family = gaussian()
)

# this will run but this is a different model
ma2 <- brm(form2, 
           data = dat, 
           chains = 2, 
           cores = 2, 
           iter = 5000, 
           warmup = 2000,
           prior = prior1,
           control = list(adapt_delta = 0.99, max_treedepth = 20)
)

#############
## 3 levels
################

# biological - meta-regression
form1 <- bf(dARR
            ~ 1  + habitat +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + habitat
)


# create prior

prior1 <- default_prior(form1, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior1$prior[5] = "constant(1)"
prior1 
# fit model

res1 <- brm(form1, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior1,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

# save this as rds

saveRDS(res1, here("Rdata", "res1.rds"))

# read in rds
res1_brms <- readRDS(here("Rdata", "res1.rds"))


summary(res1_brms)

# blmeta

res1a <- blsmeta(yi = dARR,
               vi = Var_dARR,
               es_id = es_ID,
               study_id = study_ID,
               mods = ~ 1 + habitat,
               mods_scale2 = ~ 1 + habitat,
               data = dat)
#save

saveRDS(res1a, here("Rdata", "res1a.rds"))

# read in rds

res1a <- readRDS(here("Rdata", "res1a.rds"))

print(res1a)


# try method

dat$method <- ifelse(dat$exp_design == c("A", "B", "C"), "initial", "persisitent")



form1b <- bf(dARR
             ~ 1  + method +
               (1|study_ID) + # this is u
               (1|gr(es_ID, cov = vcv)), # this is m
             sigma ~ 1 + method
)


# create prior

prior1b <- default_prior(form1b, 
                         data = dat, 
                         data2 = list(vcv = vcv),
                         family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior1b$prior[5] = "constant(1)"
prior1b 
# fit model

res2 <- brm(form1b, 
             data = dat, 
             data2 = list(vcv = vcv),
             chains = 2, 
             cores = 2, 
             iter = 3000, 
             warmup = 2000,
             #backend = "cmdstanr",
             prior = prior1b,
             #threads = threading(9),
             control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(res2)

# save this as rds
saveRDS(res2, here("Rdata", "res2.rds"))
# reading rsd

res2 <- readRDS(here("Rdata", "res2.rds"))

summary(res2)

# blsmeta


res2a <- blsmeta(yi = dARR,
                 vi = Var_dARR,
                 es_id = es_ID,
                 study_id = study_ID,
                 mods = ~ 1 + method,
                 mods_scale2 = ~ 1 + method,
                 data = dat)
#save

#save

saveRDS(res2a, here("Rdata", "res2a.rds"))

# read in rds

res2a <- readRDS(here("Rdata", "res2a.rds"))

print(res2a)

#######

### Using `blsmeta`

```{r}
#| eval: false

fit_ma10 <- blsmeta(yi = dARR,
                    vi = Var_dARR,
                    es_id = es_ID,
                    study_id = study_ID,
                    mods_scale2 = ~ 1,
                    data = dat)

saveRDS(fit_ma10, here("Rdata", "fit_ma10.rds"))
```


```{r}

saveRDS(res1a, here("Rdata", "res1a.rds"))

# read in rds

res1a <- readRDS(here("Rdata", "res1a.rds"))

print(res1a)

```
