# some tests

pacman::p_load(tidyverse, 
               here,
               patchwork,
               ggplot2,
               pander,
               brms,
               metafor,
               blsmeta,
               # reading xls
               readxl,
               metafor)

# TODO

# - phylogenetic models?


#############
# EXAMPLE 1
#############

# Pottier, P., Burke, S., Zhang, R. Y., Noble, D. W., Schwanz, L. E., Drobniak, S. M., \& Nakagawa, S. (2022). 
# Developmental plasticity in thermal tolerance: Ontogenetic variation, persistence, and future directions. 
# Ecology Letters, 25(10), 2245-2268.

# patrice's data set (categorical - biological)

dat <- read.csv(here("data", "thermal.csv"))

head(dat)

dat$si <- sqrt(dat$Var_dARR)

vcv <- vcalc(vi = Var_dARR, 
             cluster = population_ID, 
             obs = es_ID, rho = 0.5, 
             data = dat)
rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID


# use metafor

# random effect model
# note it models tau2 (simga_u^2) but brms models sigma_u 

mod <- rma(yi = dARR, 
           vi = Var_dARR, 
           test="t",
           data = dat)

summary(mod)


modb <- rma(yi = dARR, 
           vi = Var_dARR, 
           mods = ~ habitat,
           test="t",
           data = dat)

summary(modb)

modc <- rma(yi = dARR, 
            vi = Var_dARR, 
            mods = ~ habitat,
            scale = ~ habitat,
            test="t",
            data = dat)

summary(modc)

# brms

# meta-analysis 

# no correlation 

form0 <- bf(dARR
            ~ 1   +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + (1|study_ID)
)


prior0 <- default_prior(form0, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)


# fixing the varaince to 1 (meta-analysis)
prior0$prior[3] = "constant(1)"
prior0 
# fit model

fit0 <- brm(form0, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 6000, 
            warmup = 3000,
            #backend = "cmdstanr",
            prior = prior0,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

# load the rds

fit0 <- readRDS(here("Rdata", "fit0.rds"))

summary(fit0)

# save this as rds

saveRDS(fit0, here("Rdata", "fit0.rds"))

# with correlation

form0b <- bf(dARR
            ~ 1   +
              (1|p|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + (1|p|study_ID)
)



prior0b <- default_prior(form0b, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior0b$prior[5] = "constant(1)"
prior0b 
# fit model

fit0b <- brm(form0b, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 6000, 
            warmup = 3000,
            #backend = "cmdstanr",
            prior = prior0b,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit0b)

# save this as rds

saveRDS(fit0b, here("Rdata", "fit0b.rds"))



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

fit1 <- brm(form1, 
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

saveRDS(fit1, here("Rdata", "fit1.rds"))
# read in rds

fit1 <- readRDS(here("Rdata", "fit1.rds"))


summary(fit1)




# blsmeta

res0 <- blsmeta(yi = dARR,
               vi = Var_dARR,
               es_id = es_ID,
               study_id = study_ID,
               mods_scale2 = ~ 1,
               data = dat)



# response difference
# we could do a version without correlation

form1a <- bf(dARR
            ~ 1 + habitat + 
              (1|p|study_ID) + 
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + habitat + 
              (1|p|study_ID)
)


# create prior

prior1a <- default_prior(form1a, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior1a$prior[7] = "constant(1)"
prior1a


# fit model

fit1a <- brm(form1a, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 6000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior1a,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit1a)

# save this as rds

saveRDS(fit1a, here("Rdata", "fit1a.rds"))

#TODO 


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

fit1b <- brm(form1b, 
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

summary(fit1b)

# save this as rds

saveRDS(fit1b, here("Rdata", "fit1b.rds"))

# reading rsd

fit1b <- readRDS(here("Rdata", "fit1b.rds"))

summary(fit1b)

# blsmeta

fit_mr7 <- blsmeta(yi = dARR, 
                   vi = Var_dARR, 
                   es_id = es_ID,
                   mods = ~ habitat,
                   mods_scale2 = ~ habitat,
                   data = dat)
saveRDS(fit_mr7, here("Rdata", "fit_mr7.rds"))


# 
# #############################
# # genetic data set (categorical - method)
# ######################
# 
# dat <- read_excel(here("data", "genetic.xlsx"))
# 
# head(dat)
# 
# dat$study_ID <- as.factor(dat$Reference)
# dat$es_ID <- as.factor(1:nrow(dat))
# 
# 
# vcv <- vcalc(vi = Var_d, 
#              cluster = study_ID, 
#              obs = es_ID, rho = 0, 
#              data = dat)
# rownames(vcv) <- dat$es_ID
# colnames(vcv) <- dat$es_ID
# 
# # SMD = d 
# form2 <- bf(d
#             ~ 1  + method +
#               (1|study_ID) + # this is u
#               (1|gr(es_ID, cov = vcv)), # this is m
#             sigma ~ 1 + method
# )
# 
# 
# # create prior
# 
# prior2 <- default_prior(form2, 
#                         data = dat, 
#                         data2 = list(vcv = vcv),
#                         family = gaussian()
# )
# 
# # fixing the varaince to 1 (meta-analysis)
# prior2$prior[6] = "constant(1)"
# prior2 
# # fit model
# 
# fit2 <- brm(form2, 
#             data = dat, 
#             data2 = list(vcv = vcv),
#             chains = 2, 
#             cores = 2, 
#             iter = 3000, 
#             warmup = 2000,
#             #backend = "cmdstanr",
#             prior = prior2,
#             #threads = threading(9),
#             control = list(adapt_delta = 0.95, max_treedepth = 15)
# )
# 
# summary(fit2)
# 
# # save this as rds
# 
# saveRDS(fit2, here("Rdata", "fit2.rds"))


#############
# EXAMPLE 2
#############

# Midolo, G., De Frenne, P., Hölzel, N., \& Wellstein, C. (2019). 
# Global patterns of intraspecific leaf trait responses to elevation. 
# Global change biology, 25(7), 2485-2498.

#############################
# latitiude - (continious)
###########################

dat <- read.csv(here("data", "elevation.csv"))
dim(dat)
# filter data Narea
dat <- dat %>% filter(trait == "Narea")

dim(dat)

dat$study_ID <- as.factor(dat$Study_ID)
dat$es_ID <- as.factor(1:nrow(dat))

# using escalc use ROM

dat <- escalc(measure = "ROM", 
              m1i = treatment, 
              m2i = control, 
              sd1i = sd_treatment, 
              sd2i = sd_control, 
              n1i = n_treatment, 
              n2i = n_control, 
              data = dat)


vcv <- diag(dat$vi)
rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID



# random effect model
# note it models tau2 (simga_u^2) but brms models sigma_u 

mod2 <- rma(yi = yi, 
           vi = vi, 
           test="t",
           data = dat)

summary(mod2)


mod2b <- rma(yi = yi, 
            vi = vi, 
            mods = ~ elevation_log,
            test="t",
            data = dat)

summary(mod2b)

mod2c <- rma(yi = yi, 
            vi = vi, 
            mods = ~ elevation_log,
            scale = ~ elevation_log,
            test="t",
            data = dat)

summary(mod2c)

# brms

# meta-analysis 

# no correlation 

form2 <- bf(yi
            ~ 1   +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + (1|study_ID)
)



prior2 <- default_prior(form2, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)


# fixing the varaince to 1 (meta-analysis)
prior2$prior[3] = "constant(1)"
prior2 
# fit model

fit2 <- brm(form2, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 10000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior2,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2)

# save this as rds

saveRDS(fit2, here("Rdata", "fit2.rds"))

# with correlation

form2b <- bf(yi
             ~ 1   +
               (1|p|study_ID) + # this is u
               (1|gr(es_ID, cov = vcv)), # this is m
             sigma ~ 1 + (1|p|study_ID)
)



prior2b <- default_prior(form2b, 
                         data = dat, 
                         data2 = list(vcv = vcv),
                         family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior2b$prior[5] = "constant(1)"
prior2b 
# fit model

fit2b <- brm(form2b, 
             data = dat, 
             data2 = list(vcv = vcv),
             chains = 2, 
             cores = 2, 
             iter = 16000, 
             warmup = 10000,
             prior = prior2b,
             control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2b)

# save this as rds

saveRDS(fit2b, here("Rdata", "fit2b.rds"))




# SMD = d 
form3 <- bf(yi
            ~ 1  + elevation_log +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + elevation_log
)


# create prior


prior3 <- default_prior(form3, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior3$prior[6] = "constant(1)"
prior3 
# fit model

fit3 <- brm(form3, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 6000, 
            warmup = 3000,
            prior = prior3,
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit3)

# save this as rds

saveRDS(fit3, here("Rdata", "fit3.rds"))

# reading rds

fit3 <- readRDS(here("Rdata", "fit3.rds"))

summary(fit3)

# metafor 

mr3_mod <- rma.mv(yi = yi, 
                 V =  vcv,
                 mod = ~ elevation_log,
                 random = list(~1 | study_ID,
                               ~1 | es_ID),
                 #struct = "DIAG",
                 data = dat, 
                 test = "t",
                 sparse = TRUE,
                 control=list(optimizer="optim", optmethod="Nelder-Mead")
)

bubble_plot(mr3_mod, mod = "elevation_log", group = "study_ID", g = TRUE, 
            xlab = "ln(elevation difference) (moderator)",
            ylab = "SDM (effect size)",
            k.pos = "top.left",
            legend.pos = "bottom.left")

###

dat$study_ID <- as.numeric(dat$study_ID)
fit_mr33 <- blsmeta(yi = yi,
                    vi = vi,
                    es_id = es_ID,
                    study_id = study_ID,
                    mods = ~ 1 + elevation_log,
                    mods_scale2 = ~ 1 + elevation_log,
                    data = dat)
#save
saveRDS(fit_mr33, here("Rdata", "fit_mr33.rds"))

# ##################
# # grazing
# ###################
# 
# dat <- read.csv(here("data", "grazing1.csv"))
# 
# head(dat)
# 
# 
# dat$study_ID <- as.factor(dat$study)
# dat$es_ID <- as.factor(1:nrow(dat))
# dat$se <- sqrt(dat$var.eff.size)
# dat$lnRR <- dat$eff.size
# dat$cyear <- scale(dat$study.year, center = TRUE, scale = FALSE)
# 
# vcv <- diag(dat$var.eff.size)
# 
# rownames(vcv) <- dat$es_ID
# colnames(vcv) <- dat$es_ID
# 
# 
# # Zr 
# 
# form4 <- bf(lnRR
#             ~ 1  + cyear + se +
#               (1|study_ID) + # this is u
#               (1|gr(es_ID, cov = vcv)), # this is m
#             sigma ~ 1 +cyear + se
# )
# 
# # prior
# 
# prior4 <- default_prior(form4, 
#                         data = dat, 
#                         data2 = list(vcv = vcv),
#                         family = gaussian()
# )
# 
# prior4$prior[6] = "constant(1)"
# prior4 
# 
# # fit model
# 
# fit4 <- brm(form4, 
#             data = dat, 
#             data2 = list(vcv = vcv),
#             chains = 2, 
#             cores = 2, 
#             iter = 35000, 
#             warmup = 5000,
#             #backend = "cmdstanr",
#             prior = prior4,
#             #threads = threading(9),
#             control = list(adapt_delta = 0.99, max_treedepth = 20)
# )
# 
# summary(fit4)
# 
# 
# # probably not use???
# 
# ##################
# # plant data set
# ###################
# 
# dat <- read_csv(here("data", "plant_dat.csv"))
# 
# head(dat)
# 
# # filtering the one with data for r_survival
# 
# dat <- dat %>% filter(!is.na(r_survival))
# 
# # transform correlation to Zr
# 
# dat$Zr <- atanh(dat$r_survival)
# 
# dat$study_ID <- as.factor(dat$Reference)
# dat$es_ID <- as.factor(1:nrow(dat))
# 
# # leave the last 4 digit of Reference
# 
# dat$year <- as.factor(substr(dat$Reference,  nchar(dat$Reference) - 3, nchar(dat$Reference)))
# 
# # from this we need to filter out "tion" and 997a and 997b needs to be replaced by 1997
# 
# dat$year <- ifelse(dat$year == "tion", NA, dat$year)
# dat$year <- ifelse(dat$year == "997a", "1997", dat$year)
# dat$year <- ifelse(dat$year == "997b", "1997", dat$year)
# 
# # filter out the NA
# 
# dat <- dat %>% filter(!is.na(year))
# 
# dat$cyear <- scale(dat$year, center = TRUE, scale = FALSE)
# dat$sqrt.inv.n <- sqrt(1/dat$n_survival)
# dat$vi <- 1/(dat$n_survival - 3)
# 
# vcv <- diag(dat$vi)
# 
# rownames(vcv) <- dat$es_ID
# colnames(vcv) <- dat$es_ID
# 
# 
# # Zr 
# 
# form4 <- bf(Zr
#             ~ 1  + cyear + sqrt.inv.n +
#               (1|study_ID) + # this is u
#               (1|gr(es_ID, cov = vcv)), # this is m
#             sigma ~ 1 +cyear + sqrt.inv.n
# )
# 
# # prior
# 
# prior4 <- default_prior(form4, 
#                         data = dat, 
#                         data2 = list(vcv = vcv),
#                         family = gaussian()
# )
# 
# 
# # fit model
# 
# fit4 <- brm(form4, 
#             data = dat, 
#             data2 = list(vcv = vcv),
#             chains = 2, 
#             cores = 2, 
#             iter = 35000, 
#             warmup = 5000,
#             #backend = "cmdstanr",
#             prior = prior4,
#             #threads = threading(9),
#             control = list(adapt_delta = 0.99, max_treedepth = 20)
# )
# 
# summary(fit4)

#############
# EXAMPLE 3
#############

##################
# publication bias
###################

# Neuschulz, E. L., Mueller, T., Schleuning, M., \& Böhning-Gaese, K. (2016). 
# Pollination and seed dispersal are the most threatened processes of plant regeneration. 
# Scientific Reports, 6(1), 29839.

dat <- read.csv(here("data", "seed.csv"))

head(dat)

dat$study_ID <- as.factor(dat$study)
dat$es_ID <- as.factor(1:nrow(dat))
dat$se <- sqrt(dat$var.eff.size)
dat$smd <- dat$eff.size
dat$cyear <- scale(dat$study.year, center = TRUE, scale = FALSE)

vcv <- diag(dat$var.eff.size)

rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID


# no correlation 
# this runs but the one with correlation does not mix well

form4 <- bf(smd
            ~ 1   +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + (1|study_ID)
)



prior4 <- default_prior(form4, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)


# fixing the varaince to 1 (meta-analysis)
prior4$prior[3] = "constant(1)"
prior4 
# fit model

fit4 <- brm(form4, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 30000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior2,
            #threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit4)

# save this as rds

saveRDS(fit4, here("Rdata", "fit4.rds"))

# with correlation
# a good example of not working with the correlation paramter

form4b <- bf(smd
             ~ 1   +
               (1|p|study_ID) + # this is u
               (1|gr(es_ID, cov = vcv)), # this is m
             sigma ~ 1 + (1|p|study_ID)
)



prior4b <- default_prior(form4b, 
                         data = dat, 
                         data2 = list(vcv = vcv),
                         family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior4b$prior[5] = "constant(1)"
prior4b 
# fit model

fit4b <- brm(form4b, 
             data = dat, 
             data2 = list(vcv = vcv),
             chains = 2, 
             cores = 2, 
             iter =35000, 
             warmup = 5000,
             #backend = "cmdstanr",
             prior = prior4b,
             #threads = threading(9),
             control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit4b)

# save this as rds

saveRDS(fit4b, here("Rdata", "fit4b.rds"))




# meta-regression
# SMD

form5 <- bf(smd
            ~ 1  + cyear + se +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 +cyear + se
)

# prior

prior5 <- default_prior(form5, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

prior5$prior[6] = "constant(1)"
prior5

# fit model

fit5 <- brm(form5, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 35000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior5,
            #threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(fit5) 

# save fit5 as rds

saveRDS(fit5, here("Rdata", "fit5.rds"))

fit5 <- readRDS(here("Rdata", "fit5.rds"))

summary(fit5)

##

dat$study_ID <- as.numeric(dat$study_ID)
fit_mr44 <- blsmeta(yi = smd,
                    vi = var.eff.size,
                    es_id = es_ID,
                    study_id = study_ID,
                    mods = ~ 1 + se + cyear,
                    mods_scale2 = ~ 1 + se + cyear,
                    data = dat)
#save
saveRDS(fit_mr44, here("Rdata", "fit_mr44.rds"))

#

form_ma4b <- bf(smd
                ~ 1   +
                  (1|study_ID) + # this is u
                  (1|gr(es_ID, cov = vcv)), # this is m
                sigma ~ 1
)



prior_ma4b <- default_prior(form_ma4b, 
                            data = dat, 
                            data2 = list(vcv = vcv),
                            family = gaussian()
)


# fixing the varaince to 1 (meta-analysis)
prior_ma4b$prior[3] = "constant(1)"
prior_ma4b 
# fit model

fit_ma4b <- brm(form_ma4b, 
                data = dat, 
                data2 = list(vcv = vcv),
                chains = 2, 
                cores = 2, 
                iter = 30000, 
                warmup = 5000,
                #backend = "cmdstanr",
                prior = prior_ma4b,
                #threads = threading(9),
                control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit_ma4b)

# save this as rds
saveRDS(fit_ma4b, here("Rdata", "fit_ma4b.rds"))


# #-----------
# # case_441
# #-------------
# 
# dat <- read.csv(here("data", "case_441.csv"))
# 
# head(dat)
# 
# dat$study_ID <- as.factor(dat$study)
# dat$es_ID <- as.factor(1:nrow(dat))
# dat$se <- sqrt(dat$var.eff.size)
# dat$lnrr <- dat$eff.size
# 
# vcv <- diag(dat$var.eff.size)
# 
# rownames(vcv) <- dat$es_ID
# colnames(vcv) <- dat$es_ID
# 
# # lnRR
# 
# form6 <- bf(lnrr
#             ~ 1  + cyear + se +
#               (1|study_ID) + # this is u
#               (1|gr(es_ID, cov = vcv)), # this is m
#             sigma ~ 1 +cyear + se
# )
# 
# # prior
# 
# prior6 <- default_prior(form6, 
#                         data = dat, 
#                         data2 = list(vcv = vcv),
#                         family = gaussian()
# )
# 
# prior6$prior[6] = "constant(1)"
# prior6
# 
# # fit model
# 
# fit6 <- brm(form6, 
#             data = dat, 
#             data2 = list(vcv = vcv),
#             chains = 2, 
#             cores = 2, 
#             iter = 6000, 
#             warmup = 3000,
#             #backend = "cmdstanr",
#             prior = prior6,
#             #threads = threading(9),
#             control = list(adapt_delta = 0.99, max_treedepth = 20)
# )
# 
# summary(fit6)
# 
# # save fit6 as rds
# 
# saveRDS(fit6, here("Rdata", "fit6.rds"))
# 
# fit6 <- readRDS(here("Rdata", "fit6.rds"))

