# some teests 

pacman::p_load(tidyverse, 
               here,
               patchwork,
               ggplot2,
               pander,
               brms,
               metafor,
               # reading xls
               readxl)

############
# load data
############
# test
# 
# data <- data.frame(
#   yi = c(0.2, 0.5, -0.3, 0.1),  # Effect sizes
#   sei = c(0.1, 0.2, 0.15, 0.1)  # Standard errors
# )
# 
# # Add a study column
# data$study <- factor(1:nrow(data))
# 
# # Add a study column
# library(brms)
# 
# # Fit the location-scale model
# meta_model <- brm(
#   bf(yi | se(sei) ~ 1 + (1 | study), sigma ~ 1 + (1 | study)),  # Location-scale model
#   data = data,
#   family = gaussian(),
#   prior = c(
#     prior(normal(0, 1), class = Intercept),
#     prior(cauchy(0, 1), class = sd),
#     prior(cauchy(0, 1), class = sigma)
#   ),
#   iter = 4000, warmup = 1000, chains = 4, cores = 4
# )


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

# length
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
prior1$prior[6] = "constant(1)"
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

summary(fit1)

# save this as rds

saveRDS(fit1, here("Rdata", "fit1.rds"))

# response difference


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
prior1$prior[8] = "constant(1)"
prior1 


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


#############################
# genetic data set (categorical - method)
######################

dat <- read_excel(here("data", "genetic.xlsx"))

head(dat)

dat$study_ID <- as.factor(dat$Reference)
dat$es_ID <- as.factor(1:nrow(dat))


vcv <- vcalc(vi = Var_d, 
             cluster = study_ID, 
             obs = es_ID, rho = 0, 
             data = dat)
rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID

# SMD = d 
form2 <- bf(d
            ~ 1  + method +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 + method
)


# create prior

prior2 <- default_prior(form2, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

# fixing the varaince to 1 (meta-analysis)
prior2$prior[6] = "constant(1)"
prior2 
# fit model

fit2 <- brm(form2, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2)

# save this as rds

saveRDS(fit2, here("Rdata", "fit2.rds"))


#############################
# latitiude - (continious)
###########################

dat <- read.csv(here("data", "elevation.csv"))

head(dat)

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
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior3,
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit3)

# save this as rds

saveRDS(fit3, here("Rdata", "fit3.rds"))


##################
# grazing
###################

dat <- read.csv(here("data", "grazing1.csv"))

head(dat)


dat$study_ID <- as.factor(dat$study)
dat$es_ID <- as.factor(1:nrow(dat))
dat$se <- sqrt(dat$var.eff.size)
dat$lnRR <- dat$eff.size
dat$cyear <- scale(dat$study.year, center = TRUE, scale = FALSE)

vcv <- diag(dat$var.eff.size)

rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID


# Zr 

form4 <- bf(lnRR
            ~ 1  + cyear + se +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 +cyear + se
)

# prior

prior4 <- default_prior(form4, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

prior4$prior[6] = "constant(1)"
prior4 

# fit model

fit4 <- brm(form4, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 35000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior4,
            #threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(fit4)


# probably not use???

##################
# plant data set
###################

dat <- read_csv(here("data", "plant_dat.csv"))

head(dat)

# filtering the one with data for r_survival

dat <- dat %>% filter(!is.na(r_survival))

# transform correlation to Zr

dat$Zr <- atanh(dat$r_survival)

dat$study_ID <- as.factor(dat$Reference)
dat$es_ID <- as.factor(1:nrow(dat))

# leave the last 4 digit of Reference

dat$year <- as.factor(substr(dat$Reference,  nchar(dat$Reference) - 3, nchar(dat$Reference)))

# from this we need to filter out "tion" and 997a and 997b needs to be replaced by 1997

dat$year <- ifelse(dat$year == "tion", NA, dat$year)
dat$year <- ifelse(dat$year == "997a", "1997", dat$year)
dat$year <- ifelse(dat$year == "997b", "1997", dat$year)

# filter out the NA

dat <- dat %>% filter(!is.na(year))

dat$cyear <- scale(dat$year, center = TRUE, scale = FALSE)
dat$sqrt.inv.n <- sqrt(1/dat$n_survival)
dat$vi <- 1/(dat$n_survival - 3)

vcv <- diag(dat$vi)

rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID


# Zr 

form4 <- bf(Zr
            ~ 1  + cyear + sqrt.inv.n +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 +cyear + sqrt.inv.n
)

# prior

prior4 <- default_prior(form4, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)


# fit model

fit4 <- brm(form4, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 35000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior4,
            #threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(fit4)

##################
# publication bias
###################


dat <- read.csv(here("data", "case_385.csv"))

head(dat)

dat$study_ID <- as.factor(dat$study)
dat$es_ID <- as.factor(1:nrow(dat))
dat$se <- sqrt(dat$var.eff.size)
dat$smd <- dat$eff.size
dat$cyear <- scale(dat$study.year, center = TRUE, scale = FALSE)

vcv <- diag(dat$var.eff.size)

rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID

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


# case_441
#-------------

dat <- read.csv(here("data", "case_441.csv"))

head(dat)

dat$study_ID <- as.factor(dat$study)
dat$es_ID <- as.factor(1:nrow(dat))
dat$se <- sqrt(dat$var.eff.size)
dat$lnrr <- dat$eff.size

vcv <- diag(dat$var.eff.size)

rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID

# lnRR

form6 <- bf(lnrr
            ~ 1  + cyear + se +
              (1|study_ID) + # this is u
              (1|gr(es_ID, cov = vcv)), # this is m
            sigma ~ 1 +cyear + se
)

# prior

prior6 <- default_prior(form6, 
                        data = dat, 
                        data2 = list(vcv = vcv),
                        family = gaussian()
)

prior6$prior[6] = "constant(1)"
prior6

# fit model

fit6 <- brm(form6, 
            data = dat, 
            data2 = list(vcv = vcv),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior6,
            #threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(fit6)

# save fit6 as rds

saveRDS(fit6, here("Rdata", "fit6.rds"))



