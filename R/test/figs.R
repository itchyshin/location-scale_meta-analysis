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
               metafor,
               tidybayes)

# functions# Function to get variable names dynamically
get_variables_dynamic <- function(model, pattern) {
  variables <- get_variables(model)
  variables[grep(pattern, variables)]
}

rename_vars <- function(variable) {
  variable <- gsub("b_Intercept", "b_l_intercept", variable)
  variable <- gsub("b_sigma_Intercept", "b_s_intercept", variable)  
  variable <- gsub("b_habitatterrestrial", "b_l_contrast", variable)            
  variable <- gsub("b_sigma_habitatterrestrial", "b_s_contrast", variable)
  variable <- gsub("b_methodpersisitent", "b_l_contrast", variable)            
  variable <- gsub("b_sigma_methodpersisitent", "b_s_contrast", variable)
  variable <- gsub("b_elevation_log", "b_l_elevation", variable)            
  variable <- gsub("b_sigma_elevation_log", "b_s_elevation", variable)
  variable <- gsub("b_se", "b_l_se", variable) # 
  variable <- gsub("b_sigma_se", "b_s_se", variable) # 
  variable <- gsub("b_cyear", "b_l_cyear", variable) # 
  variable <- gsub("b_sigma_cyear", "b_s_cyear", variable) #
  variable <- gsub("sd_study_ID__Intercept", "sd_l_study_ID", variable)            
  variable <- gsub("sd_study_ID__sigma_Intercept", "sd_s_study_ID", variable)
  variable <- gsub("cor_study_ID__Intercept__sigma_Intercept", "cor_study_ID", variable)  
  return(variable)
}

# Function to visualize fixed effects
visualize_fixed_effects <- function(model, title) {
  fixed_effect_vars <- get_variables_dynamic(model, "^b_")
  if (length(fixed_effect_vars) == 0) {
    message("No fixed effects found")
    return(NULL)
  }
  
  tryCatch({
    fixed_effects_samples <- model %>%
      spread_draws(!!!syms(fixed_effect_vars)) %>%
      pivot_longer(cols = all_of(fixed_effect_vars), names_to = ".variable", values_to = ".value") %>%
      mutate(.variable = rename_vars(.variable))
    
    ggplot(fixed_effects_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        point_interval = "mean_qi", 
        fill = "lightcyan3", 
        color = "lightcyan4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(title = title, y = "", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_fixed_effects: ", e$message)
    return(NULL)
  })
}


###########
# Example 1
###########

dat <- read.csv(here("data", "thermal.csv"))

# selecting varaibles we need
dat <- dat %>% select(dARR, Var_dARR, es_ID, population_ID, study_ID, exp_design, habitat)

# creating SE (sqrt(Var))
dat$si <- sqrt(dat$Var_dARR)

# creating a new variable (method) according to the oringal paper
dat$method <- ifelse(dat$exp_design == c("A", "B", "C"), "initial", "persisitent")

vcv <- diag(dat$Var_dARR)
rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID


# 4 plots

# plot 1

fit_mr1 <- readRDS(here("Rdata", "fit_mr1.rds"))
p1 <- visualize_fixed_effects(fit_mr1, title = "a. Categorical biological")



# plot 2

# using metafor and use heteroscad model
# see - Nakagawa S, Lagisz M, O'Dea RE, Pottier P, Rutkowska J, Senior AM, Yang Y, Noble DW. orchaRd 2.0: An R package for visualising meta-analyses with orchard plots.
mr1_mod <- rma.mv(yi = dARR, 
                  V =  vcv,
                  mod = ~ habitat,
                  random = list(~1 | study_ID,
                                ~habitat | es_ID),
                  struct = "DIAG",
                  data = dat, 
                  test = "t",
                  sparse = TRUE,
                  control=list(optimizer="optim", optmethod="Nelder-Mead")
)

p2 <- orchard_plot(mr1_mod, mod = "habitat", 
             group = "study_ID", 
             xlab = "Effect size",
             trunk.size = 0.5,
             branch.size = 4,
             legend.pos = "top.left",
             cb = FALSE
             #legend.pos = "none"
             ) + 
  scale_colour_manual(values = rep("grey20",2))

p1 + p2


# plot 3

fit_mr2 <- readRDS(here("Rdata", "fit_mr2.rds"))

p3 <- visualize_fixed_effects(fit_mr2, title = "b. Categorical methodological")

# plot 4

mr2_mod <- rma.mv(yi = dARR, 
                  V =  vcv,
                  mod = ~ method,
                  random = list(~1 | study_ID,
                                ~method | es_ID),
                  struct = "DIAG",
                  data = dat, 
                  test = "t",
                  sparse = TRUE,
                  control=list(optimizer="optim", optmethod="Nelder-Mead")
)

p4 <- orchard_plot(mr2_mod, mod = "method", 
             group = "study_ID", 
             xlab = "Effect size",
             trunk.size = 0.5,
             branch.size = 4,
             legend.pos = "top.left",
             cb = FALSE
             #legend.pos = "none"
             ) + 
  scale_colour_manual(values = rep("grey20",2))

p3 + p4


############
# Example 2
############

dat <- read.csv(here("data", "elevation.csv"))
dim(dat)
# filter data LMA

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


# 2 plots

# plot 5

fit_mr3 <- readRDS(here("Rdata", "fit_mr3.rds"))

p5 <- visualize_fixed_effects(fit_mr3, title = "c. Continuous biological")

# plot 6

# using metafor and orchaRd
# see - Nakagawa S, Lagisz M, O'Dea RE, Pottier P, Rutkowska J, Senior AM, Yang Y, Noble DW. orchaRd 2.0: An R package for visualising meta-analyses with orchard plots.
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
            xlab = "ln(elevation) (moderator)",
            ylab = "SDM (effect size)",
            legend.pos = "bottom.left")

p6 <- bubble_plot(mr3_mod, mod = "elevation_log", 
             group = "study_ID", 
             g = TRUE, 
             xlab = "ln(elevation) (moderator)",
             ylab = "lnRR (effect size)",
             legend.pos = "bottom.left",
             k.pos = "top.left") + 
  geom_point(data = dat,
             aes(x = elevation_log, y = yi,
                 #fill = "green",
                 #colour = "green",
                 size = 1/sqrt(vi)), alpha = 0.3,
             #shape = 23
             ) #+
  # scale_color_discrete() + 
  # scale_shape_manual(values = 21) +
  # scale_fill_manual(values = "#117733") +
  # scale_colour_manual(values = "grey20")
                                   
p5 + p6

############
# Example 3
############

dat <- read.csv(here("data", "seed.csv"))

dat$study_ID <- as.factor(dat$study)
dat$es_ID <- as.factor(1:nrow(dat))
dat$se <- sqrt(dat$var.eff.size) # SE (sampling error SD)
dat$smd <- dat$eff.size # effect isze
dat$cyear <- scale(dat$study.year, center = TRUE, scale = FALSE) # centred year


# selecting varaibles we need
dat <- dat %>% select(smd, var.eff.size, es_ID, study_ID, se, cyear)

dat$vi <- dat$var.eff.size
vcv <- diag(dat$var.eff.size)
rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID
# 3 plots 

# plot 7

fit_mr4 <- readRDS(here("Rdata", "fit_mr4.rds"))

p7 <- visualize_fixed_effects(fit_mr4, title = "d.  Publication bias")


# plot 8

# using metafor and orchaRd

mr4_mod <- rma.mv(yi = smd, 
                  V =  vcv,
                  mod = ~ se + cyear,
                  random = list(~1 | study_ID,
                                ~1 | es_ID),
                  #struct = "DIAG",
                  data = dat, 
                  test = "t",
                  sparse = TRUE,
                  control=list(optimizer="optim", optmethod="Nelder-Mead")
)

p8 <- bubble_plot(mr4_mod, mod = "se", group = "study_ID", g = TRUE, 
            xlab = "Standard error (SE: moderator)",
            ylab = "SMD (effect size)",
            legend.pos = "bottom.left",
            k.pos = "top.right") + 
  geom_point(data = dat,
             aes(x = se, y =smd,
                 size = 1/sqrt(var.eff.size)), alpha = 0.3,
  ) 

# plot 9

p9 <- bubble_plot(mr4_mod, mod = "cyear", group = "study_ID", g = TRUE, 
              xlab = "Centered publicaiton year (moderator)",
              ylab = "SDM (effect size)",
              legend.pos = "bottom.left",
              k.pos = "bottom.right") + 
  geom_point(data = dat,
             aes(x = cyear, y = smd,
                 size = 1/sqrt(var.eff.size)), alpha = 0.3,
  ) 

p7 + (p8 / p9)

### arranging

layout<- "
AABB
CCDD
EEFF
GGHH
GGII
"

p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(design = layout)
