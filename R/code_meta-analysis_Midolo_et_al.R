###Global patterns of intraspecific leaf trait responses to elevation###
###Author: Gabriele Midolo (gabriele.midolo@natec.unibz.it)###
###R code for meta-regressions, analysis of publication bias and multi-model inference###

#see 'additional notes' at the end of the script

#Load data:
general_data<-read.csv("C:/Data/Midoloetal_meta-analysis_figshare_data/meta-analysis_Midolo_et_al.csv")
#please, see the description of each variable in the document available on figshare.com
table(general_data$trait) # number of observations per trait

dat<-general_data[general_data$trait=="SLA",] # select the trait you want, e.g. ("SLA")

#Code for meta-regressions####
#Calculate the variance-covariance matrix to be used in the meta-analysis:
#see Lajeunesse, M. J. (2011). Ecology, 92, 2049-2055.
#N.B. this code has been provided by Wolfgang Viechtbauer (see http://www.metafor-project.org/doku.php/analyses:gleser2009)
calc.v <- function(x) {
  v <- matrix((x$sd_control[1]^2 / (x$n_individuals[1] * x$control[1]^2)) , nrow=nrow(x), ncol=nrow(x))
  diag(v) <- x$vi
  v
} 
library(metafor)
V <- bldiag(lapply(split(dat, dat$common_id), calc.v)) # n.b. always calculate V for each dataset (trait) separately!
#N.B. the variance-covariance matrix can be estimated also via the metagear::covariance_commonControl() function.

###Run meta-regression models, examples####
#null model:
m0<-rma.mv(yi,V,data=dat,
           random =list(~ 1 | gradient_id/id,~1|species)) 
m0 # intercept is the mean pooled effect size, i.e. the overall effect of SLA variation across the data

m1<-rma.mv(yi~elevation_log,V,data=dat,
           random =list(~ 1 | gradient_id/id,~1|species))
m1 #model with moderator(s), e.g. log-transformed elevation
#N.B. random component structure is based on three-level meta-analysis approach in Konstantopoulos, S. (2011), Research Synthesis Methods, 2(1), 61–76. 


###Plot meta-analytic scatter plots, example####
#e.g. for m1:
dat$wi    <- 1/sqrt(dat$vi)
dat$size  <- .8+2*(dat$wi - min(dat$wi))/(max(dat$wi) - min(dat$wi))
preds <- predict(m1, newmods=c(log(min(dat$elevation):max(dat$elevation))))
head(preds)
plot(dat$elevation_log, dat$yi, axes=TRUE, xlab="",ylab="",type="n")
mtext("Elevation (m) [log]", side=1, line=2, cex=.85)
mtext("SLA lnRR", side=2, line=2, cex=.85)
polygon(c(log(min(dat$elevation):max(dat$elevation)), rev(log(min(dat$elevation):max(dat$elevation)))), c(preds$ci.ub, rev(preds$ci.lb)),
        col = "grey70", border = NA)
points(dat$elevation_log, dat$yi, cex=dat$size)
abline(a=0,b=0,lty="dashed")
lines(log(min(dat$elevation):max(dat$elevation)), preds$pred, lwd=3)


###Egger test and funnel plot of the null-model####
#These can be applied to other models, e.g. if the inclusion of moderators is expected to reduce asimmetry originated from true heterogeneity.
#See Nakagawa and Santos (2012) for a description of the Egger test used here: 
#https://doi.org/10.1007/s10682-012-9555-5
dat$residuals<-residuals.rma(m0)
dat$precision<-1/sqrt(dat$vi)
egger<-rma.mv(residuals~precision,vi,data=dat,random =list(~ 1 | gradient_id/id,~1|species))
egger # check p-value of the intercept term, i.e. the Egger's P
#significant P-val indicate presence of funnel plot asimmetry

#Visualize a funnel-plot:
funnel(egger, level=c(90, 95, 99), 
       shade=c("white", "gray", "darkgray"), #colors to indicate the pseudo confidence interval regions
       yaxis="seinv",#set the precision of the observation to be plotted on the y-axis 
       digits=c(2,2))
#see more at http://www.metafor-project.org/doku.php/plots:contour_enhanced_funnel_plot


###Multimodel inference analysis using 'glmulti' package ####
#see a nice example here: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti
library(glmulti)
setOldClass("rma.mv")
setMethod('getfit', 'rma.mv', function(object, ...) {
  if (object$test=="z") {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
  } else {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
  }
})

#Function for glmulti to handle multi-level mixed-effect meta-analytical models:
rma.glmulti.ran <- function (formula, data, random, ...) {
  rma.mv(as.formula(paste(deparse(formula))), V, random =list(~ 1 | gradient_id/id, ~1|species),data=dat, method="ML", ...)
}

#scale moderators and log-transform the ones with the positive skeweness
dat$dele<-scale(log(dat$elevation))
dat$ai<-scale(dat$ARIDITY_INDEX)
dat$mgst<-scale(dat$MEAN_GROWING_SEASON_TEMPERATURE)
dat$aele<-scale(dat$elevation_control)
dat$srad<-scale(log(dat$SOLAR_RADIATION))
dat$lat<-scale(dat$LAT)

dat_cor<-dat[,c("dele","ai","mgst","aele","srad","lat")]
corrplot::corrplot(cor(dat_cor),method="number") # check collinearity

res <- glmulti(yi~dele+ai+mgst+aele+srad+lat, 
               data=dat, level=1, 
               fitfunction=rma.glmulti.ran,
               crit="aicc") # rank models based on AICc
plot(res, type="s") # plot the relative importance of each moderator
print(res)
topSLA<-summary(res@objects[[1]]) #model with the lowest aicc
topSLA

output_SLA<-round(coef(res,
                       select = .95, # select models up to sum of AICc weight ≥ .95
                       varweighting = "Johnson",alphaIC=.05), #average the estimate for each moderator based on their weight
                  4) # round to 4 digits
output_SLA


### Additional Notes ####

#The dataset we used here already has lnRRs (yi) and the sampling variances (vi) calculated using metafor::escalc(), e.g.:
dat_SLA<-metafor::escalc(measure="ROM", m1i=Tr, m2i=Ct, sd1i=sdT, sd2i=sdC, n1i=n, n2i=n,data=dat_SLA)
#Where:
#Tr = mean of the trait sampled in the 'treatment' (i.e. above the site at the lowest elevation)
#Ct = mean of the trait sampled at the 'control' (i.e. the lowest elevation)
#sdT = sd of the mean sampled above the site at the lowest elevation 
#sdC = sd of the mean sampled at the lowest elevation
#n = number of individuals (sample size)(one can specify if treatment and control )
#The missing standard deviation (sd) values were obtained in each dataset via metagear::impute_sd(), e.g.:
dat_SLA<-metagear::impute_SD(dat_SLA, "sd", "mean", method = "Bracken1992")
#this requires the mean (and sd) of traits to be each on the same column or vector; then you can separate mean and sd of the control and the treatment to calculate lnRRs via escalc().

#See the 'Methods' section of the manuscript for additional details

#For any question, please contact me at gabriele.midolo@natec.unibz.it