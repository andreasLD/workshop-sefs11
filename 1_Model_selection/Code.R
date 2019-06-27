#################################################################################
# 			Script for Model selection part for the SEFS 2019 workshop			#
# 				"Advanced analysis of Ecological data with R" 					#
#							by Ralf B. Schäfer, 27.6.2019						#
#     All material comes without any guarantee and you should always reflect    #
# 		yourself what you want to do and how reliable a method is				#
#################################################################################

# You have to first set a working directory i.e. a directory where we store all files
setwd("Your_directory")
# My working directory: 
# setwd("~/Arbeit/Vortraege/2019/SEFS/Workshop")
# you can set the working directory manually in RStudio: Use the Tools | Change Working Dir... menu 
# (Session | Set Working Directory on a mac).
# if you set the directory manually you can use the following function to simplify the identification of your path
file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function
 
#' ## Introduction
#' In this session, we will conduct a multiple linear regression analysis with a strong focus on model selection. We conduct the analysis for a real dataset on marine ostracods. Our research goal is explanation. We want to identify the variables that explain the diversity of marine ostracods. We assume sparsity, i.e. that the true ostracod diversity is controlled by a few environmental variables. We seek to identify the most important variables, assuming a linear relationship of the explanatory variables with the response variable (species richness after rarefaction).
#' 
#' ## Retrieving data
#' We load the data from an article (Yasuhara et al. 2012) that is publicly available in the online repository http://datadryad.org.
data_oc <- read.csv("http://datadryad.org/bitstream/handle/10255/dryad.39576/OstraMRegS400JB.txt?sequence=1", sep = "\t")

#' We have the following variables in our data set:  
#' 
#' * E100: species richness rarefied to 100 individuals,
#' * SR: species richness,
#' * MDS 1: multidimensional scaling first axis,
#' * MDS 2: multidimensional scaling second axis,
#' * DP: water depth (m),
#' * RC: region code,
#' * BT: bottom water temperature,
#' * SA: salinity,
#' * SP: seasonality of productivity,
#' * IC: number of ice-free days per year,
#' * P: surface productivity,
#' * DCA 1: detrended correspondence analysis first axis,
#' * DCA 2: detrended correspondence analysis second axis,
#' * LAT: latitude,
#' * LON: longitude.
#' 
#' We remove the community information that is captured in the variables *DCA* and *MDS*, i.e. remove the variables *DCA 1*, *DCA 2*, *MDS 1* and *MDS 2* and store the new dataset without these information in the object `data_oc2`.
data_oc2 <- data_oc[ , !names(data_oc) %in% c("MDS1", "MDS2", "DCA1", "DCA2")]

# Usually, we would run some sort of exploratory analysis before modelling, however, we have limited time therefore we just conduct some transformation and refer the reader to the online tutorial, where details are provided: http://139.14.20.252:3838/session/5/#section-data-pre-processing
# In preparation for the modelling, we add the log-transformed variable DP to the data frame.
## ----add_variable, include = TRUE, echo = TRUE, context = "setup"--------
data_oc2$DPlog <- log10(data_oc$DP)
# In addition, we create a data set data_env in which we store all predictors except for potential responses and collinear variables.
data_env <- data_oc2[ , !names(data_oc2) %in% c("E100", "DP", "SR", "IC","SA")]

#' ## Fitting the full model
#' Before we fit the model, let us check the sample size. The function `describe()` provides several descriptive statistics for a dataset, including the number of observations and parameters, missing values, number of distinct values per variable and many more. Check out the help for the function for more information.
## ----sample_size, include = TRUE, echo = TRUE----------------------------
library(Hmisc)
describe(data_env)

#' We have 7 predictors and 136 observations in the data frame. If we included all predictors in the linear regression model, we would have 8 parameters (don't forget to count the intercept!) in the model. Interaction terms may further increase the number of parameters. We will include 1 interaction term in the model, i.e. will have 9 parameters. As rule of thumb the number of parameters should not exceed the square root of the sample size. Our number of parameters is lower and we do not need to apply techniques to reduce the number of parameters in the model. Harrell (2015: pp.72-74) suggests that the linear model results are reliable as long as the number of parameters is equal or lower than the sample size divided by 15. This more restrictive criterion would also be met. 
#' 
#'   Now we turn to fitting the full model. According to the article (Yasuhara et al. 2012), no interactions need to be considered. For the sake of demonstration we include one interaction term here: between temperature and productivity, assuming that productivity has a stronger effect on the ostracod diversity at higher temperature (irrespective whether this makes sense).
## ----full_mod, include = TRUE, echo = TRUE, context = "setup"------------
mod_1 <- lm(data_oc2$E100 ~ . + BT:SP, data = data_env, na.action = "na.fail")
summary(mod_1)

#' Most variables exhibit small *p*-values, except for LAT and the interaction term, where the *p*-values are around 0.5 and 0.7. Note that the intercept also has a *p*-value around 0.5.
#' Let us explore what happens, if we remove the `(Intercept)` from the model. We have used the `update()` function before that comes very handy to change models. `~.` means that everything remains the same, `-1` on the right hand side means that the intercept is removed.
## ----full_mod_wo_inter, include = TRUE, echo = TRUE----------------------
mod_1_wointer <- update(mod_1, . ~ . -1)
summary(mod_1_wointer)

#' Generally, you should not remove the `(Intercept)` from a linear regression model unless the true regression line is known to go through the origin and the data is in agreement with this knowledge, i.e. not forced through the origin. Be aware that the *R*^2^ is not comparable between models with and without intercept and can be high despite a poor model fit in no-intercept models (see Kvalseth 1985 for details).
#'   
#' ## Modelling approaches: Intro
#' In the following we will proceed with modelling. If our goal was prediction, we could just use the full model, unless a variable is clearly irrelevant for prediction. However, our research question is: *Which environmental variables control the diversity of marine ostracods?* and we seek to identify the most important variables, assuming a linear relationship of the
#' predictors with the response variable and sparsity, i.e. that only a few variables are relevant.  
#'   
#'   Ideally, we would construct a few alternative candidate models *a priori* based on previous studies and then select the best-fit model based on a goodness-of-fit (GOF) measure. However, the authors of the original studies did not follow this strategy and given my limited knowledge on marine ostracods (and that the aim is to become acquainted with several modelling strategies), we do not implement this strategy. In the following we go through several modelling strategies, starting with best subset selection.
#' 
#' ## Best subset selection
#' We calculate all possible models starting with *mod_1* as the most complex model. The `dredge()` function in the package [MuMin](https://cran.r-project.org/web/packages/MuMIn/index.html) computes all models and evaluates them in terms of a selected goodness of fit measure (defaults to AIC~c~).
#' 
## ----all_mods, include = TRUE, echo = TRUE-------------------------------
library(MuMIn)
allmodels <- dredge(mod_1, extra = "R^2")
print(allmodels)

#' The best fit model is given in the first line (note that the output is broken into two blocks; the first line continues in the middle of the output) and contains all variables except LAT and the interaction term BT:SP. It has an AIC~c~ of 747. Several other models have a similarly low AIC~c~ and include a similar number of variables. We could extract the top model  (and the code to this is given below) for further analysis, explanation or prediction, but the fact that it is only slightly better than other models with plausible combinations of the predictors cautions against over-interpretation. In fact, the single top model probably tends toward overfitting and has a higher prediction variance. This issue may be ameliorated through model averaging.
## ----extr_top_mod, include = TRUE, echo = TRUE---------------------------
topmod <- get.models(allmodels, subset = 1)
print(topmod)
#' 
#' ### Model averaging for several best-fit models
#' Model averaging has been criticized and should not be used without further reading on the topic (e.g. Cade 2015, Banner & Higgs 2016, Dormann et al. 2018, Fieberg & Johnson 2015). A key question is how much weight the different models obtain in the averaging process. A study by Dormann et al. (2018) recently compared different methods to determine weights and concluded that a general recommendation for a method is difficult. Nevertheless, this study provides the R code to many methods of model averaging and may give you some guidance for selecting (and implementing) a method of weighting. Here, we average all models that exhibit a maximum difference of 2 in the AIC~c~ to the best-fit model (i.e. $\Delta \text {AIC}_\text{c} < 2$). The main reason for doing this is to have a shorter output, generally you don't need to select a cutoff and see Murtaugh (2014) for a discussion of the cutoffs. 
#' 
## ----model_avg_mod, include = TRUE, echo = TRUE--------------------------
library(MuMIn)
avg_model <- model.avg(allmodels, subset = delta < 2, fit = TRUE)
summary(avg_model)

#' The output provides the component models that were used in the averaging (`Component Models:`), the averaged model coefficients (`Model-averaged coefficients:`), and a variable importance measure (`Relative variable importance:`), which is simply counting the number of times a variable is included in the component models. The averaged model coefficients are given as *full average* and *conditional average*. The *full average* uses a regression coefficient of 0 in averaging if a variable is not contained in a component model, whereas the conditional average only averages those coefficients from component models where the coefficient is unequal 0, i.e. the variable is included. Thus, the former method includes the shrinking of coefficients based on the frequency of their occurrence in the component models.  
#'   
#' Finally, we rerun the `dredge()` function, but include the argument `rank = "BIC"` for using the BIC as GOF and assign the output to a new object. Subsequently, we run the `model.avg()` function with this object and assign the output to another object, for which we call the `summary()` function afterwards. 
#' 
## ----recalc_models, include = TRUE, echo = TRUE--------------------------
allmodels2 <- dredge(mod_1, extra = "R^2", rank = "BIC")
avg_model_bic <- model.avg(allmodels2, subset = delta < 2, fit = TRUE)
summary(avg_model_bic)

#' Finally, we compute the *R*^2^ as a measure of model accuracy. As a first step we compute the fitted values using the `predict()` function, subsequently we calculate the residual sum of squares. 

## ----model_avg_accuracy_meas, include = TRUE, echo = TRUE----------------
fit_y <-  predict(avg_model)
res_ssq_avg <- sum((data_oc2$E100 - fit_y) ^ 2)
tot_ssq_avg <- sum((data_oc2$E100 - mean(data_oc2$E100)) ^ 2)
# Based on the residual sum of squares and the total sum of squares we can compute the R^2:
1 - res_ssq_avg / tot_ssq_avg
 
#' For calculation of the cross-validated root mean square prediction error (RMSPE), sometimes just termed RMSE (see for example Dormann et al. 2018), use the arguments `rank = loo` and `type = "rmse"` in the function *dredge*. This runs a leave-one-out (LOO) cross-validation. Other forms of cross-validation (e.g. 5-fold cross-validation) are currently not implemented but could be conducted manually by splitting the data in five random sets before analysis.
#' 
#' ## Stepwise model selection
#' 
#' ### Hypothesis-based stepwise model selection
#' We start with stepwise backward model selection based on the assessment of hypotheses. This requires the definition of a fixed cutoff value for the *p*-value, when to remove a variable from the model and in turn comes with all the problems of dichotomising the *p*-value scale (e.g. see my lecture notes from Session 3 of my course on data analysis https://github.com/rbslandau/Data_analysis).
#' 

#' In backward model selection, we start with the full model and remove variables. We have already inspected the regression output for the full model *mod_1* (see section *Fitting the full model*). We remove the interaction term in the first step based on its high *p*-value.
## ----stepwise_mod1, include = TRUE, echo = TRUE, context = "setup"-------
summary(mod_1)
# Remove interaction
mod_2 <- update(mod_1, ~. - BT:SP)
# Check model
summary(mod_2)

#' Note that by removing one term, the estimates and the test statistics for the other terms have changed. Actually, the estimates and standard errors are most reliable for the full model. Also note that by removing the interaction term, the *R*^2^ decreased but the *adjusted* *R*^2^ increased.  
#'   
#'   In the next step, we should remove the variable `LAT` given its high *p*-value. Note that the *p*-value for the *t*-test for `LAT` is equivalent to that from the ANOVA with the partial *F*-test, which compares the models with and without the variable. In addition, the test statistic for the regression output is identical to that of a Type 2 ANOVA (in the absence of the interaction). Finally, the function `drop1()` provides different tests for the changes in the model, in particular for the Generalized Linear Model (GLM). If called with the argument `test = "F"`, again, the output for the test statistic is identical.
## ----stepwise_mod2, include = TRUE, echo = TRUE--------------------------
# Regression output
summary(mod_2)
# Anova Type 2 output
library(car)
Anova(mod_2, type = 2)
# Output of partial F-test
mod_3 <- update(mod_2, ~. - LAT)
anova(mod_2, mod_3)
# Output of drop1
drop1(mod_2, test = "F")

#' You may ask yourself: *Why do we have all these different functions if they provide a similar output?*  
#'   
#'   Well, most functions can also be used in other contexts: `summary.lm()` can be used in the context of contrasts, `Anova()` can be used for Type 3 ANOVA, `drop1()` provides different test statistics and information-theoretic GOFs and `anova()` allows to compare nested models that differ by more than one variable. For example, we can compare *mod_1* and *mod_3* with `anova()`:
## ----stepwise_mod_anova, include = TRUE, echo = TRUE---------------------
anova(mod_1, mod_3)

#' Beside the modelling strategy to remove variables until all variables have *p*-values larger than a specified cutoff, you could also use a different strategy: remove variables until the *F*-test for the comparison of the explained variances between the simplified and the full model results in a *p*-value below a specified cutoff. 
#' 
#' **Now it is your turn.** Continue with manual backward model selection until all variables have a p-value < 0.05. Start with *mod_2*. 

#' 
#' ### Information-theoretic based stepward selection
#' As discussed, another option is stepwise selection based on information-theoretic criteria. Note that this approach also suffers from the issue of selective inference. However, we do not need to specify a cutoff value as for the *p*-value. Instead, we look for the minimum value of an information-theoretic criteria. The functions to calculate the AIC, BIC and corrected AIC are:
## ----stepwise_inform_theoretic, include = TRUE, echo = TRUE--------------
# Calculation of AIC
AIC(mod_1)
AIC(mod_2)
# Calculation of BIC
BIC(mod_1)
BIC(mod_2)
# Calculation of corrected AIC
library(MuMIn)
AICc(mod_1)
AICc(mod_2)

#' We see that the BIC provides the strongest reduction when removing a variable.
#' 
#' **Again, it is your turn.** Compare the different information-theoretic criteria for the models mod_1 to mod_5. (You need to (re)fit some of the models).

#' Note that based on the AIC, *mod_3* would be the best-fit model, whereas based on the BIC we should select the more parsimonious model *mod_4*. Given that our sample size is much larger than 50, it is no surprise that the corrected AIC is very similar to the normal AIC.
#' 

#' 
#' ## Automatic stepwise modelling
#' In this section, I provide you with the code to automatically conduct model selection based on information-theoretic criteria. We have discussed in the lecture that automated modelling should be used with extreme caution only. Notwithstanding, if specified correctly, the `step()` function yields to the same result as manual stepwise selection, but is obviously faster (and requires less lines of code). The `step()` function requires the following inputs:  
#'   
#'   * object: Initial model used to start the stepwise selection  
#' 
#' * scope: defines the sparsest ("lower") and richest ("upper") model that is acceptable
#' 
#' * direction: specifies the search direction, i.e. stepwise "forward" checks whether a variable should be added, stepwise "backward" whether a variable should be removed and "both" combines these.  
#' 
#' * trace: defines the level of information provided during the search process  
#' 
#' * k: defaults to 2 for AIC; if set to *log(n)* provides the BIC (though the output will still state "AIC")
#' 
#' For details and alternative options refer to the help pages (`?step`).   
#'   
#'   Typically, the richest model that is acceptable is the full model and the sparsest model acceptable is the intercept-only model, also called *null model*. We have fitted the full model above (*mod_1*). In the following, we first fit the *null model* and subsequently start the automatic modelling process, using the BIC as GOF measure.
## ----auto_modelling, include = TRUE, echo = TRUE-------------------------
# fit intercept-only nullmodel: no variables, only mean
nullmodel <- lm(data_oc2$E100 ~ 1, data = data_env)
# sample size
n <- nrow(data_env)
# start stepwise algorithm
step(object = mod_1, scope = list(upper = mod_1, lower = nullmodel), direction = "backward", trace = 100, k = log(n))

#' 
#' We see that the automatic backward model building checks for each step, how the BIC (output still states "AIC") changes if a specific variable or any variable is removed `<none>`. If `<none>` exhibits the lowest value, the modelling is stopped. The automated modelling with the BIC yields to the same model as our manual modelling above.  
#'   
#'   Let us determine the best-fit model if we use a forward search strategy and start with the null model. 
## ----auto_modelling2, include = TRUE, echo = TRUE------------------------
step(
  nullmodel, direction = "forward", trace = 100,
  scope = list(upper = mod_1, lower = nullmodel), k = log(n)
)

#' We end up with the same model. This increases our confidence in the modelling result. However, the outcome of the stepwise modelling algorithm is often sensitive to the selected GOF measure and the direction of the search process, i.e. different choices lead to different best-fit models. In this context, note that model averaging can represent an attractive response to deal with the uncertainty generated if the results of stepwise modelling, irrespective of whether manual or automated, are highly sensitive to the selected GOF measure and modelling direction and yield to multiple plausible models. Finally, remember that in the lecture we discussed that backward model selection should be the method of choice in most cases and the BIC if assuming a very sparse true model.
#' 
#' ### Bootstrapping of stepwise modelling 
#' A final method we run is bootstrapping of stepwise modelling. Note that Austin (2008) cautioned against the use of bootstrapping to fix the problems of stepwise model selection and see Harrell (2015:70-71) for an overview of related problems. Nevertheless, it may inform on the variability of our selection results and on the consistency in variable selection (same rationale as stability selection for the LASSO, see below).
#'   
#'   First, we run the method using the BIC as GOF.
## ----bootstrapping_BIC, include = TRUE, echo = TRUE----------------------
library(bootStepAIC)
# set seed to make analysis reproducible
set.seed(111)
# See help for details on function
# as before k can be used to select BIC
boot.stepAIC(object = mod_1, data = data_env, k = log(n), direction = "backward")

#' The function provides the following output from bootstrapping: `Covariates selected`, `Coefficients Sign` and `Stat Significance`. `Covariates selected` provides information how frequent a variable was selected for the final model in stepwise selection for the bootstrap samples. The other two elements of the output inform on the variability in the regression coefficients in terms of direction and their statistical significance at $\alpha = 0.05$. Below these elements is the output from the original (non-bootstrap) stepwise procedure. Overall, the bootstrapping results show that only few variables were selected and for most bootstrap samples the final model was the null model. However, the variables that were removed in the original procedure based on the BIC such as *RC* and *LAT* were indeed less frequently selected in the final bootstrap models than variables such as *BT*, *SP* and *DPlog*. The bootstrapping results also suggest that the selection of *LON* and *P* in the original stepwise modelling should be interpreted with caution because they were less frequently selected in the final bootstrapped models than variables such as *LAT* that, in turn, were removed in the original stepwise modelling. Conversely, the interaction term was as frequently selected as was *DPlog*, but not considered in original stepwise modelling. The changes in the coefficients of *BT*, *LON*, *DPlog* and *BT:SP* also caution against overinterpreting the direction of the regression coefficients. Overall, even if the bootstrapping can not fix the problems of stepwise selection, it can inform how robust the results of a manual or automated selection are.  
#'   
#'   **Rerun the bootstrapping procedure for the AIC, i.e. by modification of the argument `k = ...`.**

#' If you accomplished this minor exercise, you will see that for the AIC the picture is similar. *BT* and *SP* are most frequently selected, all other variables have been selected with a relatively similar frequency (i.e. 17-23%), except for *LAT* (13%) and *RC* (29%). This provides some confidence that *BT* and *SP* are indeed very important variables and that *LAT* is not very relevant. However, the removal of *RC* in the original stepwise procedure based on the AIC should be interpreted with caution and the remaining variables seem similarly important. Finally, the variability in the signs of the regression coefficients calls for caution when interpreting their direction.

#' ## Model diagnosis
#' Of course, you still need to diagnose the model assumptions for the model that results from model selection. Indeed, if you find serious violations of the model assumptions than the whole modelling process is likely not reliable. However, model diagnosis is beyond the short time we have in this workshop. Again, you will find guidance on model diagnosis in my lecture or in text books.
#' 
#' ## Relative importance of a variable
#' After model diagnosis, we may be interested in the relative importance of the variables in the model. Hierarchical partitioning and the method *PMVD* are most suitable, hence, we focus on these but compute *standardized betas* for comparison. We use the package [relaimpo](https://cran.r-project.org/package=relaimpo), note that due to patent issues the measure *PMVD* is not included in the package but [can be installed manually from the website of the maintainer](https://prof.beuth-hochschule.de/groemping/software/relaimpo/). We are not going into details of the different measures that are available, but merely focus on the demonstration in R. Check out Grömping (2015), the help for the function `calc.relimp()` and the related references for further information on relative importance measures. The hierarchical partitioning is called using the argument `type = c("lmg")`. This is because *lmg* is a special case of hierarchical partitioning for linear models based on the *R*^2^. Let us calculate both *lmg* and *standardized betas* for comparison.
#' 
## ----rela_impo, include = TRUE, echo = TRUE------------------------------
library(relaimpo)
pred_imp_lmg <- calc.relimp(mod_2, type = c("lmg"), rela = TRUE)
pred_imp_beta <- calc.relimp(mod_2, type = c("betasq"), rela = TRUE)
plot(pred_imp_lmg, main = "")
plot(pred_imp_beta, main = "")

#'   
#'   Note that you can actually compute different measures with one line of code by providing the type argument with multiple measures, for example: `type = c("lmg", "betasq")`. I did not do this in this tutorial for layout reasons. For both measures, the variables *BT*, *P* and *DPlog* are most important. In hierarchical partitioning *SP* is less relevant than when the assessment is based on its standardized regression coefficient. The standardized regression coefficients have been discussed in the lecture, the `calc.relimp()` function provides them squared to simplify their comparison.  

#' ## Shrinkage
#' As discussed, stepwise model selection is associated with several problems. 

#' One issue is that the regression coefficients are biased high. Post-selection shrinkage may provide a partial fix and we run post-selection shrinkage in R to obtain more accurate regression coefficients.
#' We need to refit the model with `lm()` and set additional arguments (`x = TRUE, y = TRUE`), because the shrinkage function requires the related information. Let us fit the model with the lowest BIC.
## ----shrink_prep, include = TRUE, echo = TRUE----------------------------
mod_5_s <- lm(data_oc2$E100 ~ BT + SP + P + LON + DPlog, data = data_env, x = TRUE, y = TRUE)

#' Now we run global shrinkage, which shrinks all parameters with the same factor, and parameterwise shrinkage. Parameterwise shrinkage combined with backward selection resulted in
#' the most accurate result in a simulation study (Houwelingen & Sauerbrei 2013). If you have categorical variables in the model use "joint" shrinkage (see Dunkler et al. 2016 for details).
#' 
## ----shrink_conduc, include = TRUE, echo = TRUE--------------------------
library(shrink)
# global shrinkage
shrink_res1 <- shrink(mod_5_s, type = "global")
shrink_res1
# reproduce results manually
coef(mod_5_s)[-1] * shrink_res1$ShrinkageFactors
# note that the intercept is removed because the intercept requires no shrinkage

# parameterwise shrinkage
shrink_res2 <- shrink(mod_5_s, type = "parameterwise")
shrink_res2

#' For our case, the results for global and parameterwise shrinkage are largely the same, the differences are for most estimated regression coefficients below 10-15%.  
#'   
#'   Bootstrapping and cross-validation of the whole stepwise selection procedure have been suggested to reduce overfitting and represent an alternative to post-selection shrinkage. The code to conduct bootstrapping has been provided above. We do not look into the details of cross-validation for stepwise regression, but the R implementation has been featured in a [blog post](https://www.r-bloggers.com/variable-selection-using-cross-validation-and-other-techniques/).
#' 
#' ## The LASSO
#' The LASSO represents an alternative to the stepwise model selection and post-selection shrinkage discussed above. It provides simultaneous variable selection and shrinkage. 

#' A detailed treatment of the LASSO (and other related methods such as ridge regression, least angle regression and the elastic net (see Hastie et al. 2015)) is beyond the level of this course, but we run an analysis for our data set. The package [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) provides the functions to run the LASSO, but also for other methods such as the elastic net. The package comes with a vignette that provides detailed information on model fitting and interpretation and I suggest to study [this vignette](https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet_beta.pdf), in case you want to fit a LASSO model to your data. Enough said/written, let's run the LASSO! The main function to do this is `glmnet()`. Note that to interpret the relative importance of the regression coefficients, you should standardize the predictors before the analysis. The data are typically centred and standardized before the LASSO is applied. Indeed, the function `glmnet()` performs these operations behind the scenes, but converts the estimated coefficients back into the original units in reporting (e.g. plots and extraction of coefficients). Hence, you can not evaluate the variable importance based on the resulting plots and coefficients (this statement relates to the evaluation of the relative importance of variables with non-zero regression coefficients. However, you can conclude that a variable that is shrunk to zero is less important than a variable that is not shrunk to zero for a given penalty), unless you provide centred and standardized data to the function.
## ----lasso_1, include = TRUE, echo = TRUE--------------------------------
library(glmnet)
# fit model with lasso, requires predictors as matrix
# and response as vector
lasso_mod <- glmnet(x = as.matrix(data_env), y = data_oc2$E100)
plot(lasso_mod, label = TRUE)

#'   
#'   The figure shows the path of regression coefficients as the $\ell_1$-Norm (equivalent to the absolute sum of regression coefficients) increases along the x-axis. The numbers on the top of the figure tell us how many regression coefficients are unequal zero. We can also plot the coefficients against the penalty term $\lambda$ or against the explained variance.
## ----lasso_2, include = TRUE, echo = TRUE--------------------------------
plot(lasso_mod, label = TRUE, xvar = "lambda")
plot(lasso_mod, label = TRUE, xvar = "dev")

#'   
#'   The interpretation of these plots with respect to the variables is largely the same. The variables 2 and 7 (numbers refer to column numbers in the data matrix, representing the variables *BT* and *DPlog*) are among the three most important ones, as their estimated betas are unequal to zero for high penalty (e.g. log $\lambda$ of 0). It is not easily discernible in the plot, which third variable is unequal zero at a log $\lambda$ of 0. However, you could identify this third variable based on the regression coefficients (see below). Anyway, the crucial question is: *Which penalty should we choose?*  
#'   A hypothesis testing approach to the LASSO is provided in Taylor & Tibshirani (2018). However, we determine the optimal $\lambda$ based on cross-validation using the function `cv.glmnet()`.
## ----lasso_cv, include = TRUE, echo = TRUE-------------------------------
# set seed to make reproducible
set.seed(111)
cvfit <- cv.glmnet(as.matrix(data_env), data_oc2$E100)
plot(cvfit)
 
#'   
#'   The figure displays the cross-validated Mean squared error (MSE) for different penalties. On top of the figure, we see again the number of variables in the models. The dashed lines indicate the $\lambda$s corresponding to the model with the lowest cross-validation MSE and to the model that is within one standard error from the minimum model. In cross-validation, often this *one-standard error rule* is applied. This rule states that the most parsimonious model that is within one standard error of the best-fit model should be selected. The rationale is that models within one standard error appear equally good and given a tendency towards overfitting, the more simple model should be selected (see James et al. 2017: 214). Anyway, we do not need to read these $\lambda$s from the plot, we can extract them from the fitted object and also use them to extract the related regression coefficients. 
## ----lasso_extract, include = TRUE, echo = TRUE--------------------------
# extract lambdas
cvfit$lambda.min
cvfit$lambda.1se
# extract regression coefficients
coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.1se")

#' We see that except for *LAT* all variables are contained in the model that has been selected using the one-standard error rule. The minimum model even includes all variables. In the section on bootstrapping of the stepwise model, we will see that indeed the removal of variables, except for *LAT*, is not as obvious as the automated stepwise modelling suggests.  
#'   
#'   As stated above, you could also use the `coef()` function to identify the third variable that has a non-zero regression coefficient at a log $\lambda$ of 0. You need to provide the value for $\lambda$ (not for log $\lambda$) to the function, as argument for `s`. Hence, to extract the coefficients for a log $\lambda$ of 0 you need to provide the argument `s = 1` (because $\lambda$ = 1 corresponds to the log $\lambda$ of 0).
## ----lasso_extract_zero, include = TRUE, echo = TRUE---------------------
coef(cvfit, s = 1)

#'   We can see that *P* is the third variable that is among the most important ones, and not *SP* (labelled *3* in plot) as one might erroneously conclude from the plots above.
#'   #' 
#' ## Extra section: Stability selection
#' We have mentioned several extensions to the LASSO . A detailed treatment of the adaptive LASSO, ridge regression and elastic net is beyond the scope, but you can find related *R* functions with examples in the following packages:  
#' 
#'   * Adaptive LASSO: function `adalasso()` (package [parcor](https://cran.r-project.org/web/packages/parcor/index.html))  
#'   
#'   * Elastic net: function `glmnet()` (package [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)) or function `enet()` (package [elasticnet](https://cran.r-project.org/web/packages/elasticnet/))  
#'   
#'   Remember that both ridge regression and the LASSO are just special cases of the elastic net, and therefore it should not surprise you that we can stick to the same function (i.e. `glmnet()`) to fit an elastic net model that we have already used for the LASSO. An example for using the elastic net is provided in [this open ebook](https://daviddalpiaz.github.io/r4sl/elastic-net.html), a comparison of ridge regression, LASSO and the elastic net along with their implementation in *R* is provided in [this blog post](https://www.datacamp.com/community/tutorials/tutorial-ridge-lasso-elastic-net).  
#'   
#'   Here, we focus on stability selection. Stability selection relies on resampling (e.g. bootstrapping) and estimates the probability that a variable is selected in a subsample related to the regularization region (i.e. value range of $\lambda$). A variable is considered more relevant, the higher its selection frequency (for details see Meinshausen & Bühlmann (2010)). We apply stability selection to our data set, for other case studies see Hothorn et al. (2011) and McArt et al. (2017). We use the function `stabsel()` from the package [stabs](https://cran.r-project.org/web/packages/stabs/). We need to set a threshold (1) for the cutoff (i.e. the selection % a variable needs to reach to be considered relevant) and (2) for the per-family error rate (see *Details* in the help of `stabsel()` for further information and Hofner et al. 2015). Typical cutoffs are between 0.6 and 0.9 (set via the argument *cutoff = *) and *PFER = 1* gives the tolerable expected number of falsely selected variables. The PFER controls the regularization region, i.e. the $\lambda$s that are used in the LASSO. Alternatively, the number of selected variables for a subsample *q* can be set directly to control the regularization region. How to use stability selection for the $\lambda$ corresponding to the minimum cross-validation error (see previous section), though not recommended, is shown in the example for the `stabsel()` function. For more information on these parameters see Meinshausen & Bühlmann (2010).
## ----lasso_stabsel, include = TRUE, echo = TRUE--------------------------
library(stabs)
## make reproducible
set.seed(1204)
(stab.lasso <- stabsel(x = as.matrix(data_env), y = data_oc2$E100,
                       fitfun = glmnet.lasso, cutoff = 0.70, PFER = 1))
# plot estimate for selection probability
plot(stab.lasso, main = "Lasso")
#'   
#'   The results suggest that only the variable *DPlog* shows a consistent relationship with the response, followed by *P* and *BT* that are slightly above a selection probability of 0.5. The result that *DPlog* is most relevant, followed by *P* and *BT*, matches with the results for the conventional LASSO model fitted before.
#'   

#' ## Model visualisation
#'  The problem with multiple regression is that often we cannot easily visualise the results. Remember that we visualised the regression with two predictors as a plane in a three-dimensional coordinate system. If we have more predictors in a model, this model defies visualisation. A naive approach would be to visualise the relationships through multiple scatter plots of the response against the different predictors, ignoring the multiple predictor context. This can lead to erroneous interpretation and conclusions, as highlighted below. The function `mcPlots` provides plots to compare the bivariate relationship between the response and a selected predictor to the relationship between the response and predictors when accounting for the other predictors. Without going into too much detail, the response is accounted for the other predictors by removing the variation explained by them, and similarly, the selected predictor is accounted for the other predictors by removing variation explained by them. Hence, the idea is to focus on the unique information in the response and the selected predictor that is not explained by or shared with other predictors. More formally, we use the residuals *e* from the model with the response *Y* explained by all other predictors, which are constituting the matrix ***Z***, as new response: *e*(*Y*|***Z***), and use the residuals from the model with the selected predictor *X* as response explained by the other predictors as new predictor: *e*(*X*|***Z***). For further details on visualisation refer to Fox and Weisberg (2019). We generate the plots for the predictor *SP*.
## ----mod_vis_pedagog, include = TRUE, echo = TRUE------------------------
library(car)
mcPlots(mod_2, ~SP, overlaid = FALSE)

#'   
#'   Note that the bivariate relationship is positive and relatively weak, whereas it is negative and stronger (indicated by a steeper slope) when accounting for the other variables. This example cautions against interpretation and pre-selecting variables for multiple regression analysis based on bivariate relationships.  
#'   
#'   *How should you visualise the results of a multiple linear regression model then?* A useful tool is the *predictor effect plot*. The main idea of this plot type is to fix all other predictors at their mean. Taking ${x}_{n}$ as our selected predictor, the regression equation becomes: $y_i = b_0 + b_1 \bar{x}_{1} + b_2 \bar{x}_{2} + ... + b_n {x}_{n}$. Thus, all other terms can be summed to a new intercept $b_{0, \text {new}}$ and our regression equation simplifies to $y_i = b_{0, \text {new}} + b_n {x}_{n}$. This relationship between the selected predictor and the response at the mean level of the other predictors can be displayed as a simple regression model. For further details I refer you to Fox and Weisberg (2019). You can create plots for your model using the code provided below that demonstrates the application for *mod_5_s* (see section *Shrinkage*).
## ----mod_vis_effects, include = TRUE, echo = TRUE------------------------
library(effects)
plot(predictorEffects(mod = mod_5_s, predictor = "SP"), ylab = "Rarefied richness (100 ind.)")
plot(predictorEffects(mod = mod_5_s, predictor = "BT"), ylab = "Rarefied richness (100 ind.)")
plot(predictorEffects(mod = mod_5_s, predictor = "DPlog"), ylab = "Rarefied richness (100 ind.)")

#'   
#'   You could create the predictor effect plots for all predictors at the same time by not specifying the `predictor` argument, i.e. for our case: `plot(predictorEffects(mod_5_s))`.
#' 
#' 
#' **References**  
#'   
#' Austin P.C. (2008) Bootstrap model selection had similar performance for selecting authentic and noise variables compared to backward variable elimination: a simulation study. *Journal of Clinical Epidemiology* 61, 1009-1017.  
#'
#' Banner K.M. & Higgs M.D. (2016) Considerations for assessing model averaging of regression coefficients. *Ecological Applications* 27, 78–93.  
#' 
#' Cade B.S. (2015) Model averaging and muddled multimodel inferences. *Ecology* 96, 2370–2382.  
#' 
#' Dormann C.F., Calabrese J.M., Guillera-Arroita G., Matechou E., Bahn V., Barton K., et al. (2018) Model averaging in ecology: a review of Bayesian, information-theoretic, and tactical approaches for predictive inference. *Ecological Monographs* 88, 485–504.  
#' 
#' Dunkler D., Sauerbrei W. & Heinze G. (2016) Global, Parameterwise and Joint Shrinkage Factor Estimation. *Journal of Statistical Software* 69, 1–19. 
#' [Free to download.](https://www.jstatsoft.org/index.php/jss/article/view/v069i08
#' 
#' Fieberg J. & Johnson D.H. (2015) MMI: Multimodel inference or models with management implications? *The Journal of Wildlife Management* 79, 708–718.  
#' 
#' Grömping, U. (2015). Variable importance in regression models. *WIREs Computational Statistics* 7, 137-152
#'
#' Harrell F.E. (2015) Regression modeling strategies: with applications to linear models, logistic regression, and survival analysis, Second edition. Springer, New York.
#'
#' Hastie T., Tibshirani R. & Wainwright M. (2015) Statistical learning with sparsity: the lasso and generalizations. CRC Press, Taylor & Francis Group, Boca Raton. [Free to download.](https://web.stanford.edu/~hastie/StatLearnSparsity_files/SLS_corrected_1.4.16.pdf)  
#'
#' Hofner B., Boccuto L. & Göker M. (2015). Controlling false discoveries in high-dimensional situations: boosting with stability selection. *BMC Bioinformatics* 16, 144. [Free to download.](https://doi.org/10.1186/s12859-015-0575-3) 
#'   
#' Hothorn T., Muller J., Schroder B., Kneib T. & Brandl R. (2011). Decomposing environmental, spatial, and spatiotemporal components of species distributions. *Ecological Monographs* 81, 329–347
#'  
#' Houwelingen H.C. van & Sauerbrei W. (2013) Cross-Validation, Shrinkage and Variable Selection in Linear Regression Revisited. *Open Journal of Statistics* 03, 79–102.
#' [Free to download.](https://file.scirp.org/pdf/OJS_2013042410543131.pdf)
#' 
#' James G., Witten D., Hastie T. & Tibshirani R. (2017) An introduction to statistical learning: with applications in R. Springer, New York.  
#'
#' Kvalseth T.O. (1985) Cautionary note about R2. *American Statistician* 39, 279–285
#'
#' McArt S.H., Urbanowicz C., McCoshum S., Irwin R.E. & Adler L.S. (2017). Landscape predictors of pathogen prevalence and range contractions in US bumblebees. *Proceedings of the Royal Society B: Biological Sciences* 284, 20172181. [Free to download (within university).]( https://royalsocietypublishing.org/doi/pdf/10.1098/rspb.2017.2181)  
#'   
#' Meinshausen N. & Bühlmann P. (2010). Stability selection. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 72, 417–473. [Free to download.](https://doi.org/10.1111/j.1467-9868.2010.00740.x)
#'
#' Murtaugh P.A. (2014) In defense of P values. *Ecology* 95, 611–617.
#' 
#' Taylor J. & Tibshirani R. (2018) Post-selection inference for -penalized likelihood models. *Canadian Journal of Statistics* 46, 41–61. [Free to download.](https://onlinelibrary.wiley.com/doi/full/10.1002/cjs.11313)
#' 
#' Yasuhara M, Hunt G, van Dijken G, Arrigo KR, Cronin TM, Wollenburg JE (2012) Patterns and controlling factors of species diversity in the Arctic Ocean.
#' *Journal of Biogeography* 39(11): 2081-2088. http://dx.doi.org/10.1111/j.1365-2699.2012.02758.x
#'  
################
# R exercise   #
################	

##############################
# Preparation for Exercise   #
##############################

# We load data that is contained within an R package
# If you cannot access the data, install the package
# via: install.packages("HSAUR2")

data("USairpollution", package = "HSAUR2")
# See package information for details on the data set
head(USairpollution)

#############################################################################
# Exercise: For an effective environmental protection,						#
# you need to know causes of pollution.										#
# In this case study (using real world data), the aim is to identify		#	
# the variables exhibiting the highest explanatory power					#
# for the SO2 air concentrations.											#
# Model the SO2 concentration as response and use the other variables		#
# as predictors. Compare the results for the following methods:				#
# 1) manual model building based on hypotheses,								#
# 2) automatic backward model selection with BIC							#
# 3) LASSO																	#
# Also compare the regression coefficients from post-selection shrinkage 	#	
# for 1) or 2) with those from the LASSO. Finally, conduct model diagnosis	#
# and plot the model with effect plots and determine the variable importance#
# for a final model that you select.										#
#############################################################################



