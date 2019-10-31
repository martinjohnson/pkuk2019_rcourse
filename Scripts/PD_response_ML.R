library(kernlab)
library(caret)
library(tidyverse)
SEED = 123456
set.seed(SEED)
source("Scripts/helpers.R")
################ Simulations ######################
N <- 2000

simdat <- tibble(
  id = paste("id", 1:N, sep = "_"),
  bw = rnorm(N, 0, 2.5) + 80,
  sex = rbinom(N, 1, .5) %>% as.factor
)

pkpop_vc <- function(bw, sex) {
  vcref <- 3.6
  bwref <- 80
  vcbw <- .597
  vcsex <- .152
  vc <- vcref * (bw / bwref)^(vcbw) * exp(vcsex)^I(sex == "1") 
  vc %>% as.double()
}

simdat <- simdat %>% 
  mutate(vc = pkpop_vc(bw, sex),
         dose = bw * 10^-3,
         cmax = (dose / vc) ) ## cmax in ng/mL


pr_bernouilli <- function(cmax) {
  beta1 <- 2
  eta <- scale_this(cmax) * beta1 ## normalise to obtain reasonable range of probabilities.
  pr <- 1 / (1 + exp(eta))
  rbinom(N, 1, pr)
}

simdat <- simdat %>%
  mutate(events = pr_bernouilli(cmax))

simdat %>% group_by(events) %>% count()

### Step-back

fitdat <- simdat %>%
  mutate(Class = if_else(events == 1, "Fail", "Response")) %>%
  select(bw, sex, dose, Class)

### Two modellers are going to be tested on this an ML vs a PK approach.
## The goal is to predict the event class.
## First, we split the data into two groups: a training set and a test set. To do this, the createDataPartition function is used:

inTrain <- createDataPartition(
  y = fitdat$Class,
  ## the outcome data
  p = .7,
  ## percentage of data in training set
  list = FALSE
)

## createDataPartition does a stratified random split of the data.
training <- fitdat[ inTrain, ]
testing <- fitdat[-inTrain, ]

nrow(training)
nrow(testing)

## To fit a model in caret use the train function.
## For example a partial least squared discriminant analysis (PLSDA) model.

plsFit <- train(
  Class ~ .,
  data = training,
  method = "pls",
  ## Center and scale the predictors for the training
  ## set and all future samples.
  preProc = c("center", "scale")
)

plsFit

# Set conditions for training model and cross-validation: 
number <- 5
repeats <- 2
control <- trainControl(method = "repeatedcv", 
                        number = number , 
                        repeats = repeats, 
                        classProbs = TRUE, 
                        savePredictions = "final", 
                        index = createResample(training$Class, repeats*number), 
                        summaryFunction = multiClassSummary, 
                        allowParallel = TRUE)

### Fit using ML method
## A naive machine learning programmer considers that using the variables sex and age and a Suppot Vector Machine model is the best option for such a problem.

# Train SVM: 
svmFit <-  ksvm(Class ~ . , 
            data = training %>%
              mutate_if(is_double, scale_this),
            kernel = "rbfdot", 
            C = 10, 
            epsilon = 0.1, 
            prob.model = TRUE, 
            cross = 10)

## Train logistic model using PK knowledge
pk_training <- training %>%
  mutate(vc = pkpop_vc(bw, sex),
         cmax = (dose / vc) )

# Train Logistic Model: 

logistic <- train(factor(Class) ~ cmax, 
                  data = pk_training, 
                  method = "glm", 
                  preProcess = c("center", "scale"),
                  trControl = control)
coef(logistic$finalModel)
#-----------------------------------
#  Compare between the two models
#-----------------------------------

# Function calculates probabilities: 
predict_prob <- function(model_selected, modeller = c("ML", "PK")) {
    predict(model_selected, 
            testing %>%
              mutate(vc = pkpop_vc(bw, sex),
                     cmax = (dose / vc) ) %>%
              mutate_if(is_double, scale_this),
            type = "prob") %>% 
      as.data.frame() %>%
      pull(Fail) %>%
      return()
}

# Use this function: 
pred_logistic <- predict_prob(logistic, "PK")
pred_svm <- predict_prob(svmFit, "ML")

# Function calculates AUC: 
library(pROC) 

test_auc <- function(prob) {
  roc(testing$Class, prob) 
}

# Use this function: 
auc_logistic <- test_auc(pred_logistic)
auc_svm <- test_auc(pred_svm)

# Create a data frame for comparing: 

df_auc <- bind_rows(tibble(TPR = auc_logistic$sensitivities, 
                               FPR = 1 - auc_logistic$specificities, 
                               Model = "Logistic, PK modeller"), 
                    tibble(TPR = auc_svm$sensitivities, 
                               FPR = 1 - auc_svm$specificities, 
                               Model = "SVM, ML modeller"))

# Plot ROC curves: 
df_auc %>% 
  ggplot(aes(FPR, TPR, color = Model)) +
  geom_line(size = 1) +
  theme_bw() +
  coord_equal() +
  geom_abline(intercept = 0, slope = 1, color = "gray37", size = 1, linetype = "dashed") + 
  labs(x = "FPR (1 - Specificity)", 
       y = "TPR (Sensitivity)", 
       title = "ROC Curve and AUC: Logistic vs SVM")


# Compare AUC: 
lapply(list(auc_logistic, auc_svm), function(x) {x[["auc"]]})


# Results based on test data: 
lapply(list(logistic, svmFit), 
       function(model) {confusionMatrix(predict(model, testing %>% 
                                                  mutate(vc = pkpop_vc(bw, sex),
                                                         cmax = (dose / vc) ) %>%
                                                  mutate_if(is_double, scale_this)
                                                  ), testing$Class, positive = "Response")})


## How do SVM work?
N <- 1000

# Training
simdat <- tibble(
  target_1 = runif(N),
  target_2 = runif(N)
)

simdat <- simdat %>%
  mutate(y = if_else( target_2 > target_1, 1, -1) )



inTrain <- createDataPartition(
  y = simdat$y,
  ## the outcome data
  p = .7,
  ## percentage of data in training set
  list = FALSE
)

## createDataPartition does a stratified random split of the data.
training <- simdat[ inTrain, ]
testing <- simdat[-inTrain, ]

p1 <- training %>%
  mutate_at("y", as.factor) %>%
  ggplot(aes(x = target_1, y = target_2, colour = y)) +
  geom_point() + theme_bw()

p1 + ggtitle("Training data")

# Vanilldadot gives a linear SVM
fit1 <- ksvm(factor(y) ~ .,
             data = training,
             kernel = "rbfdot", 
             C = 10, 
             epsilon = 0.1, 
             prob.model = TRUE, 
             cross = 10)

### Plot location of support vectors
SV1 <- training[alphaindex(fit1)[[1]],]

SV1 %>%
  mutate_at("y", as.factor) %>%
  mutate_if(is_double, scale_this) %>%
  ggplot(aes(x = target_1, y = target_2, colour = y)) +
  geom_point()  + theme_bw()  + ggtitle("Support vectors")

# Generate data for the decision line
w <- colSums(fit1@xmatrix[[1]] * fit1@coef[[1]])
b <- b(fit1)

p2 <- testing %>%
  mutate_at("y", as.factor) %>%
  mutate_if(is_double, scale_this) %>%
  ggplot(aes(x = target_1, y = target_2, colour = y)) +
  geom_point() + theme_bw() +
  geom_abline(intercept = b/w[1], slope = -w[2]/w[1]) +
  geom_abline(intercept = (b+1)/w[1], slope = -w[2]/w[1], linetype = 2) +
  geom_abline(intercept = (b-1)/w[1], slope = -w[2]/w[1], linetype = 2) 

p2 + ggtitle("Decision boundary - test data")

### Adding genomic variables

################ Simulations ######################
N <- 100
n_gene <- 150
genesim <- matrix(rbinom(N*n_gene, 1, .7) , ncol = n_gene)
colnames(genesim) <- paste("gene_", 1:n_gene, sep = "_")


simdat <- tibble(
  id = paste("id", 1:N, sep = "_"),
  bw = rnorm(N, 0, 2.5) + 80,
  sex = rbinom(N, 1, .5) %>% as.factor
)

simdat <- bind_cols(simdat, genesim %>% as.data.frame())
simdat[1:6, 1:6]

pkpop_vc <- function(bw, sex) {
  vcref <- 3.6
  bwref <- 80
  vcbw <- .597
  vcsex <- .152
  vc <- vcref * (bw / bwref)^(vcbw) * exp(vcsex)^I(sex == "1") 
  vc %>% as.double()
}

simdat <- simdat %>% 
  mutate(vc = pkpop_vc(bw, sex),
         dose = bw * 10^-3,
         cmax = (dose / vc) ) ## cmax in ng/mL


pr_bernouilli_gene <- function(dat) {
  beta1 <- 2
  gene_marker <- paste("gene_", sample(1:n_gene, 6), sep = "_")
  beta_gene <- c(0.9, -0.2, 0.5, -2.2, 0.7,  1.5)
  
  eta <- scale_this(dat$cmax) * beta1 
  eta <- eta + (  as.matrix( dat %>% select(gene_marker))  %*% beta_gene  ) ## use R matrix multiplication
  pr <- 1 / (1 + exp(eta))
  rbinom(N, 1, pr)
}

simdat <- simdat %>%
  mutate(events = pr_bernouilli_gene(.))

simdat %>% group_by(events) %>% count()

fitdat_gene <- simdat %>%
  mutate(Class = if_else(events == 1, "Fail", "Response")) %>%
  select(-id, -vc, -cmax, -events)

### Two modellers are going to be tested on this an ML vs a PK approach.
## The goal is to predict the event class.
## First, we split the data into two groups: a training set and a test set. To do this, the createDataPartition function is used:

inTrain <- createDataPartition(
  y = fitdat_gene$Class,
  ## the outcome data
  p = .7,
  ## percentage of data in training set
  list = FALSE
)

## createDataPartition does a stratified random split of the data.
training_gene <- fitdat_gene[ inTrain, ]
testing_gene <- fitdat_gene[-inTrain, ]

nrow(training_gene)
nrow(testing_gene)

### Fit using ML method
## A naive Pk modeller has forgotten the course of dimensionality
## Train logistic model using PK knowledge
pk_training_gene <- training_gene %>%
  mutate(vc = pkpop_vc(bw, sex),
         cmax = (dose / vc) )

logistic <- train(factor(Class) ~ . - dose - sex - bw - vc, 
                  data = pk_training_gene %>%
                    mutate_at("cmax", scale_this), 
                  method = "glm", 
                  trControl = control)
warnings()
cctrl1 <- trainControl(method="cv",
                       number=3,
                       returnResamp="all",
                       classProbs=TRUE,
                       summaryFunction=twoClassSummary)

logisticnet <- train(factor(Class) ~ . - dose - sex - bw - vc, 
                    data = pk_training_gene %>%
                      mutate_at("cmax", scale_this), 
                    method = "glmnet", 
                    metric = "ROC",
                    trControl = cctrl1,
                    tuneGrid = expand.grid(alpha = 1,lambda = seq(0.001,0.1,by = 0.001)))

plot(logisticnet)
# best parameter
coef(logisticnet$finalModel, logisticnet$bestTune$lambda)
glmnet::plot.glmnet(logisticnet$finalModel)

lasso_mod <- logisticnet$finalModel

plot_lasso_norm(lasso_mod)

#### DAGs

library(ggdag)
pk_dag <- dagify( Y ~ Cmax + U, Cmax ~ Dose + VC, VC ~ Sex + Weight, Dose ~ Weight)
pk_dag %>% ggdag() + theme_dag_blank()
