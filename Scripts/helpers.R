scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


plot_lasso_norm <- function(lasso_mod){
  # best coefficient
  beta=coef(lasso_mod)
  tmp <- as.data.frame(as.matrix(beta))
  tmp$coef <- row.names(tmp)
  tmp <- reshape::melt(tmp, id = "coef")
  tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
  tmp$lambda <- lasso_mod$lambda[tmp$variable+1] # extract the lambda values
  tmp$norm <- apply(abs(beta[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm
  
  ggplot(tmp[tmp$coef != "(Intercept)",], aes(norm, value, colour = coef)) + 
    geom_line() + 
    xlab("L1 norm") + 
    theme_bw() + 
    theme(legend.position = "none")
}


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

# Function calculates AUC: 
test_auc <- function(prob) {
  roc(testing$Class, prob) 
}


pr_bernouilli_gene <- function(dat) {
  beta1 <- 2
  gene_marker <- paste("gene_", sample(1:n_gene, 6), sep = "_")
  beta_gene <- c(0.9, -0.2, 0.5, -2.2, 0.7,  1.5)
  
  eta <- scale_this(dat$cmax) * beta1 
  eta <- eta + (  as.matrix( dat %>% select(gene_marker))  %*% beta_gene  ) ## use R matrix multiplication
  pr <- 1 / (1 + exp(eta))
  rbinom(N, 1, pr)
}
