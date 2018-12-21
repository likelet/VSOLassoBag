## Data simulation of biodata which has the feature of multipul relationships among training features.

# Objectives:
# From the literature, it has been reported that Lasso has the vulnerablity of dealing with the data
# whose irrelavent predictors has strong correlations to the true predictors. This vulnerablity
# can be expected especially in biomedical data where biomarkers tends to corralate in a very complex
# way to achieve a specific phenotype on petients. Hence this simulation is designed for testing the
# ability of lassobag in terms of overcoming this vulnerablity of original lasso method.


# function of generating coefficients
coeff <- function (true_predic_num, predic_total_num, distribution=c("norm", 0, 1)) {
  # return a vector of coefficients in simulation
  stopifnot(predic_total_num > true_predic_num)
  if (distribution[1] == "norm") {
    mean <- distribution[2]
    std <- distribution[3]
    true_predic <- rnorm(true_predic_num)
    irre_predic <- rep(0, predic_total_num - true_predic_num)
    return (c(ture_predic, irre_predic))
  }
}


# function of generating true prediction features
true_predictors <- function (true_predic_num, distribution=list("norm", list("random", -10, 10, 0, 5))) {
  # high-order functions to determine which method to use in generating true prediction features
  # parameters in distribution is the lower and uper bound for mean and std, the orther distribution may vary
  if (distribution[1] == "norm") {
    generate <- function(true_predic_num, distribution) {
      if (distribution[[2]][[1]] == "random") {
        mean_lowbound <- -distribution[[2]][[2]]
        mean_upbound <- distribution[[2]][[3]]
        std_lowbound <- distribution[[2]][[4]]
        std_upbound <- distribution[[2]][[5]]
        true_mean <- runif(1, mean_lowbound, mean_upbound) # choose random mean
        true_std <- runif(1, std_lowbound, std_upbound) # choose random std
        return (rnorm(true_predic_num, true_mean, true_std))
      }
    }
  }

  # generate columns and update the matrix
  true_columns <- generate(true_predic_num, distribution)
  for (column in 2:true_predic_num) {
    each_column <- generate(true_predic_num, distribution)
    true_columns <- cbind(true_columns, each_column)
  }
  return (true_columns)
}


# function of generating unrelated feature matrix
unrelated_X <- function (true_predic_num, predic_total_num, distribution=list("norm", list("random", -10, 10, 0, 5))) {
  return (true_predictors(predic_total_num, list("norm", list("random", -10, 10, 0, 5))))
}






# start with simple relationship of A + B = C
