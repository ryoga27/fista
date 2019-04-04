rm(list = ls(all = TRUE))
source("./code/fista.R")

df = generate_sample_data(n = 1000, epsilon_sd = 15)
y = df$y
x = df$x

lambdas = seq(0.01, 0.10, by = 0.01)
out_cv = cv_fista_lasso(
    y = y, x  = x, lambdas = lambdas, iter_num = 100, L = 100*max(eigen(t(x)%*%x)$values)
)
print(out_cv$lambda_opt)
out = fista_lasso(y = y, x = x, lambda = out_cv$lambda_opt, L = 100*max(eigen(t(x)%*%x)$values))
print(out$coef)
