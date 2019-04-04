lasso_penalty_function = function(x, lambda){
    out = 0
    d = length(x)
    for(j in 1:d){
        out = out + lambda*abs(x[j])
    }
    return(out)
}

lasso_penalty_solution = function(x, lambda){
    d = length(x)
    out = rep(NA, d)
    for(j in 1:d){
        if(abs(x[j]) <= lambda){
            out[j] = 0
        }
        if(abs(x[j]) >  lambda){
            out[j] = sign(x[j])*(abs(x[j]) - lambda)
        }
    }
    return(out)
}


rescale_beta_mat = function(beta_mat, x_sd, iter_num){
    d = nrow(beta_mat)
    beta_mat = beta_mat[, 1:(iter_num + 1)]
    for(i in 1:(iter_num + 1)){
        for(j in 1:d){
            beta_mat[j, i] = beta_mat[j, i]/x_sd[j]
        }
    }
    return(beta_mat)
}

fista_lasso = function(
    y, x, iter_max = 10000, eps = 0.001, lambda = 1, L = NA,
    process_save = TRUE
){

    setup = list(
        y = y, x = x,
        iter_max = iter_max,
        lambda = lambda,
        eps = eps,
        L = L
    )

    d = ncol(x)
    n = nrow(x)

    y_mean = mean(y)
    Y = y - y_mean

    X = scale(x)
    x_mean = attr(X, "scaled:center")
    x_sd   = attr(X, "scaled:scale")

    beta_mat = matrix(NA, nrow = d, ncol = iter_max + 1)
    beta_mat[, 1] = t(X)%*%y
    w_mat = matrix(NA, nrow = d, ncol = iter_max + 1)
    w_mat[, 1] = beta_mat[, 1]
    t_vec = rep(NA, iter_max + 1)
    t_vec[1] = 0
    obj = rep(NA, iter_max + 1)
    obj[1] =  (1/n)*sum((Y - X%*%beta_mat[, 1])^2) + lasso_penalty_function(beta_mat[, 1], lambda)
    # Lipschitz constant
    if(is.na(L)){
        L = 10*max(eigen(t(X)%*%X)$values)
    }

    for(i in 1:iter_max){
        v = w_mat[, i] + (1/(L*lambda))*t(X)%*%(Y - X%*%w_mat[, i])
        beta_mat[, i + 1] = lasso_penalty_solution(v, lambda)
        t_vec[i + 1] = (1/2)*(1 + sqrt(1 + 4*t_vec[i]^2))
        w_mat[, i + 1] = beta_mat[, i] + (t_vec[i] - 1)/t_vec[i + 1]*(beta_mat[, i + 1] - beta_mat[, i])
        obj[i + 1] = (1/n)*sum((Y - X%*%beta_mat[, i + 1])^2) + lasso_penalty_function(beta_mat[, i + 1], lambda)
        if(abs(obj[i + 1] - obj[i]) < eps){
            iter_num = i
            break
        }
        if(i == iter_max){
            iter_num = iter_max
            cat("not convergence \n")
            break
        }
    }

    beta_mat = rescale_beta_mat(beta_mat, x_sd, iter_num)
    nu = y_mean - x_mean%*%beta_mat[, iter_num + 1]

    if(process_save == TRUE){
        process = list(
            iter_num = iter_num,
            beta_mat = beta_mat,
            nu = nu,
            obj = obj[1:(iter_num + 1)]
        )
        result = list(
            nu = nu,
            beta = beta_mat[, iter_num + 1],
            coefficients = c(nu, beta_mat[, iter_num + 1]),
            setup = setup,
            process = process
        )
    }
    if(process_save == FALSE){
        result = list(
            nu = nu,
            beta = beta_mat[, iter_num + 1],
            coefficients = c(nu, beta_mat[, iter_num + 1]),
            setup = setup,
            process = process
        )
    }

    return(result)
}

predict = function(result, x_new){
    nu = result$nu
    beta = result$beta
    n = nrow(x_new)
    y_pred = rep(NA, length = n)
    for(i in 1:n){
        y_pred[i] = nu + beta%*%x_new[i, ]
    }
    return(y_pred)
}

cv_fista_lasso = function(
    y, x, lambdas, L = NA, K = NA, K_rate = 1/3, iter_num = 1000, display_process = TRUE
){
    n = length(y)
    lambdas_num = length(lambdas)

    if(is.na(K)){
        K = round(n*K_rate, digits = 0)
    }

    se_mat = matrix(NA, nrow = iter_num, ncol = lambdas_num)
    mse_vec = rep(NA, length = lambdas_num)

    for(i in 1:iter_num){
        for(l in 1:lambdas_num){
            id_sample = sample(x = 1:n, size = n, replace = FALSE)
            id_train = id_sample[1:(n - K)]
            id_test = id_sample[(n - K + 1):n]
            y_train = y[id_train]
            x_train = x[id_train, ]
            y_test = y[id_test]
            x_test = x[id_test, ]
            result_list = fista_lasso(y = y_train, x = x_train, lambda = lambdas[l], L = L)
            y_pred = predict(result_list, x_test)
            se_mat[i, l] = mean((y_test - y_pred)^2)
        }
        if(display_process == TRUE){
            if(i %% 10 == 0){
                cat("number of cross varidation is:", i, "\n")
            }
        }
    }
    for(l in 1:lambdas_num){
        mse_vec[l] = mean(se_mat[, l])
    }
    id_lambda_opt = min((1:L)[mse_vec == min(mse_vec)])
    lambda_opt = lambdas[id_lambda_opt]

    result = list(
        lambda_opt = lambda_opt,
        mse = mse_vec
    )

    return(result)
}

generate_sample_data = function(
    n = 100, p = 10,
    const = 1, beta = c(1, -2, 3, -4, 5, rep(0, length = p - 5)),
    epsilon_sd = 3
){
    epsilon = matrix(rnorm(n, mean = 0, sd = epsilon_sd))
    x = matrix(runif(n*p, min = -5, max = 5), nrow = n, ncol = p)
    y = const + x%*%beta + epsilon
    df = list(
        y = y,
        x = x
    )
    return(df)
}
