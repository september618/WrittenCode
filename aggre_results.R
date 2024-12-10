## Merge results ------

library(gtools)
library(purrr)

NN_values = c(500, 800, 1000, 2000, 4000, 8000)

sum_all = lapply(NN_values,function(NN){
  cat("NN = ",NN,"\n")

result_files = mixedsort(list.files(path = paste0("~/Desktop/resultsNN_",NN)))
results_list = lapply(result_files, function(file) {
  load(paste0("~/Desktop/resultsNN_",NN,"/", file))
  list(
    truth = truth,
    conv_F = fit_F$convergence,
    conv_M = fit_M$convergence,
    conv_C = fit_C$convergence,
    conv_T = fit_T$convergence,
    est_F = est_F,
    est_M = est_M,
    est_C = est_C,
    est_T = est_T,
    se_F = se_F,
    se_M = se_M,
    se_C = se_C,
    se_T = se_T,
    t_F = t_F,
    t_M = t_M,
    t_C = t_C,
    t_T = t_T
  )
})
# Combine results from all files
truth = do.call(rbind,lapply(results_list, `[[`, "truth")) %>%t()
conv_F = do.call(cbind, lapply(results_list, `[[`, "conv_F"))
conv_M = do.call(cbind, lapply(results_list, `[[`, "conv_M"))
conv_C = do.call(cbind, lapply(results_list, `[[`, "conv_C"))
conv_T = do.call(cbind, lapply(results_list, `[[`, "conv_T"))
est_F = do.call(cbind, lapply(results_list, `[[`, "est_F"))
est_M = do.call(cbind, lapply(results_list, `[[`, "est_M"))
est_C = do.call(cbind, lapply(results_list, `[[`, "est_C"))
est_T = do.call(cbind, lapply(results_list, `[[`, "est_T"))
se_F = do.call(cbind, lapply(results_list, `[[`, "se_F"))%>%apply(2,as.numeric)
se_M = do.call(cbind, lapply(results_list, `[[`, "se_M"))%>%apply(2,as.numeric)
se_C = do.call(cbind, lapply(results_list, `[[`, "se_C"))%>%apply(2,as.numeric)
se_T = do.call(cbind, lapply(results_list, `[[`, "se_T"))%>%apply(2,as.numeric)
t_F = do.call(cbind, lapply(results_list, `[[`, "t_F"))
t_M = do.call(cbind, lapply(results_list, `[[`, "t_M"))
t_C = do.call(cbind, lapply(results_list, `[[`, "t_C"))
t_T = do.call(cbind, lapply(results_list, `[[`, "t_T"))

# check convergence
conv_F_rate = rowMeans(conv_F == "Successful")
conv_M_rate = rowMeans(conv_M == "Successful")
conv_C_rate = rowMeans(conv_C == "Successful")
conv_T_rate = rowMeans(conv_T == "Successful")

idx_non_conv_F = which(conv_F!="Successful")
idx_non_conv_M = which(conv_M!="Successful")
idx_non_conv_C = which(conv_C!="Successful")
idx_non_conv_T = which(conv_T!="Successful")


# check bias se and mse
bias_F = rowMeans(est_F - truth)
bias_M = rowMeans(est_M - truth)
bias_C = rowMeans(est_C - truth)
bias_T = rowMeans(est_T - truth)
se_F_avg = rowMeans(se_F,na.rm = T)
se_M_avg = rowMeans(se_M,na.rm = T)
se_C_avg = rowMeans(se_C,na.rm = T)
se_T_avg = rowMeans(se_T,na.rm = T)
mse_F = bias_F^2 + se_F_avg^2
mse_M = bias_M^2 + se_M_avg^2
mse_C = bias_C^2 + se_C_avg^2
mse_T = bias_T^2 + se_T_avg^2

# check time
t_F_avg = rowMeans(t_F)
t_M_avg = rowMeans(t_M)
t_C_avg = rowMeans(t_C)
t_T_avg = rowMeans(t_T)


# check coverage
ci_lower_F = est_F - 1.96*se_F
ci_upper_F = est_F + 1.96*se_F
coverage_F = rowMeans(((ci_lower_F) <= truth) & (truth <= (ci_upper_F)),na.rm = T)

ci_lower_M = est_M - 1.96*se_M
ci_upper_M = est_M + 1.96*se_M
coverage_M = rowMeans(((ci_lower_M) <= truth) & (truth <= (ci_upper_M)),na.rm = T)

ci_lower_C = est_C - 1.96*se_C
ci_upper_C = est_C + 1.96*se_C
coverage_C = rowMeans(((ci_lower_C) <= truth) & (truth <= (ci_upper_C)),na.rm = T)

ci_lower_T = est_T - 1.96*se_T
ci_upper_T = est_T + 1.96*se_T
coverage_T = rowMeans(((ci_lower_T) <= truth) & (truth <= (ci_upper_T)),na.rm = T)

rst = list(
  NN = NN,
  truth = truth[,1],
  bias_F = bias_F,
  bias_M = bias_M,
  bias_C = bias_C,
  bias_T = bias_T,
  se_F_avg = se_F_avg,
  se_M_avg = se_M_avg,
  se_C_avg = se_C_avg,
  se_T_avg = se_T_avg,
  mse_F = mse_F,
  mse_M = mse_M,
  mse_C = mse_C,
  mse_T = mse_T,
  t_F_avg = t_F_avg,
  t_M_avg = t_M_avg,
  t_C_avg = t_C_avg,
  t_T_avg = t_T_avg,
  coverage_F = coverage_F,
  coverage_M = coverage_M,
  coverage_C = coverage_C,
  coverage_T = coverage_T
)
return(rst)
})




# Combine results ------
library(dplyr)
library(tidyr)

# Helper function to flatten results
flatten_result = function(rst) {
  flat_list = lapply(names(rst), function(name){
    value = rst[[name]]
    data.frame(name = value,stringsAsFactors = FALSE) %>%rename(!!name := name)
  })
  df = do.call(cbind, flat_list)
  df %>% mutate(param = row.names(df))
}

results_df = do.call(rbind, lapply(sum_all, flatten_result))

## Plot results ------
library(ggplot2)
library(ggthemes)
results_df %>% pivot_longer(cols = starts_with(c('bias', 'se', 'mse', 'coverage')),
             names_to = c('metric', 'method'),
             names_sep = '[_]',
             values_to = "value") %>%
  filter(param %in% c('mean', 'sill')) %>%
  mutate(method = factor(method, levels = c('M', 'C', 'F', 'T')),
         metric = factor(metric, levels = c('bias', 'se', 'mse', 'coverage'))) %>%
  ggplot(aes(x = NN, y = value, color = method,
             shape = method)) +
  geom_point(size = 3) +geom_line(linewidth = 0.8) +
  facet_wrap(param  ~ metric,
             scale = 'free',
             ncol = 4,
             labeller = as_labeller(c('mean' = 'mu',
                                      'sill' = 'sigma2',
                                      'bias' = 'bias',
                                      'se' = 'standard error',
                                      'mse' = 'mean squared error',
                                      'coverage' = 'coverage'))) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = 'top') + 
  scale_color_discrete('', labels = c('F' = 'Full',
                                      'M' = 'Marginal',
                                      'C' = 'Conditional',
                                      'T' = 'Tapered')) + 
  scale_shape_discrete('', labels = c('F' = 'Full',
                                      'M' = 'Marginal',
                                      'C' = 'Conditional',
                                      'T' = 'Tapered')) +
  geom_hline(data = data.frame(metric = factor('coverage',
                                               levels = c('bias', 'se', 'mse', 'coverage')),
                               value = 0.95),
             aes(yintercept = value), 
             color = 'gray', linetype = 2) + 
  xlab('Size of spatial locations') + ylab(' ') 
  # theme_get()

results_df %>% pivot_longer(cols = starts_with("t") & !matches("^truth$"),
                            names_to = c('metric', 'method'),
                            names_sep = '[_]', 
                            values_to = "value") %>%
  ggplot(aes(x = NN, y = value, color = method,
             shape = method)) +
  geom_point(size = 3) +geom_line(linewidth = 0.8) +
  facet_wrap(~ metric,
             scale = 'free',
             ncol = 4) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = 'top')

results_df%>%select(NN, t_F_avg, t_M_avg, t_C_avg, t_T_avg) %>% unique()








