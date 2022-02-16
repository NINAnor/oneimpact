# exponential decay
points <- c(0, 4.5, 7, 7.5)
plot_influence2d(points, zoi = 3, fun = exp_decay, range_plot = c(0, 12), cumulative = FALSE)
plot_influence2d(points, zoi = 3, fun = exp_decay, range_plot = c(0, 12), cumulative = TRUE)

# threshold
plot_influence2d(points, zoi = 3, fun = threshold_decay, range_plot = c(0, 12), cumulative = F)
plot_influence2d(points, zoi = 3, fun = "step_decay", range_plot = c(0, 12), cumulative = T)

# linear
plot_influence2d(points, zoi = 3, fun = bartlett_decay, range_plot = c(0, 12), cumulative = F)
plot_influence2d(points, zoi = 3, fun = "linear_decay", range_plot = c(0, 12), cumulative = T, 
                 return_df = T)

# gaussian
plot_influence2d(points, zoi = 3, fun = gaussian_decay, range_plot = c(0, 12), cumulative = F)
plot_influence2d(points, zoi = 3, fun = "half_norm_decay", range_plot = c(0, 12), cumulative = T, 
                 return_df = T)
