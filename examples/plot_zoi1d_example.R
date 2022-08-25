# one point, exponential decay
plot_zoi1d(0, radius = 3, zoi_limit = 0.05,
           fun = exp_decay, range_plot = c(0, 5), cumulative = FALSE)

# exponential decay
points <- c(0, 4.5, 7, 7.5)
plot_zoi1d(points, radius = 3, fun = exp_decay, range_plot = c(0, 12), cumulative = FALSE)
plot_zoi1d(points, radius = 3, fun = exp_decay, range_plot = c(0, 12), cumulative = TRUE)

# threshold
plot_zoi1d(points, radius = 3, fun = threshold_decay, range_plot = c(0, 12), cumulative = F)
plot_zoi1d(points, radius = 3, fun = "step_decay", range_plot = c(0, 12), cumulative = T)

# linear
plot_zoi1d(points, radius = 3, fun = bartlett_decay, range_plot = c(0, 12), cumulative = F)
plot_zoi1d(points, radius = 3, fun = "linear_decay", range_plot = c(0, 12), cumulative = T,
                 return_df = T)

# gaussian
plot_zoi1d(points, radius = 3, fun = gaussian_decay, range_plot = c(0, 12), cumulative = F)
plot_zoi1d(points, radius = 3, fun = "half_norm_decay", range_plot = c(0, 12), cumulative = T,
                 return_df = T)
