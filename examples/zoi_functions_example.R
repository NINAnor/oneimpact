library(ggplot2)

# test influence functions

# exponential decay
exp_decay(10, zoi_radius = 30)

f1 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = exp_decay, args = list(zoi_radius = 20)) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f1

# exponential decay - two sides
f1_2 <- ggplot(data.frame(x = c(-30, 30)), aes(x = x)) +
  stat_function(fun = exp_decay,
                args = list(zoi_radius = 20, oneside = FALSE)) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f1_2

# threshold
threshold_decay(5, zoi_radius = 10)
threshold_decay(10, zoi_radius = 10)

f2 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = threshold_decay,
                args = list(zoi_radius = 20), linetype = 2) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f2

# threshold - two sides
f2_2 <- ggplot(data.frame(x = c(-30, 50)), aes(x = x)) +
  stat_function(fun = threshold_decay,
                args = list(zoi_radius = 20, oneside = FALSE), linetype = 2) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f2_2

# linear, tent, or bartlett decay
bartlett_decay(5, zoi_radius = 10)
bartlett_decay(8, zoi_radius = 10)

f3 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = bartlett_decay, args = list(zoi_radius = 20), linetype = 3) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f3

# linear, two sides
f3_3 <- ggplot(data.frame(x = c(-30, 40)), aes(x = x)) +
  stat_function(fun = bartlett_decay,
                args = list(zoi_radius = 20, origin = 10, oneside = FALSE), linetype = 3) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f3_3

# guassian or half normal
gaussian_decay(5, sigma = 6)

f4 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = gaussian_decay,
                args = list(zoi_radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  geom_vline(xintercept = 20, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey") +
  theme_bw()
f4

# half normal - two sides
gaussian_decay(5, sigma = 6)

f4_2 <- ggplot(data.frame(x = c(-30, 30)), aes(x = x)) +
  stat_function(fun = gaussian_decay,
                args = list(zoi_radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  geom_vline(xintercept = c(-20, 20), linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey") +
  theme_bw()
f4_2

# plot several ZoI with the same radius
f1 +
  stat_function(fun = threshold_decay, args = list(zoi_radius = 20), linetype = 2) +
  stat_function(fun = bartlett_decay, args = list(zoi_radius = 20), linetype = 3) +
  stat_function(fun = gaussian_decay, args = list(zoi_radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
