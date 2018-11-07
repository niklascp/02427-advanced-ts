# Part 1: Simulation and discretization of diffusion processes
delta <- 2^-9
theta_1 <- .7
theta_2 <- .8
theta_3 <- 3.0
theta_4 <- -0.34
N <- 100/delta
Y_1_init <- -1.9
Y_2_init <- 1.2

run_sim <- function(sigma) {
    Y <- matrix(nrow = N, ncol = 2)
    Y[1, 1] <- Y_1_init
    Y[1, 2] <- Y_2_init
    W <- rnorm(N, 0, delta)

    for (n in 1:(N-1)) {
        Y[n + 1, 1] <- Y[n, 1] + theta_3 * (Y[n, 1] + Y[n, 2] - 1/3*(Y[n, 1]^3) + theta_4) * delta + sigma * W[n + 1]
        Y[n + 1, 2] <- Y[n, 2] - (1/theta_3) * (Y[n, 1] + theta_2 * Y[n, 2] - theta_1) * delta
    }
    df <- data.frame(Y)
    colnames(df) <- c('Y1', 'Y2')
    df$n <- seq(N)
    df$t <- df$n * delta
    df
}

# Question 1a
df_1 <- run_sim(0)

library(ggplot2)

generate_plots <- function(df, prefix) {
    gg <- ggplot(df) +
        geom_line(aes(x = t, y = Y1))
    ggsave(paste0(prefix, '-Y1.pdf'), plot = gg, width = 5, height = 3)

    gg <- ggplot(df) +
        geom_line(aes(x = t, y = Y2))
    ggsave(paste0(prefix, '-Y2.pdf'), plot = gg, width = 5, height = 3)

    gg <- ggplot(df) +
        geom_point(aes(x = Y1, y = Y2))
    ggsave(paste0(prefix, '-Y1Y2.pdf'), plot = gg, width = 5, height = 3)
}

generate_plots(df_1, 'latex/part1a-sigma0')
generate_plots(run_sim(.10), 'latex/part1a-sigma1')
generate_plots(run_sim(.20), 'latex/part1a-sigma2')
generate_plots(run_sim(.30), 'latex/part1a-sigma3')
generate_plots(run_sim(.40), 'latex/part1a-sigma4')

# Question 1b
library('reshape2')
library('data.table')
plot_grid <- function(df, prefix) {
    setDT(df)
    df$Y1_ <- cut(df$Y1, 100)
    df$Y2_ <- cut(df$Y2, 100)

    dt_agg <- df[ , .( Count = .N ) , by = list(Y1_, Y2_) ]

    p <- ggplot(data = dt_agg, aes(x = Y1_, y = Y2_)) +
        geom_tile(aes(fill = Count)) +
        scale_x_discrete(breaks = NULL) + 
        scale_y_discrete(breaks = NULL)
    ggsave(paste0(prefix, '-grid.pdf'), plot = p, width = 7, height = 7)
}
plot_grid(run_sim(.10), 'latex/part1b-sigma1')
plot_grid(run_sim(.20), 'latex/part1b-sigma2')
plot_grid(run_sim(.30), 'latex/part1b-sigma3')
plot_grid(run_sim(.40), 'latex/part1b-sigma4')
