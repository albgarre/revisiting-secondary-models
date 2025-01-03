
library(tidyverse)
library(FME)
library(cowplot)
library(ggsci)

## Define data (doi.org/10.5713/ajas.2013.13579)

d <- list(
  brie = tibble(
    temp = c(10, 15, 25, 30),
    mu = c(0.03, 0.07, 0.45, 0.94),
    lambda = c(20.13, 9.61, 5.40, 5.68)
    ),
  camembert = tibble( 
    temp = c(10, 15, 25, 30),
    mu = c(0.03, 0.09, 0.44, 1.03),
    lambda = c(24.49, 10.35, 5.74, 5.92)
    ),
  mozzarella = tibble(
    temp = c(10, 15, 25, 30),
    mu = c(0.01, 0.06, 0.25, 0.33),
    lambda = c(34.4, 9.22, 4.30, 1.79)
    ),
  cheddar = tibble( 
    temp = c(15, 25, 30),
    mu = c(0.03, 0.18, 0.28),
    lambda = c(35.95, 6.82, 6.08)
    )
  ) %>%
  map(., ~ mutate(., sqmu = sqrt(mu)))

## Visualizations

d %>%
  map(.,
      ~ ggplot(., aes(x = temp, y = sqrt(mu))) +
        geom_point() +
        geom_smooth(method = "lm")
  ) %>%
  plot_grid(plotlist = .)

d %>%
  map(.,
      ~ ggplot(., aes(x = 1/temp^2, y = lambda)) +
        geom_point() +
        geom_smooth(method = "lm")
  ) %>%
  plot_grid(plotlist = .)

## Helpers

pred_sqmu <- function(p, temp) {
  p <- as.list(p)
  p$b*(temp - p$Tmin)
}

residuals_sqmu <- function(p, my_d) {
  p <- as.list(p)
  
  # pred_sqmu(p, my_d$temp) - my_d$sqmu
  
  ## convert to ln CFU/h
  
  mu_log10 <- my_d$sqmu^2
  mu_ln <- log(10)*mu_log10
  
  pred_sqmu(p, my_d$temp) - sqrt(mu_ln)
  

}

# pred_sqmu(c(Tmin = 0, b = .1, C0 = 10), temp = 5:15)
# residuals_sqmu(c(Tmin = 0, b = .1, C0 = 5), d[[1]])

pred_lambda <- function(p, temp) {
  p <- as.list(p)
  
  B <- log( 1 + 1/(10^p$logC0) )/p$b^2
  B/(temp - p$Tmin)^2
}

residuals_lambda <- function(p, my_d) {
  p <- as.list(p)
  
  pred_lambda(p, my_d$temp) - my_d$lambda
}

# pred_lambda(c(Tmin = 0, b = .1, logC0 = 0), temp = 5:15)
# residuals_lambda(c(Tmin = 0, b = .1, logC0 = 0), d[[1]])

cost_all <- function(p, this_data,
                     weight = NULL
                     ) {
  
  r1 <- residuals_sqmu(p, this_data)
  r2 <- residuals_lambda(p, this_data)
  
  if (is.null(weight)) {
    w1 <- 1
    w2 <- 1
  } else if (weight == "sd") {
    w1 <- sd(this_data$sqmu)
    w2 <- sd(this_data$lambda)
  } else if (weight == "mean") {
    w1 <- mean(this_data$sqmu)
    w2 <- mean(this_data$lambda)
  }
  
  c(r1/w1, r2/w2)
  
}

# cost_all(c(Tmin = 0, b = .01, logC0 = 0), d[[2]],
#          weight = "sd")

## Model fitting

my_models <- d %>%
  map(.,
      ~ modFit(cost_all,
               p = c(Tmin = 0, b = .1, logC0 = 0),
               this_data = .,
               weight = "sd"
               )
      )

my_models %>% map(~ summary(.))

## Plots - Figure 2

p <- my_models %>%
  imap_dfr(.,
           ~ tibble(temp = seq(10, 30, length = 100),
                    lambda = pred_lambda(coef(.x), temp),
                    product = .y
                    )
           ) %>%
  ggplot() +
  geom_line(aes(x = temp, y = lambda, colour = product),
            linetype = 2, size = 1) +
  geom_point(aes(x = temp, y = lambda, colour = product),
             data = d %>% imap_dfr(., ~ mutate(.x, product = .y)),
             size = 5, shape = 1
             ) +
  coord_cartesian(ylim = c(0, 40)) +
  labs(x = "Storage temperature (ºC)",
       y = "Lag phase duration (h)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  scale_colour_manual(values = pal_startrek()(4)) +
  facet_wrap("product")

ggsave(p, filename = "Figure_2.png",
       width = 9, height = 6)

## Figure 3

p <- my_models %>%
  imap_dfr(.,
           ~ tibble(temp = seq(10, 30, length = 100),
                    sqmu = pred_sqmu(coef(.x), temp),
                    product = .y
           )
  ) %>%
  ggplot() +
  geom_line(aes(x = temp, y = sqmu^2, colour = product),
            linetype = 2, size = 1) +
  geom_point(aes(x = temp, y = mu*log(10),  # convert to ln
                 colour = product),
             data = d %>% imap_dfr(., ~ mutate(.x, product = .y)),
             size = 5, shape = 1
  ) +
  labs(x = "Storage temperature (ºC)",
       y = expression(Maximum~specific~growth~rate~(ln~CFU/h))
       ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  scale_colour_manual(values = pal_startrek()(4)) +
  facet_wrap("product")

ggsave(p, filename = "Figure_3.png",
       width = 9, height = 6)

## Table 2

my_models %>%
  map(., ~ summary(.)$par)

d %>%
  map2(., my_models,
       ~ cost_all(.y$par, .x, weight = "sd")
       ) %>%
  map(., ~tibble(r = ., r2 = r^2)) %>%
  map(., ~ summarize(., SSE = mean(r2))) %>%
  map(., ~ mutate(., RMSE = sqrt(SSE)))


## Fitting just the growth model

d %>%
  map(.,
      ~ modFit(residuals_sqmu,
               p = c(Tmin = 0, b = .1),
               my_d = .
      )
  ) %>% map(., ~summary(.)$par)



























