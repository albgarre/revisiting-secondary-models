
library(tidyverse)
library(readxl)
library(FME)
library(cowplot)
library(ggsci)

## Load data

d <- read_excel("./data/data_gerardo.xlsx") %>%
  filter(Parameter %in% c("SL", "logD")) %>%
  rename(temp = `Temperature (ºC)`) %>%
  select(-`Std. Error`) %>%
  pivot_wider(names_from = "Parameter", values_from = "Estimate") %>%
  mutate(D = 10^logD,
         k = log(10)/D,
         logk = log(k)) %>%
  select(-logD, -D) %>%
  mutate(SL = ifelse(is.na(SL), 0, SL))

ggplot(d) +
  geom_point(aes(x = temp, y = logk))

ggplot(d) +
  geom_point(aes(x = temp, y = SL)) +
  scale_y_log10()

## Fitting functions

pred_SL <- function(p, temp) {
  p <- as.list(p)
  log(10^p$logC0 + 1)/exp(p$a + p$b*temp)
}

residuals_SL <- function(p, my_d, logtransf = FALSE) {
  p <- as.list(p)
  
  pred <- tibble(temp = my_d$temp,
                 SL = pred_SL(p, temp = temp)
  )
  
  pred$SL - my_d$SL
  
}

pred_k <- function(p, temp) {
  p <- as.list(p)
  p$a + p$b*temp
}

residuals_logk <- function(p, my_d) {
  p <- as.list(p)
  
  pred <- tibble(temp = my_d$temp,
                 logk = pred_k(p, temp)
  )
  
  pred$logk - my_d$logk
}

cost_all <- function(p, this_data,
                     weight = NULL
                     ) {
  
  r1 <- residuals_SL(p, this_data)
  r2 <- residuals_logk(p, this_data)
  
  if (is.null(weight)) {
    w1 <- 1
    w2 <- 1
  } else if (weight == "sd") {
    w1 <- sd(this_data$SL)
    w2 <- sd(this_data$logk)
  }
  
  c(r1/w1, r2/w2)
  
}

# pred_SL(c(logC0 = 2, a = -8, b = .15), seq(50, 60))
# pred_k(c(logC0 = 2, a = -8, b = .15), seq(50, 60))
# 
# cost_all(c(logC0 = 2, a = -8, b = .15), d)
# cost_all(c(logC0 = 2, a = -8, b = .15), d, weight = "sd")

## Model fitting

coco_model <- modFit(cost_all,
                     p = c(logC0 = 2, a = -8, b = .15),
                     this_data = filter(d, Medium == "Coconut"),
                     weight = "sd"
                     )

TSB_model <- modFit(cost_all,
                     p = c(logC0 = 2, a = -8, b = .15),
                     this_data = filter(d, Medium == "TSB"),
                     weight = "sd"
                    )

summary(coco_model)
summary(TSB_model)

## Plots - Figure 1

p1 <- tibble(
  temp = seq(50, 60, length = 100),
  Coconut = pred_SL(coef(coco_model), temp),
  TSB = pred_SL(coef(TSB_model), temp)
  # logk = pred_k(coef(coco_model), temp),
) %>%
  pivot_longer(-temp) %>%
  ggplot() +
  geom_line(aes(x = temp, y = value, colour = name),
            linetype = 2, size = 1) +
  geom_point(aes(x = temp, y = SL, colour = Medium), data = d,
             size = 5, shape = 1) +
  # geom_errorbar(aes(x = `Temperature (ºC)`, 
  #                   ymin = Estimate - `Std. Error`,
  #                   ymax = Estimate + `Std. Error`, 
  #                   colour = Medium), 
  #               width = .1,
  #               data = read_excel("./data/data_gerardo.xlsx") %>% filter(Parameter == "SL")
  #               ) +
  labs(x = "Treatment temperature (ºC)",
       y = "Shoulder length (min)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  scale_colour_manual(values = pal_startrek()(2))


p2 <- tibble(
  temp = seq(50, 60, length = 100),
  Coconut = pred_k(coef(coco_model), temp),
  TSB = pred_k(coef(TSB_model), temp)
  # logk = pred_k(coef(coco_model), temp),
) %>%
  pivot_longer(-temp) %>%
  ggplot() +
  geom_line(aes(x = temp, y = value, colour = name),
            linetype = 2, size = 1) +
  geom_point(aes(x = temp, y = logk, colour = Medium), data = d,
             size = 5, shape = 1) +
  labs(x = "Treatment temperature (ºC)",
       y = expression(ln~k~(ln~min^-1))
       ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        legend.title = element_blank()
        ) +
  scale_colour_manual(values = pal_startrek()(2))

p <- plot_grid(p1, p2)

ggsave(p, filename = "Figure_1.png",
       width = 12, height = 6)






























