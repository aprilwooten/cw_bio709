#' DESCRIPTION:
#' Script for introductory work


#in-class ----------------------------------------------------------------

library(tidyverse)

# non-parametric test

x <- c(3.2, 5.0, 10.0, 100, 50)
y <- c(1.0, 3.0, 2.0, 2.1, 1.2)

df_xy <- tibble(group = c(rep("x" , length(x)),
                          rep("y", length(y))),
                value = c(x , y))

df_xy %>%
  ggplot(aes(x = group,
             y = value)) +
  geom_boxplot()

# try t-test 

t.test(x, y)


# non-parametric version of t-test
# Wilcox test

wilcox.test(x, y)

# ANOVA is for more than two groups and assuming normal distribution - parametric 

aov(weight ~ group,
    data = PlantGrowth)


# non-parametric version of ANOVA = Kruskal-Wallis Sum Test  kruskal.test()

kruskal.test(weight ~ group,
             data = PlantGrowth)

# Confidence Interval 

# insignificant 
m <- lm(Sepal.Length ~ Sepal.Width,
        data = iris)

confint(m)

# significant 
m1 <- lm(Petal.Length ~ Petal.Width,
        data = iris)

confint(m1)

## correlation ####

x <- rnorm(100, mean = 0, sd = 1)
y <- rnorm(100, mean = 0.8 * x, sd = 1)

plot(x ~ y)

## parametric - pearson 

cor.test(x, y)

## non-parametric - spearman

cor.test(x ,y, method = "spearman")

## covariance 

cov(x, y)

cov(x, y) / (sd(x) * sd(y))







