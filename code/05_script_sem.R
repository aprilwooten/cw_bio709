#' DESCRIPTION:
#' Script for SEM

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               GGally,
               vegan,
               lavaan,
               lavaanPlot)

# Specify the URL of the raw CSV file on GitHub
url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_foodweb.csv"

# Read the CSV file into a tibble
(df_fw <- read_csv(url))

# visualize data
df_fw %>%
  select(-plot_id) %>%
  ggpairs(progress = FALSE) +
  theme_bw()

# write a path diagram
m1 <- '
  mass_herbiv ~ mass_plant + cv_h_plant
  mass_pred ~ mass_herbiv'

(fit1 <- sem(model = m1, 
    data = df_fw))

summary(fit1, standardize = TRUE)

lavaanPlot(model = fit1, coefs = TRUE, stand = TRUE)

#test if a specific arrow is important or not using model comparison

m2 <- '
  mass_herbiv ~ mass_plant + cv_h_plant
  mass_pred ~ mass_herbiv + cv_h_plant'

(fit2 <- sem(m2,
             data = df_fw))

# model comparison with anova() i.e., likelihood ratio test 
anova(fit1, fit2)

# interpreting p-value here is opposite
# p<0.05 means model you developed is not great 

#structural equation modeling 
url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_herbivory.csv"

(df_herbv <- read_csv(url))

#visualize
df_herbv %>% 
  ggpairs(
    progress = FALSE,
    columns = c("soil_n",
                "sla",
                "cn_ratio",
                "per_lignin")
  ) +
  theme_bw()

m_sem <- '
# latent variable
  palatability =~ sla + cn_ratio + per_lignin
  
# regression
  palatability ~ soil_n
  herbivory ~ palatability
'

(fit_sem <- sem(m_sem,
                data = df_herbv))

sem(model = m_sem,
    data = df_herbv)

# the limitation and biggest assumption of this method is normal distribution


# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Path Analysis and Covariance Visualization
# ============================================================

library(piecewiseSEM)
data("keeley")

# The "keeley" dataset contains fire-related vegetation data
# collected from shrublands in California.
#
# ------------------------------------------------------------
# Column descriptions:
# elev  : Elevation of the site
# slope : Slope steepness
# aspect: Slope aspect (orientation)
# heat  : Heat load index (a function of slope and aspect)
# firesev: Fire severity
# age   : Time since last fire
# cover : Vegetation cover
# rich  : Plant species richness
# ------------------------------------------------------------
#
# In this exercise, you will explore relationships among variables
# using covariance and path analysis. You will replicate a published
# path model and propose an alternative.

# 1. For the variables depicted in Figure 22.1, draw a figure
#    showing the covariance between variables.

(df_keeley <- keeley %>%
  as_tibble())

df_keeley %>%
  select(-elev) %>%
  ggpairs(progress = FALSE) +
  theme_bw()

# 2. Following Figure 22.1, develop a path model using the
#    same variables and relationships. Examine if this model
#    captures the data structure using a Chi-Square test.

m3 <- '
  rich ~ hetero + cover + abiotic
  abiotic ~ distance
  hetero ~ distance 
  cover ~ firesev 
  firesev ~ age
  age ~ distance'

(fit3 <- sem(model = m3,
             data = df_keeley))

summary(fit3, standardize = TRUE)

lavaanPlot(model = fit3, coefs = TRUE, stand = TRUE)

# 3. Develop an alternative path model that you consider more
#    appropriate based on theory or observed data patterns.

m4 <- '
  rich ~ hetero + cover + abiotic + firesev
  abiotic ~ distance
  hetero ~ distance 
  cover ~ firesev + age
  firesev ~ age + abiotic
  age ~ distance'

(fit4 <- sem(model = m4,
             data = df_keeley))

# 4. Compare the performance of the published model (Figure 22.1)
#    and your alternative model.
#    - Consider fit indices, path coefficients, and interpretability.

anova(fit3, fit4)

