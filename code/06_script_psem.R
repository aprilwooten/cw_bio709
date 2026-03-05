#' DESCRIPTION:
#' Script for piecewise SEM

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               GGally,
               piecewiseSEM,
               glmmTMB)

data("keeley")

(df_keeley <- keeley %>% 
    as_tibble())

# abiotic, age, hetero, firesev, cover, rich are response variables
# we assume normal distribution
# richness is a discreet variable, and may benefit from non-normal distribution

#psem, but all normally distributed 
# Define individual models using lm function
m1 <- lm(abiotic ~ distance, data = df_keeley)
m2 <- lm(hetero ~ distance, data = df_keeley)
m3 <- lm(firesev ~ age, data = df_keeley)
m4 <- lm(cover ~ firesev, data = df_keeley)
m5 <- lm(rich ~ cover + abiotic + hetero, data = df_keeley)

# Combine into piecewise SEM
sem_model <- psem(m1, m2, m3, m4, m5)

# Evaluate
summary(sem_model, .progressBar = FALSE)

# psem with negative binomial function
# Define individual models
m1 <- lm(abiotic ~ distance, data = df_keeley)      # unchanged
m2 <- lm(hetero ~ distance, data = df_keeley)       # unchanged
m3 <- lm(firesev ~ age, data = df_keeley)           # unchanged

# m4 now includes a direct effect of hetero on cover (added path)
m4 <- lm(cover ~ firesev + hetero, data = df_keeley)  

# m5 now models richness as negative binomial (MASS::glm.nb) 
# and includes direct effect of distance on richness (added path)
m5 <- MASS::glm.nb(rich ~ cover + abiotic + hetero + distance, 
                   data = df_keeley)  

# Combine into piecewise SEM
sem_model <- psem(m1, m2, m3, m4, m5)

# Evaluate model
summary(sem_model, .progressBar = FALSE)

## fishers C accounts for ommitted variables, follows Chi squared distribution

plot(sem_model)

data("shipley")

df_shipley <- shipley %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  drop_na(growth)

df_shipley %>% 
  group_by(site) %>% 
  summarize(n_tree = n_distinct(tree))

# visualization
# Pairwise plots of key variables
df_shipley %>% 
  ggpairs(
    columns = c("dd", 
                "date",
                "growth",
                "live"),  # select variables to plot
    progress = FALSE      # disable progress bar
  ) +
  theme_bw()              # clean theme for the plot

# Model 1: date depends on dd, with random intercepts for site and tree
m1 <- glmmTMB(date ~ dd + (1 | site) + (1 | tree), 
              data = df_shipley,
              family = "gaussian")

# Model 2: growth depends on date, same random effects
m2 <- glmmTMB(growth ~ date + (1 | site) + (1 | tree), 
              data = df_shipley,
              family = "gaussian")

# Model 3: live (binary) depends on growth, logistic mixed model
m3 <- glmmTMB(live ~ growth + (1 | site) + (1 | tree), 
              data = df_shipley, 
              family = "binomial")

# Combine models into a piecewise SEM
sem_glmm <- psem(m1, m2, m3)

# Summarize SEM (paths, significance, and Shipley's test)
summary(sem_glmm, .progressBar = FALSE)


# lab ---------------------------------------------------------------------

library(piecewiseSEM)
data("meadows")

# =========================================
# EXERCISE: Piecewise SEM with Meadows Data
# =========================================
#
# ------------------------------------------------------------
# Dataset: meadows (from piecewiseSEM package)
# Variables:
#   grazed - 0 = ungrazed, 1 = grazed
#   mass   - plant biomass (g/m²)
#   elev   - plot elevation above sea level
#   rich   - plant species richness per m²
# ------------------------------------------------------------
#
# 1. Explore the dataset (structure, summary, plots).

(df_meadows <- meadows %>%
   as_tibble())

df_meadows %>% 
  ggpairs(
    columns = c("grazed", 
                "mass",
                "elev",
                "rich"),  
    progress = FALSE      
  ) +
  theme_bw()  

# 2. Develop a conceptual model: decide which variables influence others.
#    - Consider direct and indirect effects.
#    - Think about grazing as a disturbance factor.



# 3. Fit component models (e.g., lm) for each hypothesized relationship.

m1 <-lm(grazed ~ elev,
           data = df_meadows)

m2 <-lm(mass ~ grazed + rich,
           data = df_meadows)

m3 <- lm(rich ~ grazed + elev,
         data = df_meadows )


# 4. Combine models into a piecewise SEM using psem().

sem_model <- psem(m1, m2, m3)

# 5. Evaluate the SEM: path coefficients, significance, variance explained.

summary(sem_model)

# 6. Optional: try alternative models if your model deviates from the expectation.
#
# Deliverables:
# - Code for component models and combined SEM
# - Conceptual SEM diagram
# - Short reasoning about your SEM results
