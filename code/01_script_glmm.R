#' DESCRIPTION:
#' Script for GLMMs

# in-class ----------------------------------------------------------------

# sampling design defines which analysis is most appropriate 
# after reading the paper, the group structure was the nest
# because there were multiple measurements in each nest but they are not independent 
# of one another, we have to account for the known independents emerging from each data structure 

# account for group structure before examining the effect of our interest (nutrient in this case)
# multiple individuals are measured from the same nest, they were not independent of each other
# so you have to account for them being together 
# could be some shared environmental context within each nest impacting individuals 

pacman::p_load(tidyverse,
               lme4,
               glmmTMB)

## call sample data 

data(Owls)

df_owl_raw <- as_tibble(Owls)
  
(df_owl <- df_owl_raw %>%                 # Start with raw owl dataset and assign cleaned version
     janitor::clean_names() %>%              # Standardize column names (lowercase, underscores, etc.)
     mutate(across(.cols = where(is.factor), # Select all factor-type columns
                   .fns = str_to_lower))     # Convert factor levels to lowercase
  )  

# food treatment is predictor effect, and the response variable is the sibling negotiation
# whether or not they got food changes the sibling negotiation 

# want to see the general pattern in the data, number of negotiation observed
# divided by the number of chicks 

## visualization 
# x is treatment we are interested in, and y is the negotiation per chick 

df_owl %>%                          # Use cleaned owl dataset
  ggplot(aes(x = food_treatment,    # Map food treatment to x-axis
             y = neg_per_chick)) +  # Map number of negatives per chick to y-axis
  geom_boxplot(outliers = FALSE) +  # Draw boxplots without plotting outliers
  geom_jitter(alpha = 0.25) +       # Add jittered raw data points with transparency
  theme_bw()                        # Apply a clean black-and-white theme

# we can see negotiation decreased with satiated 

## ignoring nest effects (all nests respond identically ... but that's not reality)
MASS::glm.nb(sibling_negotiation ~ food_treatment + offset(log(brood_size)),
             data = df_owl)

# function to perform GLM with negative binomial distribution
# chose neg binomial distribution (glm.nb) because the response variable is count
# when count data seems really variable among data points (if variance exceeds average) then 
# poisson distribution is not appropriate, so use neg binomial dist 

# number of negotiation per chick is decreased in a satiated scenario as indicated
# by the -0.5200 under food_treatmentsatiated when ran 

# but this does not account for the differences between nests, but we have many observations
# from the same nest, likely siblings so they  have some sort of relationship with one another
# and may share environmental factors 

## adding nest effects as fixed effects (simplest way to account for nest difference)
MASS::glm.nb(sibling_negotiation ~ food_treatment + nest + offset(log(brood_size)),
             data = df_owl)
# hard to interpret the results
# yes, you want to account for group structure but you need to incorporate random effect
# you are not interested in nest as a fixed effect with un-interpretable numbers,
# but you will still be accounting for it as a random effect 
# makes more sense to include as a fixed effect if the number of nests was small, but that 
# is not the case here 

## visualization for group-level figures 
v_g9 <- unique(df_owl$nest)[1:9]    # Extract first 9 nest IDs (even though there are many more)

df_owl %>%                          # Start with owl dataset
  filter(nest %in% v_g9) %>%        # Keep only observations from selected nests
  ggplot(aes(x = food_treatment,    
             y = neg_per_chick)) +  
  geom_jitter(alpha = 0.25,         # Draw jitter plots for each treatment
              width = 0.1) +       
  facet_wrap(facets =~ nest,        # Create separate panels for each nest
             ncol = 3,              # Arrange panels into 3 columns
             nrow = 3) +            # Arrange panels into 3 rows
  theme_bw()    # Apply clean black-and-white theme

# clear that effect of satiation is different among nests 
# fixed effect compares every variable against another one (each nest to another) (26 parameters)
# random effects determines how variable the nests are from one another (1 parameter, as 
# standard deviation) increases statistical power to determine the effect of the factor you 
# are interested in 

## use glmmTMB to include random effects 
## regular glm = y ~ x 
## with random effect = y ~ x + (1 | col name that indicates group structure)


m_ri <- glmmTMB(
  sibling_negotiation ~       # Response variable
    food_treatment +          # Fixed effect of food treatment
    (1 | nest) +              # Random intercept for each nest
    offset(log(brood_size)),  # Offset to model rate per chick
  data = df_owl,              
  family = nbinom2()          # Negative binomial distribution
)

summary(m_ri)                 # Display model results and parameter estimates

## get random intercept values
# control group is represented by the intercept value in this case
# we are expected different level of negotiation by nest, and you change it by what 
# value was observed in the data set and compare it to the control 
# each nest has its own unique value, instead of using an average
# accounting for this unique value one by one before you evaluate treatment, so you can
# better account for the effect of the treatment on the nest instead of using an average 
# of all the nest's values 

head(coef(m_ri)$cond$nest)

# we assumed that the intercept was changing, or expected number of negotiation without 
# additional food, so that is why the intercept changed but the food_treatmentsatiated
# value did not change 
# the intercept reflects the negotiation regardless of food, some nests are just more 
# active than others 

# ----------not expected to perform this on my own, super in-depth-------------
# 1. Global (population-level) intercept on the response scale
# ------------------------------------------------------------
# The model uses a log link (negative binomial),
# so the intercept is on the log scale.
# We exponentiate it to return to the original response scale.
g0 <- exp(fixef(m_ri)$cond[1])


# ------------------------------------------------------------
# 2. Select a subset of nests to visualize
# ------------------------------------------------------------
# Plotting all 27 nests would be visually overwhelming,
# so we randomly select 9 nests for this example.
set.seed(123)  # ensures reproducibility
v_g9_ran <- sample(unique(df_owl$nest),
                   size = 9)


# ------------------------------------------------------------
# 3. Extract nest-specific coefficients (random intercept model)
# ------------------------------------------------------------
# coef(m_ri) returns the sum of fixed + random effects for each nest.
# These values are still on the log scale.
df_g9 <- coef(m_ri)$cond$nest %>% 
  as_tibble(rownames = "nest") %>%      # convert to tibble and keep nest ID
  filter(nest %in% v_g9_ran) %>%        # keep only the selected nests
  rename(
    log_g = `(Intercept)`,              # nest-specific intercept (log scale)
    b = food_treatmentsatiated          # fixed slope for food treatment
  ) %>% 
  mutate(
    g = exp(log_g),                     # intercept on response scale
    s = exp(log_g + b)                  # predicted value under satiated treatment
  )


# ------------------------------------------------------------
# 4. Create the figure
# ------------------------------------------------------------
df_owl %>% 
  filter(nest %in% v_g9_ran) %>%        # plot only the selected nests
  ggplot(aes(x = food_treatment,
             y = neg_per_chick)) +
  # Raw data points (jittered to reduce overlap)
  geom_jitter(width = 0.1,
              alpha = 0.5) +
  # Dashed horizontal lines:
  # nest-specific intercepts (baseline differences among nests)
  geom_hline(data = df_g9,
             aes(yintercept = g),
             alpha = 0.5,
             linetype = "dashed") +
  # Solid line segments:
  # predicted change from unfed to satiated treatment
  # using a common (fixed) slope across nests
  geom_segment(data = df_g9,
               aes(y = g,
                   yend = s,
                   x = 1,
                   xend = 2),
               linewidth = 0.5,
               linetype = "solid") +
  # Solid blue horizontal line:
  # global (population-level) intercept
  geom_hline(yintercept = g0,
             alpha = 0.5,
             linewidth = 1,
             linetype = "solid",
             color = "steelblue") +
  # Facet by nest to emphasize group-level structure
  facet_wrap(facets =~ nest,
             nrow = 3,
             ncol = 3) +
  theme_bw()

## we estimated the variable intercept by group
# the intercept value depicts the level of negotiation per group 
# expected level of negotiation is drawn as dotted line for each nest 
# blue line is overall average of the data, estimated from the model 
# deviation of dotted line from blue line accounts for nest specific context to better
# account for food association 
# variation among nests is determined with a standard deviation 


# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: GLMM Exercise - `sleep` Data Set
# ============================================================

# ------------------------------------------------------------
# sleep dataset (built-in R dataset)
#
# This dataset contains measurements of increased sleep time
# after administering two different drugs.
#
# Structure:
# - 20 observations total
# - 10 subjects (each measured twice)
# - Paired / repeated-measures design
#
# Variables:
#   extra : Increase in hours of sleep compared to baseline
#   group : Indicates which drug was administered (factor with two levels ("1", "2"))
#   ID    : factor identifying each subject; Each subject appears once in each group
# ------------------------------------------------------------

# Q1 – Visualization:
# Compare sleep increase ("extra") between the two drug groups.

print(sleep)

sleep %>%
  mutate(group = as.numeric(group)) %>%
  ggplot(aes(x = group,
             y = extra,
             color = ID)) +
  geom_line() +
  geom_point() +
  facet_wrap(facet=~ID)

# Goals:
# - Show individual-level responses
# - Highlight paired structure (same subject in both groups)
# - Use color to identify subjects
# - facet by individual #Connect observations from the same subject using lines

# Q2 - Model development:

## choose a gaussian model 
## ID was the random effect column 


m_sleep <- glmmTMB(extra ~ group + (1 | ID),
        data = sleep,
        family = "gaussian")

m_sleep_wo_r <- glmmTMB(extra ~ group,
data = sleep,
family = "gaussian")

# Goal:
#   Examine how drug administration affects sleep duration.
#
# Key considerations:
#   - Response variable (extra) is continuous
#   - Drug (group) represents the treatment of interest
#   - Subject ID represents repeated measurements on the same
#     individuals

# ============================================================
# EXERCISE: GLMM Exercise - `grouseticks` Data Set
# ============================================================

library(lme4)
data("grouseticks")

# ------------------------------------------------------------
# grouseticks dataset (from lme4 package)
#
# This dataset contains counts of parasitic ticks
# collected from red grouse chicks across multiple years
# and locations in Scotland.
#
# Structure:
# - 403 observations
# - Repeated measurements across broods and years
# - Count data with hierarchical (nested) structure
#
# Variables:
#   TICKS : Number of ticks counted on a chick
#   YEAR  : Sampling year
#   HEIGHT: height above sea level (meters)
#   LOCATION : Sampling site (grouping variable)
#   INDEX : Observation-level identifier
#
# Key features of the dataset:
# - Response variable is count data
# - Observations are grouped by brood and year
# ------------------------------------------------------------

# Q1 – Visualization:
#
# Goal:
#   Examine the relationship between parasite load (ticks) at the brood level and height above sea level.
#
# Key considerations:
# - Calculate average tick counts for each brood
# - Plot mean ticks vs. height
# - Color points by sampling year

# Q2 – Model development:
#
# Goal:
#   Develop a model to examine the relationship between parasite load (ticks) at the brood level and height above sea level.
#
# Key considerations:
#   - Response variable (TICKS) is count
#   - HEIGHT represents the variable of interest
#   - BROOD represents a grouping factor of repeated measurements
#   - YEAR represents another grouping factor of repeated measurements
