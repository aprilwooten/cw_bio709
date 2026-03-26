# NOTE:
# When instructed to "test XXX", you must report the outcome as comments
# that clearly summarize the relevant statistical results (e.g., effect size,
# direction, significance, and interpretation).
# Providing code alone without documenting and interpreting the results
# in comments will result in point deductions.

# dataset 1 ---------------------------------------------------------------

link1 <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_insect_emergence.rds"
df_emg <- readRDS(url(link1, "rb"))

# This dataset ('df_emg') contains daily measurements of aquatic insect emergence
# from two wetland sites over a full calendar year (Jan 1–Dec 31).

# Data structure:
# t           : Day of the year (integer), where 1 = January 1 and 365 = December 31
# site        : Site identifier (factor), with "s1" and "s2" representing the two wetlands
# emergence   : Emergence flux of aquatic insects (g/day)

# Q1. Visualize seasonal patterns in emergence flux at both sites
#     (e.g., plot emergence vs. day of year, with separate lines or colors for each site).
#     [1 point]

pacman::p_load(tidyverse,
               ggeffects,
               mgcv)

df_emg %>%
      ggplot(aes(x = t, 
           y = emergence, 
           color = site)) +
  geom_line(alpha = 0.5) +
  labs(x = "Day of year", 
       y = "Emergence (g/day)") +
  theme_bw()

# Q2. Test whether emergence flux differs significantly between the two sites,
#     while appropriately accounting for seasonal variation
#     [4 points]

emg_gam <- gam(emergence ~ site + s(t), 
               data = df_emg)

summary(emg_gam)

## Interpretation: 
## The smooth term s(t) in the GAM depicts the nonlinear seasonal pattern
## of insect emergence at the two wetland sites throughout the calendar year.
## The effective degrees of freedom (8.66) indicates that the relationship is 
## flexible, rather than a straight line. This makes sense in the context of 
## insect emergence during the warmer seasons, and a decrease in emergence
## during the colder seasons.
## The low p-value of (2e-16) also indicates that the seasonal pattern
## is incredibly significant, meaning the smooth term explains a large
## portion of the variation in temp. that the linear term couldn't capture alone.

# dataset 2 ---------------------------------------------------------------

link2 <- "https://raw.githubusercontent.com/aterui/cw_bio709/master/data_fmt/data_lake_invert.rds"
df_inv <- readRDS(url(link2, "rb"))

# This dataset 'df_inv' contains 100 observations from 10 lakes.
# Within each lake, 10 plots were established, spaced ~500 m apart.
# At each plot, the following variables were measured:

# s          : Species richness of invertebrates associated with aquatic plants at each plot
# hb         : Standing biomass of invertebrates associated with aquatic plants at each plot
# prod       : Production rate of aquatic plants (macrophytes), measured as g/month
# substrate  : Median diameter of substrate materials (mm)
# cond       : Water electrical conductivity (µS/cm);
#              a proxy for ionized nutrient levels (higher values may indicate eutrophication)
# lake       : lake ID

# Researcher's hypothesis was that: 
# (a) conductivity influences the productivity of macrophyes.
# (b) macrophyte's production rate ('prod') dictates invertebrate biomass ('hb') through bottom-up effects
# (c) macrophyte's production rate ('prod') dictates invertebrate richness ('s') through bottom-up effects 

# Q1. Create a scatter plot of macrophyte production ('prod', y-axis)
#     versus water conductivity ('cond', x-axis), with points colored by lake identity.
#     [1 point]

df_inv %>%
  ggplot(aes(x = cond,
             y = prod,
             color = factor(lake))) +
  geom_point() + 
  theme_classic()

# Q2. Create a scatter plot of raw invertebrate biomass ('hb', y-axis)
#     versus macrophyte production ('prod', x-axis), with points colored by lake identity.
#     [1 point]

df_inv %>%
  ggplot(aes(x = prod, 
             y = hb, 
             color = factor(lake))) +
  geom_point() + 
  theme_classic()

# Q3. Create a scatter plot of "log-transformed" invertebrate biomass ('hb', y-axis)
#     versus macrophyte production ('prod', x-axis), with points colored by lake identity.
#     [1 point]

df_inv %>%
  ggplot(aes(x = prod,
             y = log(hb), 
             color = factor(lake))) +
  geom_point() + 
  theme_classic()

# Q4. Test hypothesis (a) by modeling macrophyte production while
#     statistically controlling for potential confounding variables ('substrate', 'lake').
#     [3 points]

pacman::p_load(tidyverse,
               lme4,
               glmmTMB)

hyp_a <- lmer(prod ~ cond + substrate + (1 | lake), 
              data = df_inv)
summary(hyp_a)

## Interpretation: 
## The conductivity and substrate variables are not signficant, indicating that 
## conductivity, indicating that they do not have a measurable impact
## on macrophyte production.  


# Q5. Test hypotheses (a–c) simultaneously using a unified modeling framework.
#     Based on the resulting statistical tests, determine whether the overarching
#     hypothesis (a–c, combined) is supported or rejected.
#     - Use appropriate probability distributions.
#     - Use variable transformation if appropriate given the data.
#     [4 points]

library(piecewiseSEM)

pacman::p_load(tidyverse,
               GGally,
               piecewiseSEM)

hyp_a <- lmer(prod ~ cond + substrate + (1 | lake), 
              data = df_inv)

hyp_b   <- lmer(log(hb) ~ prod + (1 | lake), 
                data = df_inv)

hyp_c    <- lmer(s ~ prod + (1 | lake), 
                 data = df_inv)

sem_mod <- psem(hyp_a, hyp_b, hyp_c)
summary(sem_mod)

##Interpretation:
## (a) is statistically significant with a p-value of (0.0)
## (b) is not statistically significant with a p-value of (0.17)
## (c) is statistically significant with a p-value of (0.0)
## The hypotheses are partially supported, with two of the three 
## paths returning statistical significance 


# dataset 3 ---------------------------------------------------------------

link3 <- "https://raw.githubusercontent.com/aterui/cw_bio709/master/data_fmt/nutrient.rds"
nutrient <- readRDS(url(link3, "rb"))

print(trees)

# This dataset ('trees') contains measurements of 31 felled black cherry trees.
# The three variables represent tree diameter, height, and timber volume.
# Note: the variable 'Girth' is actually the diameter measured at 4 ft 6 in above ground.

# Data structure:
# Girth   : Numeric, tree diameter in inches (mislabelled as girth)
# Height  : Numeric, tree height in feet
# Volume  : Numeric, timber volume in cubic feet

# Q1. Visualize relationships among tree diameter ('Girth'), height ('Height'),
#     and timber volume ('Volume') (e.g., using scatterplot matrix or pairwise scatter plots).
#     [1 point]

trees%>%
ggpairs(c("Girth", 
          "Height",
          "Volume"))

# Q2. Perform an appropriate ordination or dimension reduction method to 
#     summarize these three variables into fewer composite axes.
#     Then, identify and retain axes that explain meaningful variation in the original variables
#     [3 points]

pca_trees <- prcomp(
  trees[, c("Girth", "Height", "Volume")], scale. = TRUE)
summary(pca_trees)

##Interpretation:
## The PC1 variance is is >80%, meaning that it is a biologically important
## axis and should be retained. 

# Q3. If justified, test whether the retained axis (or axes) is significantly 
#     related to "nutrient"; 
#     skip regression if the ordination does not support meaningful interpretation.
#     [1 point]

pc1 <- pca_trees$x[,1]
pc_trees <- lm(pc1 ~ nutrient)
summary(pc_trees)

##Interpretation: 
## There is a strong and significant positive relationship between the 
## nutrient levels and PC1. The slope of 0.210 indicates that for each 
## unit increase in nutrient, PC1 increases by 0.210 units. The p-value
## is significant (7.70e-06) indicating that this relationship is significant.

# dataset 4 ---------------------------------------------------------------

df_nile <- dplyr::tibble(
  year = time(Nile), # observation year
  discharge = as.numeric(Nile) # discharge
)

df_sunspot <- dplyr::tibble(
  year = time(sunspot.year), # observation year
  sunspots = as.numeric(sunspot.year) # the number of sunspots
)

# These datasets contain:
# - df_nile    : Annual discharge of the Nile River (Nile dataset)
# - df_sunspot : Annual sunspot counts (sunspot.year dataset)

# Q1. Create a combined data frame aligning the observation years
#     (i.e., only include years present in both datasets)
#     [1 point]

library(dplyr)

df_nile_sun <- inner_join(
                df_nile, 
                df_sunspot, 
                by = "year")

# Q2. Test whether the number of sunspots is significantly related to Nile's discharge
#     [4 points]

mod_nile_sun <- 
  lm(discharge ~ sunspots, 
     data = df_nile_sun)
summary(mod_nile_sun)

## Interpretation 
## There is no meaningful relationship between sunspots and Nile discharge. 
## Even though there is a small negative effect (-0.05654),
## the result is not statistically significant (p = 0.887). 
## Sunspot activity does not predict Nile discharge. 
