#' DESCRIPTION:
#' Script for Constrained Ordination

# in-class ----------------------------------------------------------------


pacman::p_load(tidyverse,
               GGally,
               vegan)
# Load example datasets from the vegan package

data("varespec")  # Species abundance matrix: plant species measured at 44 sites
data("varechem")  # Environmental variables (soil chemistry) for the same sites

# Response matrix: species abundances
m_y <- varespec
colnames(m_y) <- str_to_lower(colnames(varespec))   # Convert column names to lowercase for consistency
                                                    # This avoids issues with formulas or later data manipulation
# Clean and prepare environmental data
# - Convert to tibble for easier use with dplyr
# - Use janitor::clean_names() to make column names lowercase, consistent, and safe for programming
df_env <- as_tibble(varechem) %>% 
  janitor::clean_names()

# Pairwise scatterplot matrix of the first 3 species columns in m_y
m_y %>%
  ggpairs(
    progress = FALSE,      # Disable the progress bar during plotting
    columns = 1:3,         # Only use the first 3 species for the plot
    aes(
      alpha = 0.5           # Make points semi-transparent to reduce overplotting
    )
  ) +
  theme_bw()               # Use a clean theme with white background

# Run Redundancy Analysis (RDA)
# - Response: m_y (species abundance matrix, varespec)
#   Note: Hellinger transformation is recommended for raw species counts.
# - Predictors: n, p, ca (numeric soil chemistry variables from df_env)
#   These represent Nitrogen (n), Phosphorus (p), and Calcium (ca) levels.
# - Data: df_env (cleaned environmental variables)
# The parentheses around obj_rda print the RDA summary immediately
(obj_rda <- rda(m_y ~ n + p + ca,
                data = df_env))

# Perform permutation test on the RDA object
# - obj_rda: RDA object containing species (response) and soil chemistry (predictors)
# - by = "margin": test the marginal effect of each predictor
#   This evaluates how much variation each predictor explains **after accounting for all other predictors**
# - permutations = 999: number of random permutations used to generate the null distribution for the F-statistic
anova.cca(obj_rda, 
          by = "margin", 
          permutations = 999)
# Extract site (sample) scores from the RDA object
# - display = "sites": returns the coordinates of sampling sites in RDA space
# - scaling = 2: correlation scaling, which emphasizes relationships between sites
#   and explanatory variables (i.e., how site positions correlate with predictors)
# The site scores are then combined with the environmental data for visualization

df_rda <- scores(obj_rda, 
                 display = "sites",
                 scaling = 2) %>% 
  bind_cols(df_env) %>%           # append environmental variables (e.g., soil chemistry)
  janitor::clean_names()          # standardize column names for tidy workflows

# RDA vectors for environmental predictors
df_bp <- scores(obj_rda, 
                display = "bp", 
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

# Create a ggplot2 ordination plot
# - Points represent sites positioned by their constrained community composition
# - Color gradient reflects the nitrogen (n) concentration at each site
df_rda %>% 
  ggplot(aes(x = rda1,
             y = rda2)) +        # color sites by nitrogen level
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = rda1 * 10, # 10 is arbitrary scaling for visualization
                   y = 0, yend = rda2 * 10),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = rda1 * 10.5,    # slightly beyond arrow tip
                y = rda2 * 10.5,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "RDA1",
       y = "RDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()

# Perform distance-based Redundancy Analysis (dbRDA)
# - Response: m_y (species abundance matrix)
# - Predictors: n, p, ca (environmental variables: nitrogen, phosphorus, potassium, calcium)
# - data = df_env: provides the environmental variables used as predictors
# - distance = "bray": computes a Bray-Curtis dissimilarity matrix from the raw species abundances
#   internally, before performing the constrained ordination
# The result is a dbrda object representing the constrained ordination in distance space
(obj_db <- dbrda(m_y ~ n + p + ca,
                 data = df_env,
                 distance = "bray"))

# Perform a permutation test on the dbRDA object
# - obj_db: the dbRDA (dbrda) object containing species distances and predictors
# - by = "margin": tests the **marginal effect** of each predictor, i.e., the unique contribution
#   of a variable after accounting for all other predictors in the model
# - permutations = 999: randomly permutes the data 999 times to generate a null distribution
#   for the F-statistic, allowing assessment of statistical significance
anova.cca(obj_db,
          by = "margin",
          permutations = 999)

# Extract site scores from the dbRDA object
df_db <- scores(obj_db, 
                display = "sites",
                scaling = 2) %>% 
  as_tibble() %>%              
  bind_cols(df_env) %>%        
  janitor::clean_names()       

# dbRDA vectors for environmental predictors
df_bp <- scores(obj_db, 
                display = "bp", 
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

# Create a ggplot2 ordination plot
df_db %>% 
  ggplot(aes(x = db_rda1,
             y = db_rda2)) +        # color sites by nitrogen level
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = db_rda1,
                   y = 0, yend = db_rda2),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = db_rda1 * 1.1,    # slightly beyond arrow tip
                y = db_rda2 * 1.1,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "dbRDA1",
       y = "dbRDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()


# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Community Ordination and Environmental Gradients
# ============================================================

library(vegan)
data("mite", "mite.env")

# The mite datasets contain information on Oribatid mite communities
# sampled from a small peatland area (2.5 m × 10 m).
#
# There are linked datasets:
# ------------------------------------------------------------
# mite     : Species abundance data (35 mite species × 70 sites)
# mite.env : Environmental variables measured at the same sites
# ------------------------------------------------------------
#
# Environmental variable descriptions (mite.env):
# ------------------------------------------------------------
# SubsDens : Substrate density (g/L)
# WatrCont : Water content of the substrate (g/L)
# Substrate: Substrate type (factor with multiple levels)
# Shrub    : Shrub density (ordered factor: low → high)
# Topo     : Microtopography (Blanket vs Hummock)
# ------------------------------------------------------------

# 1. Explore and visualize interrelationships among species abundances.
#    - Examine patterns of co-occurrence.
#    - Assess whether relationships among species appear linear or nonlinear.

ggpairs(wisconsin(mite),
        columns = 1:5, 
        progress = FALSE)

# 2. Fit a redundancy analysis (RDA) model using environmental variables of your choice.
#    - Visualize the ordination results.
#    - Examine gradients and species–environment relationships.
#    - Evaluate whether the assumptions of RDA are appropriate for these data.

(mite_rda <- rda(mite ~ SubsDens + WatrCont + Substrate,
                 data = mite.env))

df_env <- mite.env %>%
  janitor::clean_names() %>%
  mutate(num_shurb = as.numeric(shrub))



# 3. Apply alternative ordination methods.
#    - Canonical correspondence analysis (CCA; see ?cca()).
#    - Distance-based RDA (dbRDA).

m_y <- mite

obj_rda_mite <- rda(m_y ~ subs_dens + watr_cont + shrub,
                    data = df_env)

scores(obj_rda_mite,
       display = "sites",
       scaling = 2) %>%
  bind_cols(df_env) %>%
  as_tibble() %>%
  janitor::clean_names () %>%
  ggplot(aes(x = rda1, 
             y = rda2)) +
  geom_point()

obj_db_mite <- dbrda(m_y ~ subs_dens + watr_cont + shrub,
                     data = df_env,
                     distance = "bray")

obj_cca_mite <- cca(m_y ~ subs_dens + watr_cont,
                    data = df_env)

scores(obj_cca_mite,
       display = "sites",
       scaling = 2) %>%
  bind_cols(df_env) %>%
  as_tibble() %>%
  janitor::clean_names () %>%
  ggplot(aes(x = cca1, 
             y = cca2)) +
  geom_point()

# 4. Compare RDA, CCA, and dbRDA.
#    - Perform permutation analysis to examine the significance of predictor variables
#    - Discuss which method is most appropriate for these data and why.

anova.cca(obj_cca_mite,
          by = "margin",
          permutations = 999)

anova.cca(obj_db_mite,
          by = "margin",
          permutations = 999)

anova.cca(obj_rda_mite,
          by = "margin",
          permutations = 999)

# try variable transformation 
m_mite_trans <- vegan::wisconsin(mite)
