#' DESCRIPTION:
#' Script for Unconstrained Ordination

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               GGally,
               vegan)

## PCA
## use iris data
df_iris <- iris %>%
  as_tibble() %>%
  janitor::clean_names()

df_iris %>%
  ggpairs(
    progress = FALSE,                         # Disable progress bar during plotting
    columns = c("sepal_length",               # Select the four numeric columns to include
                "sepal_width",
                "petal_length",
                "petal_width"),
    aes(
      color = species,                        # Color points by species
      alpha = 0.5)                             # Make points semi-transparent for better visibility 
      ) +
  theme_bw()                                  # Use a clean theme for better readability
  theme_bw()                                  # Use a clean theme for better readability

# Extract only the petal measurements from the iris dataset
  
  df_petal <- df_iris %>% 
 select(starts_with("petal_"))     # Select columns whose names start with 
                                   #"petal_" (petal length and petal width)
  
  obj_pca <- prcomp(               # Perform Principal Component Analysis (PCA) 
    x = df_petal,                  #on the selected petal data
    center = TRUE,                # Input data for PCA
    scale = TRUE)                 # Subtract the mean of each column so variables are centered at 0
                                  # Divide by the standard deviation of each column so variables have unit variance
  # obj_pca now contains:
  # - obj_pca$x        : PCA scores (coordinates of each sample in PC space)
  # - obj_pca$rotation : Eigenvectors (directions of principal components)
  # - obj_pca$sdev     : Standard deviations of the principal components (square = variance explained)
  
  print(obj_pca)
  
## get more information 
  summary(obj_pca)
  
## get PC axes values
  df_pca <- df_iris %>%
    bind_cols(obj_pca$x)
  
## draw a figure comparing PC1 values between species 
  df_pca %>%
    ggplot(aes(x = species,
               y = PC1)) +
    geom_boxplot()
  
  
## NMDS

## call sample data from vegan package
  data(dune)
  
#visual
  # Take the 'dune' dataset from vegan
  dune %>% 
    as_tibble() %>%         # Convert 'dune' from matrix/data frame to tibble for tidyverse compatibility
    select(1:3) %>%         # Select the first three species (columns) for visualization
    ggpairs() +             # Create a pairwise scatterplot matrix (plots all pairwise relationships)
    theme_bw()              # Apply a clean, minimal theme with a white background
  
  
  # Compute pairwise dissimilarities between sites using Bray-Curtis distance
  # column is species, row is site (or similar)
 
  m_bray <- vegdist(dune,               # Input community matrix (sites × species)
                 method = "bray")       # Use Bray-Curtis dissimilarity, commonly used in ecology
 
 obj_nmds <- metaMDS(comm = m_bray,
         k = 2)

#visualize NMDS
 data(dune.env)
  
 # Combine environmental data with NMDS coordinates
 df_nmds <- dune.env %>%           # Start with the environmental data for each site
   as_tibble() %>%                  # Convert to tibble for tidyverse-friendly operations
   bind_cols(obj_nmds$points) %>%   # Add NMDS coordinates (site scores) as new columns
   janitor::clean_names()           # Clean column names (lowercase, replace spaces/special characters)
 
 # Visualize NMDS site scores with points and 95% confidence ellipses by land-use intensity
 df_nmds %>% 
   ggplot(aes(
     x = mds1,           # NMDS axis 1 (first dimension from metaMDS)
     y = mds2,           # NMDS axis 2 (second dimension from metaMDS)
     color = use         # Color points by land-use intensity (or other grouping variable)
   )) +
   geom_point(size = 3) +               # Plot each site as a point, slightly larger for visibility
   stat_ellipse(level = 0.95,           # Draw 95% confidence ellipses around each group
                linetype = 2) +        # Dashed line for ellipse
   theme_bw() +                         # Apply a clean black-and-white theme
   labs(color = "Land-use intensity",   # Add a legend title
        x = "NMDS1",                    # Label x-axis
        y = "NMDS2")                    # Label y-axis
 
 # Perform PERMANOVA to test whether plant community composition differs by land-use intensity
 adonis2(
   m_bray ~ use,   # Model formula: Bray-Curtis dissimilarities (m_bray) explained by 'use' factor
   data = df_nmds  # Data frame containing the grouping variable 'use' aligned with the distance matrix
 )
 
# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: PCA using the iris dataset
# ============================================================

# In this exercise, you will perform a Principal Component
# Analysis (PCA) using all morphological measurements in the
# iris dataset and visualize multivariate trait patterns
# among species.

# 1. Using all four morphological variables
#    (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width),
#    perform a PCA.
  
   df_ps <- df_iris %>%
    select(-species)
  
   obj_pca4 <-  prcomp(
     x = df_ps,                 
    center = TRUE,                
    scale = TRUE)
  
  summary(obj_pca4)

# 2. Visualize flower morphology in PC axes whose cumulative
#    contribution exceeds 90%; color points by species.
  
  df_pca4 <- df_iris %>%
    bind_cols(obj_pca4$x)
  
  df_pca4 %>% 
    ggplot(aes(
      x = PC1,  
      y = PC2,       
      color = species 
    )) +
    geom_point() 

# 3. Which morphological traits contribute most strongly to
#    the first and second principal components? How?

# Petal length impacts PC1 the most, and Petal length impacts PC2 the most.
# I hope I interpreted that correctly (the greatest value under each PC column)
  
# ============================================================
# EXERCISE: NMDS using the BCI dataset
# ============================================================

# In this exercise, you will perform a Non-metric Multidimensional
# Scaling (NMDS) using the BCI tree community dataset and explore
# patterns in species composition among sites.

data("BCI", "BCI.env")

# 1. Using the BCI dataset, calculate a dissimilarity matrix
#    (e.g., Bray-Curtis) and perform NMDS.
  
  m_bray_bci <- vegdist(BCI,               
                    method = "bray") 
  
  obj_nmds_bci <- metaMDS(comm = m_bray_bci,
                      k = 2)
  
# 2. Visualize sites in NMDS space.
#    - How are sites positioned relative to each other?
#    - Color or shape points by environmental groups or site
#      characteristics of your choice.
  
 data(BCI.env)
 
 df_nmds_bci <- BCI.env %>%           
   as_tibble() %>%                  
   bind_cols(obj_nmds_bci$points) %>%  
   janitor::clean_names() 
 
 df_nmds_bci %>% 
   ggplot(aes(
     x = mds1,           
     y = mds2,          
     color = habitat         
   )) +
   geom_point(size = 3)  
 
# 3. Perform PERMANOVA to examine if communities are grouped
#    by the environmental variable you selected.
 
 adonis2(formula = m_bray_bci ~ habitat,   
   data = df_nmds_bci)
