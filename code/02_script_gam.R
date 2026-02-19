#' DESCRIPTION:
#' Script for GAMs

rm(list = ls())

# in-class ----------------------------------------------------------------
# we have used linear model up until this point
# using a straight line to describe x and y 
# GAM is an extension to a non-linear model, describes non-linear
# additives in the data 

pacman::p_load(tidyverse,
               ggeffects,
               mgcv)

link <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_water_temp.csv"

(df_wt_raw <- read_csv(link))

## check data format 
sapply(df_wt_raw, class)

# date time is characterized as a character, which is not great 
# need to modify to the date time column 

# Start from the raw water-temperature dataset
(df_wt <- df_wt_raw %>% 
    mutate(date = as.Date(date_time, 
                     format = "%m/%d/%Y"),    # Convert the character datetime column into a Date object.
                     year = year(date),       # The format argument specifies month/day/year.
                    month = month(date)) %>%  # Extract the calendar year from the Date object
                    filter(year == 2022,
                    between(month,
                    left = 3, right = 10)))   # Extract the month (1–12) from the Date object
     
                                              # Subset the data to observations from the year 2022 only
                                              # March through October

# calculating the mean and average temperature by date, time,and site
# data came from two sites (peabody park and near residential buildings)
# wanted to compare daily average temp between two sites (every 15 mins)
# summarize by daily average use group_by date and site, combo of two items
# specify which datasets are summarized 

## summarize by date and site 
df_wt_daily <- df_wt %>%
  group_by(date,
           site) %>%
  summarize(temp = mean(temp, na.rm = TRUE) %>%
  round(3), # 3 is how many decimal points you want to keep 
  .groups = "drop")

## visualize data 
df_wt_daily %>%
  ggplot(aes(x = date,
             y = temp,
             color = site)) +
  geom_point(alpha = 0.25) +
  theme_bw() +
  labs(x = "Date",
       y = "Water Temperature",
       color = "Wetland Type")

# this pattern is not linear, on the graph
# seasonal data is one of the most difficult to deal with
# because it is nonlinear 
# we are interested in seasonal pattern in water temperature 

# Add new variables and ensure proper data types
df_wt_daily <- df_wt_daily %>% 
  mutate(j_date = yday(date),           # Convert the date column to Julian day (day of year)
  site = factor(site))                  # Useful for modeling seasonal trends
    
                                        # Ensure 'site' is treated as a factor (categorical variable)
  
# describe how water temperature is different between sites and season
# control when the water temp was measured first before comparing btwn two sites
# linear modeling approach isnt appropriate, but we can check 

# Fit a Generalized Linear Model (GLM) with Gaussian family
m_glm <- glm(
  temp ~ j_date + site,  # Model daily temperature as a function of Julian day and wetland site
  data = df_wt_daily,    # Dataset to use
  family = "gaussian")   # Specify Gaussian distribution (standard linear regression)


summary(m_glm)
# strong effect of julian date
# since linear, this assumes julian date has monotonic influence 
# model prediction would suggest water temp is higher in october, which is not 
# the case 
# water temp is lower in woody than in open canpoy, which makes sense 
# but there was no stat sig


# Generate model predictions across all Julian days and wetland sites
df_pred <- ggpredict(m_glm,
                     terms = c("j_date [all]",    # Use all observed values of Julian day
                         "site [all]" )) %>%      # Generate predictions for all levels of the factor 'site'
                   as_tibble() %>%
  rename(j_date = x, 
         site = group)                            # 'group' from ggpredict() corresponds to the factor variable 'site'
                                                  # 'x' from ggpredict() corresponds to the predictor 'j_date'
          

# Rename the default columns to match the original dataset
  

# Plot daily water temperature and overlay model predictions
df_wt_daily %>% 
  ggplot(aes(
    x = j_date,   # Julian day on x-axis
    y = temp,     # Observed daily temperature on y-axis
    color = site  # Color points by wetland type (factor)
  )) +
  geom_point(alpha = 0.25) +
  # Overlay predicted values from the model
  # df_pred contains predictions from ggpredict()
  # aes(y = predicted) maps the model's predicted temperature to y
  geom_line(data = df_pred,
            aes(y = predicted)) +
  theme_bw() +
  labs(x = "Julian Date",         # x-axis label
       y = "Water Temperature",   # y-axis label
       color = "Wetland Type"     # Legend title for site color
  )

## GAM application - gam() is from mgcv package 
m_gam <- gam(temp ~ site + s(j_date),
             data = df_wt_daily,
             family = "gaussian")

summary(m_gam)

## visualize GAM predicton 
df_pred_gam <- ggpredict(m_gam,
                         terms = c(
                           "j_date [all]", 
                           "site [all]")
) %>% 
  as_tibble() %>%
  rename(site = group,
         j_date = x)

df_wt_daily %>% 
  ggplot(aes(
    x = j_date,
    y = temp, 
    color = site
  )) +
  geom_point(alpha = 0.25) +
  # Overlay predicted values from the GAM
  geom_line(data = df_pred_gam,
            aes(y = predicted)) +
  theme_bw() +
  labs(x = "Julian Date",         # x-axis label
       y = "Water Temperature",   # y-axis label
       color = "Wetland Type"     # Legend title for site color
  )

# lab ---------------------------------------------------------------------

# 1. Read directly from the raw GitHub URL
url <- "https://raw.githubusercontent.com/aterui/public-proj_restore-aqua-complex/v.1.0/data_raw/data_bat.csv"

# Try reading normally
df_bat <- read_csv("https://raw.githubusercontent.com/aterui/public-proj_restore-aqua-complex/v.1.0/data_raw/data_bat.csv",
                   show_col_types = FALSE)



# ============================================================
# DATA GUIDE: Bat Detector Data
# ============================================================

# ----------------------------
# Raw data columns
# ----------------------------

# Site
#   Location where bat detectors are deployed.
#   Levels:
#     "RECCON"  = prairie site without wetland
#     "RECWET"  = prairie site with constructed wetland
#     "WOODCON" = woody site without wetland
#     "WOODWET" = woody site with constructed wetland

# DATE
#   Calendar date of each bat pass record.
#   Expected format: YYYY-MM-DD (verify and standardize).

# TIME
#   Time of bat pass detection.
#   Expected format: HH:MM:SS (verify and standardize).

# AUTO ID*
#   Automatically identified bat species.
#   Species IDs may contain misclassifications or unknown labels
#   that should be carefully reviewed during data cleaning.

# ============================================================
# GOAL 1: Clean data
# ============================================================

# 1. Format column names
#   - Convert column names to a clean format

df_bat <- read_csv("https://raw.githubusercontent.com/aterui/public-proj_restore-aqua-complex/v.1.0/data_raw/data_bat.csv",
                   show_col_types = FALSE) %>%
  janitor::clean_names()

# 2. Examine each column carefully
#   - Check for missing values, inconsistent formats, and typos
#   - Confirm DATE and TIME are properly parsed as date/time objects
#   - Inspect AUTO ID values for NA
#   - Remove or correct invalid or unusable records as needed

sapply(df_bat, FUN = function(x) sum(is.na(x)))

# New derived columns to create:
# Site-level categories:
#   Prairie sites: "RECCON", "RECWET"
#   Woody sites:   "WOODCON", "WOODWET"

# 3. habitat_type
#   Broad site classification:
#     "prairie" = RECCON, RECWET
#     "woody"   = WOODCON, WOODWET

# 4. wetland_status
#   Presence/absence of wetland:
#     "no_wetland" = RECCON, WOODCON
#     "wetland"    = RECWET, WOODWET

df_bat <- df_bat %>%
  mutate(date = as.Date(date, format = "%m%d%Y"),
         habitat_type = case_when(site %in% c("RECCON", "RECWET") ~ "prairie",
                                  site %in% c("WOODCON", "WOODWET") ~ "woody"),
         wetland_status = case_when(site %in% c("RECCON", "WOODCON") ~ "no_wetland",
                                    site %in% c("RECWET", "WOODWET") ~ "wetland")) %>%
  drop_na(auto_id)


# ============================================================
# GOAL 2: Visualize daily bat activity
# ============================================================

# Objective:
#   Quantify and visualize bat activity as the number of bat passes per day.

# Steps:
#   - Aggregate data to calculate daily bat passes
#   - Convert DATE to Julian date
#   - Plot number of bat passes as a function of Julian date
#   - Optionally:
#       * Color or facet plots by site
#       * Smooth trends to visualize seasonal patterns

df_n <- df_bat %>%
  group_by(date,
           site,
           habitat_type,       # to keep this column's info
           wetland_status)%>%  # to keep this column's info
  summarize(pass = n(),
            .groups = "drop") %>%
  mutate(month = month(date),
        year = year(date),
         j_date = yday(date))
   

## plot by habitat type and wetland status 
df_n %>%
  ggplot(aes(x = date,
             y = pass,
             color = wetland_status)) +
  geom_point() +
  facet_wrap(facets =~ habitat_type) +
  theme_bw()


            

# ============================================================
# GOAL 3: Model differences among sites
# ============================================================

# Objective:
#   Test whether bat activity differs among the four detector sites.
#   Does the presence of wetland affect bat activity?
#   Is the effect of wetland is site-dependent?

# Modeling considerations:
#   - Response variable: daily bat passes
#   - Predictors may include:
#       * habitat_type
#       * wetland_status
#       * site (four-level factor)
#       * Julian date (to account for seasonality)
#   - Consider appropriate count models

mean(df_n$pass)
var(df_n$pass)

m_bat<- gam(pass ~ habitat_type + wetland_status + s(j_date),
    family = "nb",
    data = df_n)

summary(m_bat)