#' DESCRIPTION:
#' Script for time-series

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               forecast,
               lterdatasampler,
               daymetr,
               glarma)

url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_ts_anormaly.csv"
(df_ts <- read_csv(url))

# Plot time-series anomalies
df_ts %>% 
  ggplot(aes(x = year,          # Map 'year' to the x-axis
             y = anormaly)) +  # Map 'anormaly' to the y-axis
  geom_line() +                  # Add a line connecting the points
  geom_point() +                 # Add points at each observation
  theme_bw() +                   # Use a clean black-and-white theme
  labs(
    x = "Year",                     # Label x-axis
    y = "Anomaly"                   # Label y-axis
  )


# Fit a simple linear model: anomaly as a function of year
m_ts <- lm(anormaly ~ year, 
           data = df_ts)

# Show the model summary
summary(m_ts)

## figure 
df_ts %>% 
  ggplot(aes(x = year, 
             y = anormaly)) +  # Set up the plot: x = year, y = anomaly
  geom_line(linetype = "dotted") +     # Draw a dotted line connecting the points
  geom_point(alpha = 0.25) +           # Add points with transparency (alpha = 0.25)
  geom_abline(intercept = coef(m_ts)[1],  # Add regression line from linear model
              slope = coef(m_ts)[2]) +
  theme_bw()                             # Use a clean black-and-white theme

y <- NULL
y[1] <- 0

for (i in 1:99) {
  y[ i + 1] <- y[i] + rnorm(1, mean = 0, sd = 1)}

tibble(y = y,
       x = 1:length(y)) %>%
  ggplot(aes(x = x, 
             y = y)) +
  geom_point() +
  geom_line()

# Create a tibble (modern data frame) from the LakeHuron time series
df_huron <- tibble(
  year = time(LakeHuron),                # Extracts the time component (years) from the LakeHuron ts object
  water_level = as.numeric(LakeHuron)    # Converts LakeHuron values to numeric (from ts class)
) %>% 
  arrange(year)                           # Ensures the data is ordered by year

# Plot Lake Huron time series with a linear trend
df_huron %>% 
  ggplot(aes(x = year, y = water_level)) +
  geom_point(alpha = 0.25) +       # Semi-transparent points
  geom_line(linetype = "dotted") + # Dotted line connecting points
  geom_smooth(method = "lm",       # Linear trend line
              color = "black",
              linewidth = 0.5) +
  theme_bw() +
  labs(x = "Year", y = "Water Level")

## MAR Model 
(m_ar1 <- Arima(
  df_huron$water_level,       # The time series data we want to model
  order = c(0, 0, 1)          # ARIMA model orders: c(p, d, q)
  ))
## the first number reprsents what order you want to use in AR process

## fitted values
df_huron_ar1 <- df_huron %>%
  mutate(fit = fitted(m_ar1) %>%
           as.numeric())

df_huron_ar1 %>%
  ggplot() +
  geom_point(aes( x = year,
                  y = water_level),
             alpha = 0.25) +
  geom_line(aes(x = year,
                y = fit),
            color = "steelblue") +
  theme_bw()


## ARMA Model
(m_ma1 <- Arima(
  df_huron$water_level,       # The time series data we want to model
  order = c(1, 0, 1)          # ARIMA model orders: c(p, d, q)
))

## ARIMA Model 
(m_arima <- Arima(df_huron$water_level,
                  order = c(1, 1, 0)))

## model selection 
auto.arima(
  y = df_huron$water_level,
  stepwise = FALSE, 
  ic = "aic")

## ARIMAX 
data("ntl_icecover")

#convert to tibble and filter for a specific lake and year >= 1980
df_ice <- ntl_icecover %>%
  as_tibble() %>%
  filter(between(year, 1980, 2014),
         lakeid == "Lake Mendota") %>%
  arrange(year)

#create a line + point plot of ice duration over time 
df_ice %>%
  ggplot(aes(x = year,
             y = ice_duration)) +
  geom_line(linetype = "dotted") +
  geom_point(alpha = 0.25) +
  theme_bw() +
  labs (y = "Ice duration",
        x = "Year")

# Download daily climate data from Daymet for Lake Mendota
list_mendota <- download_daymet(
  site = "Lake_Mendota",   # Arbitrary name you assign to this site
  lat = 43.1,              # Latitude of the lake
  lon = -89.4,             # Longitude of the lake
  start = 1980,            # Start year
  end = 2024,              # End year
  internal = TRUE          # Return the data as an R object rather than saving to disk
)

df_temp <- list_mendota$data %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(
    date = as.Date(paste(year, yday, sep = "-"),
                   format = "%Y-%j"),
    month = month(date)
  ) %>%
  arrange(year, yday) %>%
  group_by(year) %>%
  summarize(temp_min = round(mean(tmin_deg_c), 2))

df_ice <- df_ice %>%
  left_join(df_temp, by = "year")

(obj_arima <- auto.arima(
  df_ice$ice_duration,
  xreg = df_ice$temp_min,
  stepwise = FALSE
))

df_ice %>%
  ggplot(aes(x = temp_min,
             y = ice_duration)) +
  geom_point()

confint(obj_arima, level = 0.95)

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Bison Body Mass, Climate, and Time-Series Analysis
# ============================================================

library(lterdatasampler)

# The "knz_bison" dataset contains long-term monitoring data
# on bison captured at Konza Prairie Biological Station.
#
# ------------------------------------------------------------
# Key columns may include:
# rec_year      : Year of capture
# animal_sex    : Sex of the individual (e.g., female, male)
# animal_weight : Body mass of bison
# ------------------------------------------------------------
#
# In this exercise, you will explore long-term trends in bison
# body mass and evaluate how climate variability may influence
# weight dynamics over time.

# 1. Explore the structure of the knz_bison dataset.
#    - Inspect variable types and missing values.
#    - Reformat variables as needed for analysis.

sapply(knz_bison,
       FUN = function(x) {
         sum(is.na(x))
       })

df_knz <- knz_bison %>% 
  drop_na(animal_weight)

# 2. Subset the data to include observations from 1994–2012.

df_knz <- df_knz %>% 
  filter(between(rec_year, 1994, 2012)) %>% 
  rename(year = rec_year)


# 3. Calculate the average body mass for female and male bison
#    for each year in the selected time period.

df_mu <- df_knz %>% 
  group_by(year,
           animal_sex) %>% 
  summarize(w = mean(animal_weight),
            .groups = "drop")

# 4. Obtain climate data from the daymetr dataset.
#    - Identify relevant climate variables (e.g., temperature,
#      precipitation).
#    - Associate climate data with knz_bison by year.
#    - Coordinates: Lat 39.09300	Lon -96.57500

df_clim_knz <- download_daymet(
  site = "Konza",
  lat = 39.09300,
  lon = -96.57500,
  start = 1994,            # Start year
  end = 2012,              # End year
  internal = TRUE          # Return the data as an R object rather than saving to disk
) %>% 
  {.[["data"]]} %>% 
  janitor::clean_names() %>% 
  as_tibble() %>% 
  mutate(
    date = as.Date(paste(year, yday, sep = "-"),
                   format = "%Y-%j")
  ) %>% 
  relocate(date)

df_clim_mu <- df_clim_knz %>% 
  group_by(year) %>% 
  summarize(cprcp = sum(prcp_mm_day))

df_mu <- df_mu %>% 
  left_join(df_clim_mu,
            by = "year")

df_mu %>% 
  ggplot(
    aes(x = year,
        y = w,
        color = animal_sex)
  ) +
  geom_point()

# 5. Perform a time-series analysis to examine whether selected
#    climate variables influence annual bison body mass.
#    - Consider temporal autocorrelation and lag effects.
#    - Model males and females separately

m_male <- df_mu %>% 
  filter(animal_sex == "M") %>% 
  arrange(year) %>% 
  { auto.arima(y = .$w,
               xreg = .$cprcp,
               stepwise = FALSE, 
               d = 0) }

confint(m_male)

m_female <- df_mu %>% 
  filter(animal_sex == "F") %>% 
  arrange(year) %>% 
  { auto.arima(y = .$w,
               xreg = .$cprcp,
               stepwise = FALSE, 
               d = 0) }

confint(m_female)

# 6. Using your fitted model, compare observed bison body mass
#    with predicted values for the period 2014–2020.
#    - Evaluate model performance and discuss sources of uncertainty.


