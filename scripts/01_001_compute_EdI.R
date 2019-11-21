##' ---
##' title: "Compute EdI"
##' date: "`r format(Sys.time(), '%d %B, %Y')`"
##' author: Damien G.
##' licence: GPL-3
##' ---
##' 
##' ## Short Description --------------------------------------------------------------------------
##' 
##' HDI is based on life expectency it can't so it is a tricky index to be used to model health
##' based outcomes such as Cancer Mortality rates. To overcome this issue, we want to compute a
##' custom index similar to HDI that will only consider the Education and Standard of living 
##' components.
##' 
##' Our methodology is based on HDI calculation methodology from the dedicated (technical note)[http://hdr.undp.org/sites/default/files/hdr2015_technical_notes.pdf].
##' 
##' EdI is computed using the (2018_all_indicators.xlsx)[http://hdr.undp.org/sites/default/files/2018_all_indicators.xlsx] 
##' files downloaded from hdr.undp.org on November 19th 2019
##' 
##' The 4 indices of interest are:
##'   - Expected years of schooling (years) (id: 69706)
##'   - Gross national income (GNI) per capita (2011 PPP$) (id: 141706)
##'   - Life expectancy at birth (years) (id: 69206)
##'   - Mean years of schooling (years) (id: 103006)
##'   


##' ## Main script --------------------------------------------------------------------------

rm(list = ls())

## working directory and libraries
setwd("/mnt/data/georgesd/_PROJECTS/edi_globocan/ediglobocanhub/workdir/")
library('dplyr')
library('tidyr')
library('readxl')
library('readr')

out.dir <- paste0("../data")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

## load all indicators data
dat.undp <- readxl::read_xlsx('../data/undp/2018_all_indicators.xlsx')

## extract the expected years of schooling
dat.eyos <-
  dat.undp %>%
  filter(indicator_id %in% 69706) %>%
  select(iso3, country_name, eyos = `9999`)

## extract GNI  
dat.gni <-
  dat.undp %>%
  filter(indicator_id %in% 141706) %>%
  select(iso3, country_name, gni = `9999`)

## extract life expectancy at birth
dat.leab <-
  dat.undp %>%
  filter(indicator_id %in% 69206) %>%
  select(iso3, country_name, leab = `9999`)

## extract mean years of schooling
dat.myos <-
  dat.undp %>%
  filter(indicator_id %in% 103006) %>%
  select(iso3, country_name, myos = `9999`)

## compute EdI
dat.EdI <-
  dat.eyos %>%
  full_join(dat.gni, by = c('iso3', 'country_name')) %>%
  full_join(dat.leab, by = c('iso3', 'country_name')) %>%
  full_join(dat.myos, by = c('iso3', 'country_name')) %>%
  mutate(
    Health_index = (leab - 20) / (85 - 20),
    Mean_years_of_schooling_index = myos / 15,
    Expected_years_of_schooling_index = eyos / 18,
    Education_index = (Mean_years_of_schooling_index + Expected_years_of_schooling_index) / 2,
    Income_index = (log(gni) - log(100)) / (log(75000) - log(100)),
    HDI = round((Health_index * Education_index * Income_index)^(1/3), digits = 3),
    EdI = round((Education_index * Income_index)^(1/2), digits = 3)
  )

## check coutries without values
dat.EdI %>% filter(is.na(EdI)) %>% pull('country_name')
# [1] "Korea (Democratic People's Rep. of)" "Nauru" "San Marino"                         
# [4] "Tuvalu" "Somalia"   

## keep only usefull data for our analysis
dat.EdI.mini <-
  dat.EdI %>%
  filter(!is.na(EdI)) %>%
  select(iso3, country_name, EdI)

## save EdI table
write_csv(dat.EdI.mini, path = '../data/EdI.csv')
