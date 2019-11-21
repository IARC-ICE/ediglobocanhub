##' ---
##' title: "Compute globocan ASR"
##' date: "`r format(Sys.time(), '%d %B, %Y')`"
##' author: Damien G.
##' licence: GPL-3
##' ---
##' 
##' ## Short Description ---------------------------------------------------------------------------
##' Here we compute from 5 years age group cancer incidence (resp. mortality) data per cancer 
##' type, country and sex the per 100.000 person.year age-standardized rates (ASR).
##' 
##' Raw data can be extracted from [GCO cancer today website](https://gco.iarc.fr/today/home).
##' 
##' Data extraction has been done on September 5th 2018.
##' 

##' ## Main script ---------------------------------------------------------------------------------

rm(list = ls())

## working directory and libraries
setwd("/mnt/data/georgesd/_PROJECTS/edi_globocan/ediglobocanhub/workdir/")
library('readr')
library('dplyr')
library('tidyr')
library('rjson')

out.dir <- paste0("../data")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

## load 5 years cancer incidence and mortality data per cancer site, country and sex
dat.im <- read_csv('../data/globocan/IMP_05092018_IM.csv') %>% rename_all(tolower)
## load 5 years population data per country and sex
dat.pop <- read_csv('../data/globocan/IMP_05092018_Pops.csv') %>% rename_all(tolower)
## load the country code correspondance table
dat.count <- 
  fromJSON(file = '../data/globocan/GLOBOCAN_COUNTRY_REF.json') %>% 
  purrr::map_df(
    ~ .x %>% unlist() %>% t() %>% as_tibble()
  )
dat.count.valid.name <-
  read_xlsx('../data/globocan/DataEntryTemplateCountryLevel.xlsx', sheet = 'Data') %>% 
  rename_all(tolower)

## load the cancer code correspondance table
dat.corcanc <- read_xlsx('../data/globocan/cancer_corresp_table.xlsx') %>% rename_all(tolower)
## load the Segi (1960) and modified by Doll and al. (1966) standard population table
dat.sw <- 
  read_csv("../data/globocan/ASR_population.csv") %>%
  mutate(ageg = paste0('w', sub('-', '_', sub('\\+', '', ageg)))) %>%
  spread(key = ageg, value = wi)

## combine tables and compute ASR
dat.asr <-
  dat.im %>% 
  rename(cancer_id = cancer, ntotal = total) %>%
  ## add cancer label (globocan 2012 style)
  left_join(dat.corcanc, by = 'cancer_id') %>%
  mutate(cancer = cancer_globocan_2012) %>%
  ## group the cancer the globocan 2018 cancer types to match 2012 nomenclature 
  dplyr::select(country, type, sex, cancer, ntotal, n0_4:n85)%>%
  tidyr::gather(key = ageg, value = n, ntotal, n0_4:n85) %>%
  group_by(country, type, sex, cancer, ageg) %>%
  summarise(n = sum(n, na.rm = TRUE)) %>%
  ungroup() %>%
  spread(ageg, n) %>%
  ## add population demographic data
  left_join(dat.pop %>% rename(ptotal = total), by = c("country", "sex")) %>%
  ## add countries label
  left_join(
    dat.count %>% 
      dplyr::select(country, iso_3_code) %>% 
      mutate(country = as.numeric(country)),
    by = 'country'
  ) %>%
  left_join(
    dat.count.valid.name %>% 
      dplyr::select(label = cntry_terr, iso_3_code),
    by = 'iso_3_code'
  ) %>%
  rename(country_id = country, country = label, iso3 = iso_3_code) %>%
  ## remove aggregated countries
  filter(country_id < 900) %>%
  ## compute the ASR
  mutate(
    ASR = 
      n0_4   / p0_4   * dat.sw$w0_4   +
      n5_9   / p5_9   * dat.sw$w5_9   +
      n10_14 / p10_14 * dat.sw$w10_14 +
      n15_19 / p15_19 * dat.sw$w15_19 +
      n20_24 / p20_24 * dat.sw$w20_24 +
      n25_29 / p25_29 * dat.sw$w25_29 +
      n30_34 / p30_34 * dat.sw$w30_34 +
      n35_39 / p35_39 * dat.sw$w35_39 +
      n40_44 / p40_44 * dat.sw$w40_44 +
      n45_49 / p45_49 * dat.sw$w45_49 +
      n50_54 / p50_54 * dat.sw$w50_54 +
      n55_59 / p55_59 * dat.sw$w55_59 +
      n60_64 / p60_64 * dat.sw$w60_64 +
      n65_69 / p65_69 * dat.sw$w65_69 +
      n70_74 / p70_74 * dat.sw$w70_74 +
      n75_79 / p75_79 * dat.sw$w75_79 +
      n80_84 / p80_84 * dat.sw$w80_84 +
      n85    / p85    * dat.sw$w85,
    p_tot = 
      p0_4 +
      p5_9 +
      p10_14 +
      p15_19 +
      p20_24 +
      p25_29 +
      p30_34 +
      p35_39 +
      p40_44 +
      p45_49 +
      p50_54 +
      p55_59 +
      p60_64 +
      p65_69 +
      p70_74 +
      p75_79 +
      p80_84 +
      p85
  ) %>%
  ## recorde sex and type as a factorial variable, remove useless space in country names
  mutate(
    sex = factor(sex, levels = c(0, 2, 1), labels = c("Both sexes", "Female", "Male")),
    type = plyr::mapvalues(type, from = c(0, 1), to = c("cancer incidence", "cancer mortality")),
    country = sub(" +$", "", country)
  )

## keep only usefull data for our analysis
dat.asr.mini <-
  dat.asr %>%
  filter(!is.na(ASR)) %>%
  select(iso3, country, sex, cancer, type, ptotal, ntotal, ASR)

## save ASR from globocan table
write_csv(dat.asr.mini, path = '../data/ASR_globocan.csv')



