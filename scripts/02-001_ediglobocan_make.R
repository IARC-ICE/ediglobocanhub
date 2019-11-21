##' ---
##' title: "Profiling global cancer incidence and mortality by socioeconomic 
##' development main script"
##' date: "`r format(Sys.time(), '%d %B, %Y')`"
##' author: Damien G.
##' licence: GPL-3
##' ---

##' ## Short description 
##' 
##' Here is the main script tpo produce the "Profiling global cancer incidence and mortality by
##' socioeconomic development" paper outputs.
##' 
##' Results should be fully reproducible.
##' 
##' The workflow is fully handled through drake package.
##' Sourcing this file should produce all the tables and figures from the article.
##' 
##' Note: fig3 has to be manually edit to be in the final version to add group label in the legend
##' key wich is a pain to be coded
##' 
##' The main steps are :
##'  1. read.dat.plan: read cancer incidence and morality ASR and EdI data
##'  2. reshape.dat.plan: combine and reshape data
##'  3. rank.cancer.plan: order cancer types according to their burden
##'  4. gam.fit.plan: produce models of cancer incidence/mortality ASR vs EdI per cancer type
##'     and sex
##'  5. tab.prod.plan: produce tabS1 and tabS2
##'  6. scatter.plot.plan: produce fig1, fig2, figS1 and figS2
##'  7. gam.stack.plot.plan: produce fig3
##'  

##' ## Main script 

## 0. setup ----
rm(list = ls())

## define the path to project and define a working directory
path.to.hub <- '/mnt/data/georgesd/_PROJECTS/edi_globocan/ediglobocanhub/'
path.to.wd <- file.path(path.to.hub, 'workdir')
setwd(path.to.wd)

## get the same version of R packages that have been used for the analysis 
## it takes a while at the first execution
if(!require("checkpoint")) {install.packages("checkpoint"); require("checkpoint")}
checkpoint(
  "2019-11-21", 
  project = '../scripts', 
  R.version = '3.6.1', 
  checkpointLocation = '.',
  forceProject = FALSE
)

## define a seed for reproducibility
set.seed(21112019)

## require drake package and load the main functions contains in 02-002_ediglobocan_functions.R
require('drake')
require('visNetwork')
source('../scripts/02-002_ediglobocan_functions.R')

## define the output version to keep a track of previous versions
out.version <- 'lgh-001'
dir.create(file.path('../results', out.version), recursive = TRUE, showWarnings = FALSE)

## 1. read.dat.plan: read cancer incidence and morality ASR and EdI data ----
read.dat.plan <-
  drake_plan(
    dat.asr = read_csv(file_in('../data/ASR_globocan.csv')),
    dat.edi = read_csv(file_in('../data/EdI.csv')),
    countries.to.highlight = read_csv(file_in("../data/countries_to_highlight.csv"))
  )

## 2. reshape.dat.plan: combine and reshape data ----
reshape.dat.plan <-
  drake_plan(
    edi.breaks = setNames(c(0, 0.55, 0.7, 0.8), c('Low', 'Medium', 'High', 'Very High')),
    cancer.grp.1m = c("Prostate", "Breast", "Colorectum","Melanoma of skin", "Thyroid", "Kidney"), 
    cancer.grp.1s = c("Bladder", "Corpus uteri", "Non-Hodgkin lymphoma", "Testis", "Leukaemia", "Multiple myeloma", "Hodgkin lymphoma"),
    cancer.grp.2m = c("Lung", "Pancreas", "Brain, nervous system"), 
    cancer.grp.2s = c(), 
    cancer.grp.3m = c("Cervix uteri", "Liver", "Stomach"),
    cancer.grp.3s = c("Larynx", "Kaposi sarcoma"),
    cancer.grp.4m = c(),
    cancer.grp.4s = c("Lip, oral cavity", "Ovary", "Oesophagus", "Gallbladder", "Nasopharynx", "Other pharynx"),
    cancer.grp.5m = c(),
    cancer.grp.5s = c(),
    cancer.sex =
      c(
        "All cancers excl. non-melanoma skin cancer" = "Both sexes",
        "Prostate" = "Male", 
        "Breast" = "Female", 
        "Colorectum" = "Both sexes",
        "Melanoma of skin" = "Both sexes", 
        "Thyroid" = "Both sexes", 
        "Kidney" = "Both sexes",
        "Bladder" = "Both sexes", 
        "Corpus uteri" = "Female", 
        "Non-Hodgkin lymphoma" = "Both sexes", 
        "Testis" = "Male", 
        "Leukaemia" = "Both sexes", 
        "Multiple myeloma" = "Both sexes", 
        "Hodgkin lymphoma" = "Both sexes",
        "Lung" = "Both sexes", 
        "Pancreas" = "Both sexes", 
        "Brain, nervous system" = "Both sexes",
        "Cervix uteri" = "Female", 
        "Liver" = "Both sexes", 
        "Stomach" = "Both sexes",
        "Larynx" = "Both sexes", 
        "Kaposi sarcoma" = "Both sexes",
        "Lip, oral cavity" = "Both sexes", 
        "Ovary" = "Female", 
        "Oesophagus" = "Both sexes", 
        "Gallbladder" = "Both sexes", 
        "Nasopharynx" = "Both sexes", 
        "Other pharynx" = "Both sexes"
      ),
    gg.dat = format_data_for_plot(dat.asr, dat.edi)
  )

## 3. rank.cancer.plan: order cancer types according to their burden ----
rank.cancer.plan <-
  drake_plan(
    asr.tot = 
      compute_asr_tot(
        gg.dat %>% 
        mutate(cancer.sex.sel = plyr::mapvalues(cancer, names(cancer.sex), cancer.sex)) %>%
        filter(cancer %in% names(cancer.sex), cancer.sex.sel == sex)
      ),
    cancer.group.1.label.ranked = rank_cancer(c(cancer.grp.1m, cancer.grp.1s), asr.tot),
    cancer.group.2.label.ranked = rank_cancer(c(cancer.grp.2m, cancer.grp.2s), asr.tot),
    cancer.group.3.label.ranked = rank_cancer(c(cancer.grp.3m, cancer.grp.3s), asr.tot),
    cancer.group.4.label.ranked = rank_cancer(c(cancer.grp.4m, cancer.grp.4s), asr.tot),
    cancer.group.1.hsv = colorRampPalette(RColorBrewer::brewer.pal(4, "Reds")) (length(cancer.group.1.label.ranked)),
    cancer.group.2.hsv = colorRampPalette(RColorBrewer::brewer.pal(4, "Purples")) (length(cancer.group.2.label.ranked)),
    cancer.group.3.hsv = colorRampPalette(RColorBrewer::brewer.pal(4, "Greens")) (length(cancer.group.3.label.ranked)),
    cancer.group.4.hsv = colorRampPalette(RColorBrewer::brewer.pal(4, "Blues")) (length(cancer.group.4.label.ranked))
  )



## 4. gam.fit.plan produce models of cancer incidence/mortality ASR vs EdI per cancer type -----
##    and sex 
gam.fit.plan <-
  drake_plan(
    gam.fit.tab = 
      gam_fit_tab(
        gg.dat,
        group.extrem.index = .05
      )
  )

## 5. tab.prod.plan: produce tabS1 and tabS2 -----
tab.prod.plan <-
  drake_plan(
    ## tab1 will show by country the EdI vakue and classification and associated incidence 
    ## and mortality ASR for both sex and all cancers 
    tab1 = 
      gg.dat %>% 
      filter(sex == 'Both sexes', cancer == 'All cancers excl. non-melanoma skin cancer') %>%
      dplyr::select(country, type, sex, cancer, EdI, ASR) %>%
      group_by(sex, cancer) %>%
      spread(type, ASR) %>% 
      ungroup() %>%
      mutate(
        EdI_cat = cut(EdI, breaks = c(edi.breaks, 1), labels = names(edi.breaks)),
        `cancer incidence` = round(`cancer incidence`, digits = 1), 
        `cancer mortality` = round(`cancer mortality`, digits = 1)
      ) %>% 
      dplyr::select(
        `EdI category` = EdI_cat,
        Country = country,
        EdI = EdI,
        `Incidence rate all cancers` = `cancer incidence`,
        `Mortality rate all cancers` = `cancer mortality`
      ) %>%
      arrange(EdI),
    tab1.write = write.csv(
      tab1, 
      file = file_out(!!file.path('../results', out.version, 'SupplementTable1.csv'))
    ),
    ## tab2 is the mean value of mortality and incidence ASR by cancer and EdI class
    tab2 = 
      gg.dat %>% 
      mutate(cancer.sex.sel = plyr::mapvalues(cancer, names(cancer.sex), cancer.sex)) %>%
      filter(cancer %in% names(cancer.sex), cancer.sex.sel == sex) %>%
      dplyr::select(country, type, sex, cancer, EdI, ASR) %>%
      mutate(
        EdI_cat = cut(EdI, breaks = c(edi.breaks, 1), labels = names(edi.breaks))
      ) %>%
      group_by(sex, cancer, EdI_cat, type) %>%
      summarise(ASR_mean = mean(ASR) %>% round(digits = 1)) %>%
      ungroup() %>%
      mutate(ASR_label = paste('EdI', EdI_cat, type)) %>%
      dplyr::select(cancer, sex, ASR_label, ASR_mean) %>%
      group_by(cancer, sex) %>%
      spread(ASR_label, ASR_mean) %>% ungroup %>%
      dplyr::select(
        Cancer = cancer, Sex = sex, `EdI Low cancer incidence`, `EdI Low cancer mortality`,
        `EdI Medium cancer incidence`, `EdI Medium cancer mortality`, `EdI High cancer incidence`,
        `EdI High cancer mortality`, `EdI Very High cancer incidence`, `EdI Very High cancer mortality`
      ) %>%
      ungroup() %>%
      arrange(Cancer, Sex),
    tab2.write = 
      write.csv(
        tab2, 
        file = file_out(!!file.path('../results', out.version, 'SupplementTable2.csv'))
      )
  )

## 6. scatter.plot.plan: produce fig1, fig2, figS1 and figS2 -----
scatter.plot.plan <-
  drake_plan(
    scatter.outliers = read_csv(file_in('../data/outliers.table.csv')),
    scatter.plots = 
      plot_scatter(
        gam.fit.tab,
        edi.breaks = edi.breaks,
        scatter.outliers = scatter.outliers
      ),
    fig1 = 
      arrange_scatter(
        scatter.plots,
        cancer.grp.1m, cancer.grp.1s,
        cancer.grp.2m, cancer.grp.2s,
        cancer.grp.3m, cancer.grp.3s,
        cancer.grp.4m, cancer.grp.4s,
        cancer.grp.5m, cancer.grp.5s,
        type = 'fig1',
        countries.to.highlight = countries.to.highlight,
        cancer.sex = cancer.sex
      ),
    fig1.save = 
      ggsave(
        plot = fig1, 
        filename = file_out(!!file.path('../results', out.version, 'fig-1.pdf')),
        width = 210, height = 297 / 2, units = "mm"
      ),
    fig2.4x3 = 
      arrange_scatter(
        scatter.plots,
        cancer.grp.1m, cancer.grp.1s,
        cancer.grp.2m, cancer.grp.2s,
        cancer.grp.3m, cancer.grp.3s,
        cancer.grp.4m, cancer.grp.4s,
        cancer.grp.5m, cancer.grp.5s,
        type = 'fig2-4x3',
        cancer.sex = cancer.sex
      ),
    fig2.4x3.save = 
      ggsave(
        plot = fig2.4x3, 
        filename = file_out(!!file.path('../results', out.version, 'fig-2.pdf')),
        width = 297, height = 210, units = "mm", scale = 1
      ),
    fig.s1.4x3 = 
      arrange_scatter(
        scatter.plots,
        cancer.grp.1m, cancer.grp.1s,
        cancer.grp.2m, cancer.grp.2s,
        cancer.grp.3m, cancer.grp.3s,
        cancer.grp.4m, cancer.grp.4s,
        cancer.grp.5m, cancer.grp.5s,
        type = 'fig-s1-4x3',
        cancer.sex = cancer.sex
      ),
    fig.s1.4x3.save = 
      ggsave(
        plot = fig.s1.4x3, 
        filename = file_out(!!file.path('../results', out.version, 'fig-s1.pdf')),
        width = 297, height = 210, units = "mm"
      ),
    fig.s2.4x3 = 
      arrange_scatter(
        scatter.plots,
        cancer.grp.1m, cancer.grp.1s,
        cancer.grp.2m, cancer.grp.2s,
        cancer.grp.3m, cancer.grp.3s,
        cancer.grp.4m, cancer.grp.4s,
        cancer.grp.5m, cancer.grp.5s,
        type = 'fig-s2-4x3',
        cancer.sex = cancer.sex
      ),
    fig.s2.4x3.save = 
      ggsave(
        plot = fig.s2.4x3, 
        filename = file_out(!!file.path('../results', out.version, 'fig-s2.pdf')),
        width = 297, height = 210, units = "mm"
      )
  )

## 7. gam.stack.plot.plan: produce fig3 -----
gam.stack.plot.plan <-
  drake_plan(
    fig3 = 
      gam_stack_plot(
        gam.fit.tab %>% filter(sex %in% c("Male", "Female")),
        edi.breaks,
        cancer.group.1.label.ranked,
        cancer.group.2.label.ranked,
        cancer.group.3.label.ranked,
        cancer.group.4.label.ranked,
        cancer.group.1.hsv, 
        cancer.group.2.hsv, 
        cancer.group.3.hsv, 
        cancer.group.4.hsv,
        split.group = FALSE,
        add.others = TRUE,
        add.all.line = FALSE
      ) + 
      facet_wrap(
        ~ sex + type,
        dir = 'v',
        ncol = 2,
        labeller = label_wrap_gen(multi_line = FALSE, width = 200)
      ),
    fig3.save =
      ggsave(
        fig3,
        filename = file_out(!!file.path('../results/', out.version, 'fig-3.svg')),
        width = 297, height = 210, units = 'mm'
      )
  )

## 8.1 combined plan ----
plan.full <- 
  bind_plans(
    read.dat.plan,
    reshape.dat.plan,
    rank.cancer.plan,
    gam.fit.plan,
    tab.prod.plan,
    scatter.plot.plan,
    gam.stack.plot.plan
  )

## 8.2 configure plan ----
config.full <- drake_config(plan.full)

## visualisation of code strucure
save.vis <- FALSE ## if TRUE save the graph on the hard drive else display on the graph windows
vis_drake_graph(
  config.full, 
  file = if(save.vis) '../docs/ediglobocan-graph.html' else character(0), 
  selfcontained = if(save.vis) TRUE else FALSE
)

## 8.3 run plan ----
# clean(destroy = TRUE)
make(plan.full)

