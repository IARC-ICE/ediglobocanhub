## load the rquired libraries
library(readr)
library(dplyr)
library(tidyr)
library(rlang)
library(lubridate)
library(svglite)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

#' Combine Cancer Incidence and Mortality ASR data with EdI 
#'
#' @param dat.asr data.frame, cancer incidence and mortality ASR per country, cancer type and sex
#' @param dat.edi data.frame, EdI data per country 
#'
#' @return data.frame, the combined dataset use in the rest of the analysis
format_data_for_plot <-
  function(
    dat.asr, 
    dat.edi
  ){
    inner_join(
      dat.asr,
      dat.edi,
      by = 'iso3'
    ) %>%
      mutate(sex = factor(sex, levels = c('Male', 'Female', 'Both sexes')))
  }

#' Fit gam model (ASR vs EdI) and store the projections in a suitable data.frame 
#'
#' This is the key function of the analysis. 
#' Here is fit a quasi-poisson cubic spline model on the log of the rates (incidence and 
#' mortality ASR) of the countries according to their EdI.
#' A model per cancer type, type of rate (incidence/mortality) and sex (both sexes/male/female) 
#' is computed. 
#' Models evalution metrics and models projection over a gradient of EdI is also retuned.
#' No weight is used because we are interested in the overall trends accross EdI gradient (
#' all countries have the same weights).
#' Extrem EdI values are 'recentered' to prevent for border effect (group.extrem.index param).
#' Only a subselection of all computed models will be displayed in the main paper.
#' 
#' @param gg.dat data.frame, the output of format_data_for_plot function
#' @param group.extrem.index the percentage of EdI values considered as extrem and recentered to
#'        the closest EdI quantile value.
#'
#' @return nested data.frame, contains models evaluation metrics and projection over a gradient  
#'         of EdI values.
gam_fit_tab <-
  function(
    gg.dat,
    group.extrem.index = .05
  ){
    require(mgcv)
    ## construct the gam ang get the rsquared associated
    gamFun <-
      function(
        dat, 
        group.extrem.index
      ){
        df_ <- 
          dat %>%
          mutate(
            EdI =
              EdI %>%
              replace(
                EdI < quantile(EdI, probs = group.extrem.index / 2), 
                quantile(EdI, probs = group.extrem.index / 2)
              ) %>%
              replace(
                EdI > quantile(EdI, probs = 1 - (group.extrem.index / 2)), 
                quantile(EdI, probs = 1 - (group.extrem.index / 2))
              )
            
          )
        
        ## note: after testing we choose a quasipoisson family model with a log link
        ## no weights have assigned to countries since we want to model general trends 
        ## accross EdI gradient
        b <- gam(ASR ~ s(EdI, bs = 'cr'), data = df_, family = quasipoisson) 
        # AER::dispersiontest(glm(ASR ~ splines::bs(EdI), data = df_, family = poisson), trafo = 1)
        s <- summary(b)
        df_pred_ <- 
          tibble(
            EdI = 
              seq(min(df_$EdI, na.rm = TRUE), 
                  max(df_$EdI, na.rm = TRUE), 
                  length.out = 300)
            )
        ## make model prediction along the EdI gradient
        b_pred_ <- predict(b, df_pred_, type = 'response', se.fit = TRUE)
        ## store all the results in a list
        l <- 
          tibble(
            EdI = df_pred_$EdI, 
            ASR = b_pred_$fit,
            ASR.se = b_pred_$se.fit,
            mod.name = "gam",
            mod.disp = s$dispersion,
            mod.rsq = s$r.sq,
            mod.devexpl = s$dev.expl,
            mod.edf = s$edf,
            mod.aic = AIC(b),
            hdi.effect.pval = s$s.table['s(EdI)', 'p-value']
          )
        return(l)
      }
    
    ## apply the function to the per cancer, sex and type nested dataset
    gg.dat %>% 
      group_by(cancer, sex, type) %>% 
      nest() %>%
      mutate(
        gam_pred_tab = 
          purrr::map(
            data,
            ~ gamFun(.x, group.extrem.index)
          )
      )
  }

#' Produce single incidence and mortality ASR vs EdI scatter plots panel
#' 
#' individual panel for Fig1, Fig2, FigS1 and FigS2 of the manuscript are produced 
#' through this function.
#'
#' @param gam.fit.tab nested data.frame, output of gam_fit_tab function
#' @param edi.breaks numeric, EdI break points
#' @param scatter.outliers data.frame, the oulier ASR to be displayed as text on plots
#' @param text.size integer, the default text size
#' @param log.scale logical, should log scale transformation of the axis be performed?
#'
#' @return nested data.frame, scatter plot and associated legend is produced for each
#'         cancer type and sex
#' 
plot_scatter  <- 
  function(
    gam.fit.tab,
    edi.breaks = NULL,
    scatter.outliers = NULL,
    text.size = 10,
    log.scale = FALSE
  ){

    gg <- 
      gam.fit.tab %>%
      ungroup() %>%
      # filter(
      #   !(cancer %in% c('All cancers excl. non-melanoma skin cancer'))
      # ) %>%
      ## to prevent merging warning messages
      mutate(
        cancer = as.character(cancer),
        sex = as.character(sex)
      ) %>%
      group_by(cancer, sex) %>%
      nest() %>%
      mutate(
        plot = 
        purrr::pmap(
          list(.data = data, .cancer = cancer, .sex = sex),
          function(.data, .cancer, .sex) {
            
            ggplot(
              data =
                data.frame(
                  .data %>% select(type, data) %>% unnest(cols = c(data)),
                  cancer = .cancer,
                  sex = .sex,
                  stringsAsFactors = FALSE
                ) %>%
                anti_join(
                  scatter.outliers,
                  by = c('cancer', 'sex', 'type', 'country')
                ),
              aes(
                x = EdI,
                y = ASR,
                colour = type
              )
            ) +
            geom_point(
              aes(
                size = ptotal
              ),
              alpha = .5,
              shape = 21
            ) +
            xlab("") + #xlab(gg_hdi_index)  +
            ylab("") + #ylab("Age-Standardized Rate (per 100 000)") +
            scale_size_continuous(
              name = "population size (millions):",
              breaks = 10^(7:9),
              labels = function(x) scales::comma(x, big.mark = ',', scale = 10^(-6))
            ) +
            scale_colour_manual(name = "", values = c("#4682b4", "#dc143c")) +
            scale_x_continuous(limits = c(0.25, 1), breaks = seq(.35, .95, .15)) +
            scale_y_continuous(
              trans = if(log.scale) 'log1p' else 'identity',
              breaks = if(.cancer %in% c('Brain, nervous system', 'Larynx', 'Nasopharynx')){ if(log.scale) log1p(c(0, 5, 10)) else c(0, 5, 10)} else waiver()
            ) + 
            geom_vline(xintercept = edi.breaks[2:4], lty = 2, colour = "lightgrey") +
            geom_text(
              data =
                data.frame(
                  .data %>% select(type, data) %>% unnest(cols = c(data)),
                  cancer = .cancer,
                  sex = .sex,
                  stringsAsFactors = FALSE
                ) %>%
                inner_join(
                  scatter.outliers,
                  by = c('cancer', 'sex', 'type', 'country')
                ),
              aes(y = ASR_fake, label = label, hjust = hjust),
              size = 2.5,
              show.legend = FALSE
            ) +
            geom_line(
              data =
                data.frame(
                  .data %>% select(type, gam_pred_tab) %>% unnest(cols = c(gam_pred_tab)),
                  cancer = .cancer,
                  sex = .sex
                ),
              lty = 1,
              show.legend = FALSE
            ) +
            facet_wrap(
              ~ cancer + sex,
              ncol = 2,
              labeller = label_wrap_gen(multi_line=FALSE, width = 200),
              scale = "free_y"
            ) +
            theme(
              text = element_text(size = text.size),
              strip.background = element_rect(fill = NA, color = NA, linetype = 1),
              strip.text = element_text(hjust = 0.02),
              panel.background = element_rect(fill = NA),
              legend.position = "bottom",#none",
              legend.text = element_text(size = rel(1.2)),
              legend.title = element_text(size = rel(1.2)),
              plot.margin = unit(c(0.3, 0.1, 0.3, 0), "cm"),
              axis.text = element_text(size = rel(.8)),
              axis.title.x = element_blank(),
              axis.line.x.bottom = element_line(colour = 'black'),
              axis.line.y.left = element_line(colour = 'black')
            )
          }
        ),
        legend = 
          purrr::map2(plot, cancer,
            ~ if(.y == 'All cancers excl. non-melanoma skin cancer') {
              cowplot::get_legend(
                (.x + 
                  theme(
                    legend.title = element_text(size = 12),
                    legend.text = element_text(size = 12),
                    legend.key = element_rect(fill = NA)
                  )
                )
              )
              } else {
                cowplot::get_legend(
                  (.x + theme(legend.key = element_rect(fill = NA)))
                )
              }),
        plot = purrr::map(plot, ~ .x + theme(legend.position = "none"))
      )
  }


#' Combine incidence and mortality ASR vs EdI scatter plots
#' 
#' Fig1, Fig2, FigS1 and FigS2 of the manuscript are produced through this function.
#'
#' @param scatter.plots nested data.frame, output of plot_scatter function
#' @param cancer.grp.1m character, group A cancer type to be plotted in main graphs
#' @param cancer.grp.1s character, group A cancer type to be plotted in supp graphs
#' @param cancer.grp.2m character, group B cancer type to be plotted in main graphs
#' @param cancer.grp.2s character, group B cancer type to be plotted in supp graphs
#' @param cancer.grp.3m character, group C cancer type to be plotted in main graphs
#' @param cancer.grp.3s character, group C cancer type to be plotted in supp graphs
#' @param cancer.grp.4m character, group D cancer type to be plotted in main graphs
#' @param cancer.grp.4s character, group D cancer type to be plotted in supp graphs
#' @param cancer.grp.5m character, extra group cancer type to be plotted in main graphs
#' @param cancer.grp.5s character, extra group cancer type to be plotted in supp graphs
#' @param type character, the type of graph to be produced ('fig1', 'fig2-4x3', 'fig-s1-4x3'
#'        'fig-s2-4x3')
#' @param countries.to.highlight data.frame, the countries to be highlited in 'fig1'
#' @param text.size integer, default font size
#' @param cancer.sex character, the sex to be considered for each cancer type
#'
#' @return ggplot object, 'fig1', 'fig2', 'fig-s1' or 'fig-s2' is returned
arrange_scatter <-
  function(
    scatter.plots,
    cancer.grp.1m, cancer.grp.1s,
    cancer.grp.2m, cancer.grp.2s,
    cancer.grp.3m, cancer.grp.3s,
    cancer.grp.4m, cancer.grp.4s,
    cancer.grp.5m, cancer.grp.5s,
    type = 'fig2-4x3',
    countries.to.highlight = NULL,
    text.size = 14,
    cancer.sex = NULL
  ){
    library(cowplot)
    
    ## conbine individual scater plots in a strip
    merge_plot <- 
      function(cancer.grp, nrow = 1, ncol = NULL){
        cowplot::plot_grid(
          plotlist = 
            scatter.plots %>% 
            ungroup() %>%
            mutate(cancer.sex.sel = plyr::mapvalues(cancer, names(cancer.sex), cancer.sex)) %>%
            filter(
              cancer %in% cancer.grp,
              cancer.sex.sel == sex
            ) %>%
            mutate(cancer = forcats::fct_relevel(cancer, cancer.grp)) %>%
            arrange(cancer) %>%
            pull(plot),
          nrow = nrow,
          ncol = ncol
        ) + 
          theme(
            plot.background = element_blank() # element_rect(color = 'lightgrey', linetype = 1)
          )
      }
    
    ## produce independently x and y axis legend
    axis_label_2 <-
      function(label, angle, size){
        grid::textGrob(label, just = "center", rot = angle, gp = grid::gpar(fontsize = size))
      }
    y.lab <- 
      axis_label_2(
        label = "Age-Standardized Rate (per 100,000)", 
        angle = 90, size = text.size
      )
    x.lab <- 
      axis_label_2(
        label = "Education and Income index (EdI)", 
        angle = 0, size = text.size
      )
    legend.glob <- scatter.plots$legend[[2]]

    ## produce different plots according to type arg
    if(type == 'fig2-4x3')
    {
      p1m.r <- merge_plot(cancer.grp.1m, nrow = 2, ncol = 3)
      p2m.r <- merge_plot(cancer.grp.2m, nrow = 1, ncol = 3)
      p3m.r <- merge_plot(cancer.grp.3m, nrow = 1, ncol = 3)
      
      fig2l.1 <-
        cowplot::plot_grid(
          p1m.r, 
          p2m.r, 
          p3m.r,
          labels = c('Group A', 'Group B', 'Group C'),
          label_y = 1.01, label_x = - 0.025,
          label_size = .8 *  text.size,
          rel_heights = c(2, 1, 1),
          nrow = 3
        ) 
      
      fig2l.2 <-
        plot_grid(
          plot_grid(y.lab, fig2l.1, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, x.lab, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, legend.glob, nrow = 1, rel_widths = c(.05, .95)),
          ncol = 1, rel_heights = c(.9, .05, 0.05), scale = .98
        )
      
      fig.out <- fig2l.2
    } 
    else if(type == 'fig-s1-4x3') 
    {
      p1s.r <- merge_plot(cancer.grp.1s, nrow = 3, ncol = 3)
      p3s.r <- merge_plot(cancer.grp.3s, nrow = 1, ncol = 3)
      figs2l.1 <-
        cowplot::plot_grid(
          p1s.r, 
          p3s.r,
          nrow = 2,
          labels = c('Group A', 'Group C'),
          label_y = 1.01, label_x = - 0.025,
          label_size = .8 *  text.size,
          rel_heights = c(3, 1)
        ) 
      
      figs2l.2 <-
        plot_grid(
          plot_grid(y.lab, figs2l.1, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, x.lab, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, legend.glob, nrow = 1, rel_widths = c(.05, .95)),
          ncol = 1, rel_heights = c(.9, .05, 0.05), scale = .98
        )
      
      fig.out <- figs2l.2
      
    } else if(type == 'fig-s2-4x3') 
    {
      p4s.r <- merge_plot(cancer.grp.4s, nrow = 2, ncol = 3)
      
      figs3l.1 <-
        cowplot::plot_grid(
          p4s.r,
          NULL,
          nrow = 2,
          labels = c('Group D', ''),
          label_y = 1.01, label_x = - 0.025,
          label_size = .8 *  text.size,
          rel_heights = c(2, 2)
        ) 
      
      figs3l.2 <-
        plot_grid(
          plot_grid(y.lab, figs3l.1, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, x.lab, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, legend.glob, nrow = 1, rel_widths = c(.05, .95)),
          ncol = 1, rel_heights = c(.9, .05, 0.05), scale = .98
        )
      
      fig.out <- figs3l.2
      
    } else if(type == 'fig1') 
    {
      fig1.dat <- 
        scatter.plots %>% 
        filter(
          cancer == 'All cancers excl. non-melanoma skin cancer',
          sex == 'Both sexes'
        ) %>%
        select(cancer, sex, data) %>%
        unnest(cols = c(data)) %>%
        select(cancer, sex, type, data) %>%
        unnest(cols = c(data)) %>%
        select(cancer, sex, country, type, EdI, ASR) %>%
        inner_join(
          countries.to.highlight %>%
            select(cancer, sex, country, country_lab, hjust, nudge_x, vjust, nudge_y),
          by = c('cancer', 'sex', 'country')
        )
      
      fig1 <-
        scatter.plots %>% 
        filter(
          cancer == 'All cancers excl. non-melanoma skin cancer'
        ) %>%
        select(cancer, sex, plot) %>%
        pull(plot) %>% .[[1]] +
        geom_label(
          aes(
            label = country_lab
          ),
          data = fig1.dat,
          size = 2,
          show.legend = FALSE
        )
      
      legend.glob <- scatter.plots$legend[[1]]
        
      
      figs1.m <-
        plot_grid(
          plot_grid(y.lab, fig1, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, x.lab, nrow = 1, rel_widths = c(.05, .95)),
          plot_grid(NULL, legend.glob, nrow = 1, rel_widths = c(.05, .95)),
          ncol = 1, rel_heights = c(.9, .05, 0.05)
        )
      
      fig.out <- figs1.m
    }  else stop('unknown type')
    
    return(fig.out)
  }


#' Rank cancer based on sum of ASR over countries
#'
#' @param cancer.group character, a list of cancer of interest
#' @param asr.tot data.frame, the output of compute_asr_tot function
#'
#' @return character, the list of cancer order according to their overall burden (sum of incidence
#'         and mortality ASR)
#'
rank_cancer <- 
  function(
    cancer.group, 
    asr.tot
  ){
    
    cancer.group.label.ranked <-
      asr.tot %>%
      filter(
        cancer %in% cancer.group
      ) %>%
      arrange((ASR.rank)) %>%
      pull(cancer)
    
    cancer.group.label.ranked
  }


#' Compute the sum of all ASR or a set of countries
#' 
#' This will help us to order cancer according to their overall burden
#'
#' @param gg.dat data.frame, Incidence and mortality ASR per cancer type , per country 
#'        (opt. per sex)
#'
#' @return data.frame, a data.frame with sum of cancer incidence and mortality ASR and 
#'         a syntetic rank based on the sum of incidence and mortality ASR
#'
compute_asr_tot <-
  function(gg.dat){
    gg.dat %>%
    ungroup() %>%
    filter(
      cancer != 'All cancers excl. non-melanoma skin cancer'
      # sex == 'Both sexes'
    ) %>%
    group_by(cancer, type) %>%
    summarise(
      ASR.tot = sum(ASR, na.rm = TRUE)
    ) %>%
    spread(type, ASR.tot) %>%
    ungroup() %>%
    mutate(
      ASR.rank = rank((`cancer incidence` + `cancer mortality`))
    ) %>% 
    arrange(ASR.rank)
  }



#' Produce stack graph (e.g. fig3)
#' 
#' Stack all the incidence/mortality ASR relativce to EdI models produces in a previous stage.
#' Several option will help to find an optimal design.
#'
#' @param gam.fit.tab data.frame containing per country cancer type (opt. per sex) the incidence
#'                    and mortality ASR projected along EdI gradient.
#' @param edi.breaks  numeric, EdI break points
#' @param cancer.group.1.label.ranked character, the cancer in group A (ordered)
#' @param cancer.group.2.label.ranked character, the cancer in group B (ordered)
#' @param cancer.group.3.label.ranked character, the cancer in group C (ordered)
#' @param cancer.group.4.label.ranked character, the cancer in group D (ordered)
#' @param cancer.group.1.hsv character, the colors associated to cancer in group A
#' @param cancer.group.2.hsv character, the colors associated to cancer in group B
#' @param cancer.group.3.hsv character, the colors associated to cancer in group C
#' @param cancer.group.4.hsv character, the colors associated to cancer in group D
#' @param split.group logical, should cancer group be facetted
#' @param add.others logical, should other cancers (not in group A,B,C,D) be displayed
#' @param add.all.line logical, add all cancers dashed line
#'
#' @return a ggplot object, an ASR vs EdI projections stacked graph
gam_stack_plot <-
  function(
    gam.fit.tab,
    edi.breaks,
    cancer.group.1.label.ranked,
    cancer.group.2.label.ranked,
    cancer.group.3.label.ranked,
    cancer.group.4.label.ranked,
    cancer.group.1.hsv,
    cancer.group.2.hsv,
    cancer.group.3.hsv,
    cancer.group.4.hsv,
    split.group = TRUE,
    add.others = TRUE,
    add.all.line = TRUE
  ){
    
    gg.stack.dat <-
      ## get model predictions
      gam.fit.tab %>%
      select(cancer, sex, type, gam_pred_tab) %>%
      unnest(cols = c(gam_pred_tab)) %>%
      ungroup() %>%
      ## reorder cancers
      mutate(
        ASR = pmax(ASR, 0),
        cancer =
          factor(
            replace(
              cancer, 
              !cancer %in% 
                c(
                  cancer.group.1.label.ranked, cancer.group.2.label.ranked, cancer.group.3.label.ranked,
                  cancer.group.4.label.ranked, 'All cancers excl. non-melanoma skin cancer', 'All cancers', 
                  'Non-melanoma skin cancer'
                ),
              'Others/Unspec'
            ),
            levels = 
              c(
                cancer.group.1.label.ranked, cancer.group.2.label.ranked, cancer.group.3.label.ranked,
                cancer.group.4.label.ranked, 'All cancers excl. non-melanoma skin cancer', 'All cancers', 
                'Non-melanoma skin cancer', 'Others/Unspec'
              )
          ),
        cancer.grp1 =
          plyr::mapvalues(
            cancer,
            from = levels(cancer),
            to =
              c(
                rep('group A', length(cancer.group.1.label.ranked)),
                rep('group B', length(cancer.group.2.label.ranked)),
                rep('group C', length(cancer.group.3.label.ranked)),
                rep('group D', length(cancer.group.4.label.ranked)),
                rep('group 7', 4)
              )
          )
      ) %>%
      group_by(cancer, sex, type, EdI, cancer.grp1) %>%
      ## handle duplicates (if any)
      summarise(ASR = sum(ASR)) %>%
      ungroup()
    
    ## detect if we are working with Both sexes combined or with separated Male/Female 
    muti.sex.plot <- 
      gg.stack.dat %>%
      select(cancer, sex, type) %>%
      distinct() %>%
      group_by(cancer, type) %>%
      summarise(n = n()) %>%
      mutate(multi.sex = (n > 1)) %>%
      pull('multi.sex') %>%
      any()

    ## produce the plot
    gg.stack.001 <-
      gg.stack.dat %>%
      ## remove not relevant cancer types
      filter(
        cancer != 'All cancers excl. non-melanoma skin cancer',
        cancer != 'All cancers', 
        cancer != 'Non-melanoma skin cancer'
      ) %>%
      {
        if(!add.others) filter(., cancer != 'Others/Unspec') else .
      } %>%
      ggplot(
        aes(
          x = EdI,
          y = ASR,
          fill = factor(cancer)
        )
      ) +
      scale_fill_manual(
        values =
          c(cancer.group.1.hsv, cancer.group.2.hsv, cancer.group.3.hsv, cancer.group.4.hsv, if(add.others) 'grey90')
      ) +
      geom_area(position = 'stack', colour = 'lightgrey', size = .01) +
      facet_wrap(
        ~ type,
        ncol = 2,
        labeller = label_wrap_gen(multi_line=FALSE, width = 200)
      ) +
      labs(
        y = "Age-Standardized Rate (per 100 000)",
        x = 'Education and Income Index (EdI)',
        fill = 'cancer type' 
      ) +
      scale_y_continuous(labels = scales::comma) +
      theme(
        strip.background = element_rect(fill = NA, color = NA, linetype = 1),
        strip.text = element_text(hjust = 0.02),
        panel.background = element_rect(fill = NA),
        legend.position = "right",
        plot.margin = unit(c(0,0,0,0), "cm"),
        text = element_text(size = 14)#,
      ) +
      geom_vline(xintercept = edi.breaks[2:4], lty = 2, colour = "lightgrey") +
      guides(fill = guide_legend(ncol = 1))
    
    gg.stack.out <- 
      if(!split.group){
        gg.stack.001 + 
          if(add.all.line){
            geom_line(
              aes(x = EdI, y = ASR),
              colour = 'black',
              linetype = 2,
              data = 
                gg.stack.dat %>%
                filter(cancer == 'All cancers excl. non-melanoma skin cancer'),
              show.legend = FALSE,
              inherit.aes = FALSE
            )
          }
        
      } else{
        gg.stack.001 +
          facet_wrap(
            ~ cancer.grp1 + type,
            ncol = 2,
            labeller = label_wrap_gen(multi_line=FALSE, width = 200)
          ) 
      }
    
    gg.stack.out
  }

