library(tidyverse)
library(apc)
library(ggplot2)
library(reshape2)

function(input, output, session) {
  output$display <- renderPlot({
    disease <- input$DV_cancer
  })
  
  output$fitting <- renderPlot({
    disease <- input$Fitting_cancer
    type <- switch(input$Fitting_comp,
                   "Age Effects" = "age",
                   "Period Effects" = "period",
                   "Cohort Effects" = "cohort")
    # paste0("you have select ", disease, " and ", type, ".")
    directory <- paste0("data/",disease,"/")
    myList <- readRDS(paste0(directory,disease,".plotlist.rds"))
    lens <- length(myList)/3
    age.idx <- 1:lens
    per.idx <- (lens+1):(2*lens)
    coh.idx <- (2*lens+1):length(myList)
    
    idx <- function(type){
      switch(type,
             "age" = age.idx,
             "period" = per.idx,
             "cohort" = coh.idx
      )
    }
    
    glist <- lapply(myList[idx(type)], ggplotGrob)
    marrangeGrob(glist, ncol = 2, nrow=lens/2,top="")
  }, width=1000, height=600)
  
  observeEvent(input$update,{
    disease <- input$Forecast_cancer
    directory <- paste0("data/",disease,"/")
    files <- dir(directory)[grep(paste0(disease,".RData"), dir(directory))]
    # env <- environment()
    lapply(paste0(directory,files), load, .GlobalEnv)
    updateSliderInput(session, "cutoff", min=min.age, max=max.age, value=(min.age+max.age)/2)
    updateSliderInput(session, "Forecast_per", max=input$cutoff-min.age, value=(input$cutoff-min.age)/2)
  })
  
  output$forecast <- renderPlot({
    disease <- input$Forecast_cancer
    per <- input$Forecast_per
    cutoff.age <- input$cutoff
    
    ymax <- switch(disease,
                   "Lung Cancer"=400,
                   "Pancreatic Cancer"=40,
                   "Prostate Cancer"=140,
                   "Liver Cancer"=140
                   )
    if(is.null(ymax)) ymax <- 80
    
    from <- 1998
    to <- 2016
    pred.from <- to + 1
    pred.to <- to + per
    
    directory <- paste0("data/",disease,"/")
    files <- dir(directory)[grep(paste0(disease,".RData"), dir(directory))]
    # env <- environment()
    lapply(paste0(directory,files), load, .GlobalEnv)
    # rm(myList)     # this is for fitting
    # source("boot_apc.R")
    
    agegroup <- matrix(c(min.age, cutoff.age, 
                         cutoff.age+1,max.age, 
                         min.age, max.age), nrow=3, byrow=TRUE)
    
    
    # re-fit by different age group
    
    death.names <- ls(.GlobalEnv)[grep("mat.death",ls(.GlobalEnv))]
    pop.names <- ls(.GlobalEnv)[grep("mat.pop",ls(.GlobalEnv))]
    
    # needed output
    # case.sum, case.pred
    # pop, pop.pred
    
    # pop -> smooth.spline() -> pop.pred
    # data -> apc.list() -> apc.sums() -> sums
    # data + pop -> apc.list() -> apc.fit.model() -> fit.apc -> apc.forecast() -> pred
    
    
    scale <- 1e5
    
    pred.pop <- NULL
    sums <- NULL
    pred <- NULL
    myList <- list()
    List <- list()
    for (i in seq_along(death.names)) {
      apc.data.list(
        response=get(death.names[i]),
        dose=get(pop.names[i]),
        data.format="AP",
        age1=min.age,
        per1=from,
        unit=1,
        coh1=NULL,
        per.zero=NULL,
        per.max=NULL,
        time.adjust=0
      ) %>%
        apc.fit.model("poisson.response", "AC") -> fit.apc
      
      for (j in 1:nrow(agegroup)) {
        # 1. forecast pop data
        pred.pop <- c(pred.pop, paste0(pop.names[i], ".", agegroup[j,1],".",agegroup[j,2]))
        assign(pred.pop[(i-1)*nrow(agegroup)+j],
               get(pop.names[i])[agegroup[j,1]:agegroup[j,2]-min.age+1,] %>% 
                 colSums() %>% 
                 smooth.spline(x=from:to,spar=0.53) %>%
                 predict(pred.from:pred.to) %>% .$y)
        
        # 2. get apc.sums
        sums <- c(sums, paste0(substr(death.names[i], 11, 14), ".", agegroup[j,1],".",agegroup[j,2], ".sums"))
        assign(sums[(i-1)*nrow(agegroup)+j],
               apc.data.list(
                 response=get(death.names[i])[agegroup[j,1]:agegroup[j,2]-min.age+1,],
                 dose=NULL,
                 data.format="AP",
                 age1=agegroup[j,1],
                 per1=from,
                 unit=1,
                 coh1=NULL,
                 per.zero=NULL,
                 per.max=NULL,
                 time.adjust=0
               ) %>%
                 apc.data.sums(data.type="r")
        )
        
        # 3. get the pred
        pred <- c(pred, paste0(substr(death.names[i], 11, 14), ".", agegroup[j,1],".",agegroup[j,2], ".pred"))
        assign(pred[(i-1)*nrow(agegroup)+j],
               fit.apc %>%
                 apc.forecast.ac(
                   sum.per.by.age=c(ifelse(agegroup[j,1]==min.age,2,agegroup[j,1]-min.age+1), 
                                    agegroup[j,2]-min.age+1)
                 ) %>%
                 .$response.forecast.per.by.age %>%
                 `*`(scale) %>%
                 as.data.frame() %>%
                 filter(between(row_number(), pred.from-pred.from+1, pred.to-pred.from+1)) %>%
                 `/`(get(paste0(pop.names[i],".",agegroup[j,1],".", agegroup[j,2]))) %>%
                 cbind(year=pred.from:pred.to) %>%
                 mutate(hi=forecast+se, lo=forecast-se, variable=paste0("V",j)) %>%
                 select(-c(se.proc, se.est, se)) %>%
                 rename(value=forecast)
        )
      }
      myList[[i]] <- 
        apply(agegroup,1,function(x){
          get(paste0(substr(death.names[i],11,14),".",x[1],".",x[2],".sums")) %>%
            .$sums.per %>%
            `/`(get(pop.names[i])[x[1]:x[2]-min.age+1,] %>% colSums) %>%
            `*`(scale)
        }) %>%
        cbind(year=from:to) %>%
        as.data.frame() %>%
        melt(id.vars="year") %>%
        mutate(hi=NA, lo=NA) %>% 
        rbind(
          get(pred[(i-1)*nrow(agegroup)+1]), get(pred[(i-1)*nrow(agegroup)+2]), get(pred[(i-1)*nrow(agegroup)+3])
        )
      
      myList[[i]] %>%
        filter(year %in% from:to) %>%
        ggplot(aes(x=year,y=value, colour=variable)) +
        geom_point() +
        scale_x_continuous(limits=c(from, pred.to)) +
        scale_y_continuous(limits=c(0,ymax)) +
        geom_smooth(method="loess", span=0.7, linetype=2, se=FALSE) +
        geom_line(data=myList[[i]] %>% filter(year %in% pred.from:pred.to),
                  aes(x=year, y=value, colour=variable), size=1) +
        geom_ribbon(data=myList[[i]] %>% filter(year %in% pred.from:pred.to),
                    aes(ymin=lo, ymax=hi, fill=variable), alpha=0.2) +
        labs(x="",y="annual deaths per 100,000 population", 
             title=paste0(toupper(substr(death.names[i],13,14)),"-born ",
                          ifelse(substr(death.names[i],11,11)=="F", "Female", "Male")),
             colour="", fill="") +
        theme_bw() +
        theme(legend.position = c(0.2,0.9),
              legend.background = element_rect(fill = alpha("white", 0.0)))+
        scale_colour_discrete(
          labels=c(
            paste0("under ",cutoff.age),
            paste0("over ",cutoff.age),
            "All population")
        ) + 
        scale_fill_discrete(
          labels=c(
            paste0("under ",cutoff.age),
            paste0("over ",cutoff.age),
            "All population")
        ) -> List[[i]]
    }
    rm(list=ls(.GlobalEnv)[ls(.GlobalEnv)!="List"])
    
    glist <- lapply(List, ggplotGrob)
    
    # pdf("test.pdf", width=8, height=6)
    marrangeGrob(glist, ncol = 2, nrow=length(pop.names)/2,top="")
    # dev.off()
  }, width=1000, height=600)
}