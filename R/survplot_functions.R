#Jenny Smith

# March 30, 2017
#Purpose: Create a list of survival analyses and Given survival package Survfit() object, plot the survival curves and calculate the P values for differences in survival.



#' Create Kaplan-Meier plots
#'
#' @param fit the resulting object from surival::survfit
#' @param LegendTitle a string
#' @param timeUnit a string (years, months, days, etc)
#' @param colors a character vector of the same length as the number of groups in the KM plot
#' @param pval an option string wit the pvalue
#' @param max.year a numeric value to set x-axis limits.
#'
#' @return
#' @export
#'
#' @examples
#' fit <-  survival::survfit(survival::Surv(time, status) ~ x, data = survival::aml)
#' SurvivalPlot(fit=fit, LegendTitle="AML", timeUnit="months", colors=c("red","blue"))
SurvivalPlot <- function(fit, LegendTitle, timeUnit, colors,pval=NULL,max.year=NULL){
  #fit the output of the survFit() function.
  #legend title is a character vector which will be the plot legend title.
  #timeUnit is a character vector for the follow-up time in days, years, months etc.
  #colors is a character vector the same length as the number of groups.
  #max.year is the time (singe numeric) in years to end the plot at.

  if(is.null(pval)){
    pval <- ""
  }else{
    pval <- paste0("p = ",pval)
  }

  if(!is.null(max.year)){
    # if(max(fit$time) < max.year){
    pos.x <- max.year
    # }
  }else{
    pos.x <- max(fit$time)
  }

  if(is.null(fit$strata)){
    lm <- 1 #left margin (lm)
    rm <- 4  #right margin(rm)
    num.cols <- 1
    num.rows <- 1

    #main survival plot
    p <- GGally::ggsurv(fit, surv.col = colors,
                cens.col=colors, CI=FALSE,
                lty.est = 1, size.est = 1.25,
                cens.size = 2.0,order.legend=FALSE)

  }else{

    m <- max(nchar(gsub("^.+\\=(.+)", "\\1", names(fit$strata)))) #maximum number of characters
    d <- length(fit$strata) #how many groups
    lm <- dplyr::case_when(m <= 2 ~ 0.90,
                    m <= 5 ~ 1.25,
                    m > 5 & m <= 8 ~ 1.75,
                    m > 8 & m < 13 ~ 2.75,
                    m >= 13 & m < 17 ~ 3.4,
                    m >= 17 ~ 4.0)#left margin (lm)
    rm <- dplyr::case_when(d <= 4 ~ 1.5, d > 4 ~ 2.0) #right margin(rm)
    num.cols <- 2
    num.rows <- 2

    if(d > 4 & d < 7){
      num.rows <- 3
    }else if(d > 7 & d <= 18){
      num.cols <- 3
      num.rows <- 6
    }else if(d > 18){
      num.cols <- 5
      num.rows <- 5
    }

    if(num.cols*num.rows < d){
      num.rows <- 10
    }

    #main survival plot
    p <- GGally::ggsurv(fit, surv.col = colors, CI=FALSE,
                lty.est = 1, size.est = 1.25,
                cens.size = 2.0,order.legend=FALSE)
  }

  #customize the plot
  p <- p +
    ggplot2::labs(y= "Fraction Surviving", x = paste(timeUnit, "Follow-up")) +
    ggplot2::scale_y_continuous(limits = c(0,1),
                       breaks = seq(0,1,0.2),
                       minor_breaks = NULL) +
    ggplot2::scale_x_continuous(limits = c(0,pos.x),
                       breaks = seq(0,pos.x, 2)) +

    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=18),
          panel.background = ggplot2::element_rect(fill="white"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour = "black", fill= NA, size=1.0),

          axis.text = ggplot2::element_text(colour = "black"),
          axis.text.y = ggplot2::element_text(size = 18),
          axis.text.x = ggplot2::element_text(size=18),
          axis.title = ggplot2::element_text(size = 18),
          legend.position = "top",
          legend.direction  = "horizontal",
          legend.margin = ggplot2::margin(0,0,0,0, unit = "pt"),
          legend.text = ggplot2::element_text(size=12),
          legend.title = ggplot2::element_text(size=10),
          legend.key=ggplot2::element_blank()) +

    ggplot2::theme(plot.margin = ggplot2::margin(0, rm,.5,lm,unit="cm")) +
    ggplot2::annotate(geom="text", x=1, y=0.05, label=pval, size=5) +
    ggplot2::guides(color=ggplot2::guide_legend(ncol=num.cols, nrow=num.rows))

  if (length(colors) > 1){
    p$scales$scales[[1]]$name <- LegendTitle
    p$scales$scales[[2]]$name <- LegendTitle
  }

  return(p)
}


#risk table dataframe
#' Title
#'
#' @param survFit.obj the resulting object from  survival::survfit()
#' @param f.levels factor levels of the groups
#' @param times a numeric vector of the summary times for the KM estimates
#' @param col a string indicating the text color
#' @param ret.data boolean - return the summary table or not
#' @param max.year a numeric of the the x-axis limit
#'
#' @return
#' @export
#'
#' @examples
#' fit <-  survival::survfit(survival::Surv(time, status) ~ x, data = survival::aml)
#' risk.table(survFit.obj=fit)
risk.table <- function(survFit.obj, f.levels=NULL, times=NULL,col="black",
                       ret.data=FALSE, max.year=NULL){
  #survFit.obj is the results  survfit()
  #f.levels is an optional character vector to relevel the order of the groups, if desired. else levels is alphabetical
  #times would be a seperate numeric vector, with years/x-axis breaks points, if desired.

  if(!is.null(max.year)){
    # if (max(survFit.obj$time) < max.year){
      times <- seq(0,max.year, by= 2)
      xlims <- max.year
    # }
  }else{
    times <- seq(0,max(survFit.obj$time), by= 2)
    xlims <- max(survFit.obj$time)
  }

  #Summary dataframe with the selected time points to get % OS/ EFS for plot x-axis
  summ <- summary(survFit.obj, times=times, extend=TRUE)


  if(is.null(f.levels)){
    f.levels <- unique(as.character(summ$strata)) %>%
      gsub("^.+\\=(.+)", "\\1", .)
  }

  #Create a new dataframe to hold the time,KM estimate and strata for plotting
  risk.df <- dplyr::bind_cols(list(time=summ$time,
                            surv=summ$surv,
                            n.risk=summ$n.risk,
                            Group=as.character(summ$strata))) %>%
    magrittr::set_colnames(c("time","surv","N.Risk","Group")) %>%
    dplyr::mutate(Group=gsub("^.+\\=(.+)", "\\1", Group),
           Group=factor(Group, levels = f.levels))


  if(ret.data){
    return(risk.df)
  }

  #variable for deciding length of margins. longer group names would need more space.
  m <- max(nchar(levels(risk.df$Group))) # print(c("length:", m))
  d <- length(levels(risk.df$Group))
  lm <- dplyr::case_when(m <= 5 ~ 1.75,
                  m > 5 & m <= 8 ~ 1.35,
                  m > 8 & m <= 10 ~ 1.25,
                  m > 10 & m < 13 ~ 1.10,
                  m >= 13 ~ 1.20)

  rm <- dplyr::case_when(d <= 5 ~ 1.5,
                  m > 5 ~ 2.0)


  #gglot the risk table (risk.df)
  tbl <- ggplot2::ggplot(risk.df, ggplot2::aes(x = time, y = Group, label=N.Risk)) +
    ggplot2::geom_text(size = 5)  +
    ggplot2::theme_bw()  +
    ggplot2::scale_x_continuous("Number at risk",
                       limits = c(0,xlims),
                       breaks = times) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          legend.position="none",
          axis.line = ggplot2::element_line(size=0.7, color="black"),
          axis.title.x = ggplot2::element_text(size=18),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(face="bold",
                                     size=14,
                                     color= col[levels(risk.df$Group)],
                                     hjust=1)) +
    #of note: this margin is still not perfect..... the size of window in Rstudio is different than what actually shows up on the tiff...
    ggplot2::theme(plot.margin = ggplot2::margin(-1.5, rm, 0.1, lm, unit="cm"))
  return(tbl)
}


#' Create a dataframe output to summarize the results of Cox proportional hazards regression
#'
#' @param coxph.res the results of survival::coxph()
#' @param Colname a character string of the column name
#'
#' @return
#' @export
#'
#' @examples
#' coxreg <- survival::coxph(survival::Surv(time, status) ~ x, data = survival::aml)
#' coxSummaryTable(coxreg)
coxSummaryTable <- function(coxph.res,Colname="Group"){
  c.mod <- coxph.res
  c.summ <- summary(coxph.res)
  c.table <- as.data.frame(c.summ$coefficients) %>% cbind(.,as.data.frame(c.summ$conf.int))

  c.HR <- round(c.table$`exp(coef)`, digits=3)
  c.CI <- paste(round(c.table$`lower .95`,digits=3),
                round(c.table$`upper .95`, digits = 3), sep="-")
  c.pVal <- round(c.table$`Pr(>|z|)`, digits=3)

  comp <- gsub("^.+(Yes|No|Unknown|High|Low)", "\\1",names(c.mod$coefficients))
  c.stats <- data.frame(names(c.mod$coefficients) ,comp, c.HR, c.CI, c.pVal)
  colnames(c.stats) <- c("Variable", Colname,"Hazard Ratio", "95% Confidence Interval", "P-value")

  return(c.stats)
}



#' Calculate the p-value for differnces in KM curves from survdiff() function
#'
#' @param survdiff the results from survival::survdiff()
#' @param digits a numeric value of the number of digits for the p-value
#'
#' @return
#' @export
#'
#' @examples
#' diff <- survival::survdiff(survival::Surv(time, status) ~ x, data = survival::aml)
#' calc_KMcurve_pvalues(diff)
calc_KMcurve_pvalues <- function(survdiff, digits=3){
  p <-  pchisq(survdiff$chisq, length(survdiff$n)-1,
              lower.tail = FALSE) %>%
    round(.,  digits = digits)

  return(p)
}



#' Create a table with OS/EFS percent at specified time point
#'
#' @param fit results of survival::survfit()
#' @param time the time point at which to summarize OS % and EFS %
#' @param pvalues a string of the p-value to include in the table
#'
#' @return
#' @export
#'
#' @examples
#' fit <-  survival::survfit(survival::Surv(time, status) ~ x, data = survival::aml)
#' outcome_table(fit)
outcome_table <- function(fit, time=5, pvalues=NULL){

  summ <- summary(fit, times=time)
  type <- ifelse(any(grepl("OS", summ$call)),"OS","EFS")
  summ <- summ[!grepl("table|std.chaz|call", names(summ))]
  summ <- summ[sapply(summ,length)==length(summ$strata)]

  surv_col <- paste("Percent", type)
  n <- length(summ$strata)

  tab <- tibble::tibble(as.data.frame(summ)) %>%
    dplyr::mutate_at(vars(surv,lower,upper,std.err), ~round(.*100, digits = 2)) %>%
    dplyr::mutate(Group=stringr::str_split_fixed(strata, pattern = "=", n=2)[,2])

  tab <- tab %>%
    dplyr::select(Group,
                  "Number Patients per Group"=n,
                  "Time (years)"=time,
                  !!surv_col := surv,
                  "Lower 95% CI"=lower,
                  "Upper 95% CI"=upper,
                  "Standard Error"=std.err)


  if(!is.null(pvalues) & length(pvalues) < n){
    pval_ref <- NA #ref is the first position in the strata
    pvalues=c(pval_ref, pvalues) #must be in same order as strata in the survift()!

    tab <- tab %>%
      tibble::add_column("p-value"=pvalues, .after=7)
  }

  return(tab)
}

