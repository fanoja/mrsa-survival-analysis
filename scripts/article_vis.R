source("scripts/utilities.R")
source("scripts/network_vis.R")

load_library("ggplot2")
load_library("gridExtra")
load_library("cowplot")
load_library("bayesplot")
load_library("reshape2")
load_library("rstanarm")
load_library("igraph")

theme_set(theme_cowplot(font_family = "Times", font_size=12)) # Font size and type

## Exploratory Visualization ##

get_binom_ci <- function(n_subset, n_all, cl = 0.95){
    #' Returns the min and max confidence intervals based on n_subset and n_all.
    #'
    #' If confidence interval becomes negative, make it zero instead. Use with lapply().
    #'
    #' @param n_subset N subset.
    #' @param n_all N all values.
    #' @param cl Confidence level, by default 0.95.

    # Correct for negative CI:
    if (n_subset != 0 & n_all != 0){
        bintest <- binom.test(n_subset, n_all, conf.level = cl)
        ci <- as.numeric(bintest$conf.int)
        minci <- ci[1]
        maxci <- ci[2]
    }
    else {
        minci <- 0
        maxci <- 0
    }

    return(c(minci, maxci))
}

get_resclear <- function(data, treatment){
    #' Calculate the number of resistant, non-resistant, cleared+resistant and cleared+non-resistant observations from data for a given treatment.
    #'
    #' @param data Survival data frame.
    #' @param treatment Covariate of interest as a string.
    #' @return resmat A matrix of the following counts: resistant & cleared, non-resistant & cleared, resistant, non-resistant.
    
    # Caclulate the number of resistant/non-resistant cases and the number of non-R/R & C cases
    n_res <- dim(data[which(data[,treatment] == 1),])[1]
    n_nonres <- dim(data[which(data[,treatment] == 0),])[1]
    n_res_clear <- dim(data[which(data[,treatment] == 1 & data$y != 0),])[1]
    n_nonres_clear <- dim(data[which(data[,treatment] == 0 & data$y != 0),])[1]
    
    # Collect these to a matrix:
    resmat <- matrix(c(n_res_clear, n_nonres_clear, n_res, n_nonres), ncol = 4, nrow = 1)
    colnames(resmat) <- c("rc", "nrc", "res", "non_res")
    
    return(resmat)

}

get_resclear_barplot <- function(data, treatment, get_numeric = FALSE, additional_title = ""){
    #' Create a barplot of the resistance and clearance probability of a treatment from survival data.
    #'
    #' @param data Survival data frame.
    #' @param treatment One of the antibiotics of interest, for example "Mupirocin".
    #' @return p_bar A ggplot2 bar plot.
    
    resmatD <- get_resclear(data[which(data$ARM == 1),], treatment)
    resmatE <- get_resclear(data[which(data$ARM == 0),], treatment)
    
    if (get_numeric){
        print("Decolonization:")
        print(resmatD)
        print("Education:")
        print(resmatE)
    }
    
    minci_res <- c(get_binom_ci(resmatD[,"rc"], resmatD[,"res"])[1], get_binom_ci(resmatE[,"rc"], resmatE[,"res"])[1])
    maxci_res <- c(get_binom_ci(resmatD[,"rc"], resmatD[,"res"])[2], get_binom_ci(resmatE[,"rc"], resmatE[,"res"])[2])
    minci_nonres <- c(get_binom_ci(resmatD[,"nrc"], resmatD[,"non_res"])[1], get_binom_ci(resmatE[,"nrc"], resmatE[,"non_res"])[1])
    maxci_nonres <- c(get_binom_ci(resmatD[,"nrc"], resmatD[,"non_res"])[2], get_binom_ci(resmatE[,"nrc"], resmatE[,"non_res"])[2])

    plot_result <- data.frame("p_rc" = c(resmatD[,"rc"]/resmatD[,"res"], resmatE[,"rc"]/resmatE[,"res"]),
                                "p_nrc" =c(resmatD[,"nrc"]/resmatD[,"non_res"], resmatE[,"nrc"]/resmatE[,"non_res"]),
                                "ARM" = c("D", "E")
                             )

    df <- melt(plot_result)
    df$minci <- c(minci_res, minci_nonres)
    df$maxci <- c(maxci_res, maxci_nonres)

    p_bar <- ggplot(df, aes(x=variable, y=value, fill = ARM, group = ARM)) +
                geom_bar(stat = "identity", position = position_dodge())+
                ylim(0,1)+
                # Use the three lines below, if you wish to have error bars
                #geom_point(position = position_dodge(0.9), size = 5) +
                #geom_line(position = position_dodge(0.9)) +
                geom_errorbar( aes(x=variable, ymin=minci, ymax=maxci), width=.5, size=1.0, position=position_dodge(0.9)) + 
                labs( x = "", y = expression(paste("Probability of clearance at ", v[1])),
                     title = paste(additional_title, treatment), fill = "ARM") +
                scale_x_discrete(labels = c("p_rc" = expression(paste("Resistant at ", v[0])), "p_nrc"=expression(paste("Non-resistant at ", v[0]))),limits=c("p_rc", "p_nrc")) +
                scale_fill_manual(values=c(color_scheme_get("blue")$light_highlight, color_scheme_get("pink")$light_highlight)) 
                #theme_minimal()
    #print(p_bar)
    
    return(p_bar)
    
    
}


## Bayesian Inference Visualizations ###

# Quantile calculations for xlim:

get_param_quantiles <- function(fit, q=0.95){
    #' Get quantiles from fit.
    #'
    #' @param fit Fitted model.
    #' @param q By default 0.95 quantiles.
    #' @return qs Quantiles of interest.
    
    minq <- (1-q)/2
    maxq <- 1-(1-q)/2
    qs <- apply(fit, 2, quantile, probs = c(minq, maxq))
    
    return(qs)
}

get_mods_quantile_xlims <- function(list_of_mods, params, q=0.95){
    #' Get quantile-based x-axis limits for a list of models that are to be plotted
    #'   in the same figure. For example, a list of mods by site etc.
    #'
    #' Returns a list with min and max x-axis limit values for plotting with the mcmc_intervals -function. This function is especially useful when you want to plot multiple mcmc_intervals figures with the same x-axis limits - this automatically finds the largest x-axis range from the fitted models.
    #'
    #' @param list_of_mods A list of the models of interest.
    #' @param params A vector of covariates of interest, such as "Mupirocin" and "Clindamycin".
    
    minx <- c()
    maxx <- c()
    
    for (mod in list_of_mods){
        pars <- names(mod$coefficients)[names(mod$coefficients) %in% params]
        fit <- as.matrix(mod, pars = pars)
        qs <- get_param_quantiles(fit, q)
        minx <- c(minx, min(qs))
        maxx <- c(maxx, max(qs))

    }

    return(list("minx" = min(minx), "maxx" = max(maxx)))
}

custom_mcmc_intervals <- function(mod, pars, prob_outer = 0.95, prob_inner = 0.5){
    #' Customized mcmc_intervals function, where the zero line is below the intervals.
    #' 
    #' @param mod Model of interest
    #' @param params Parameters to plot
    #' @param prob_outer Probability of the outer interval
    
    ## Create df with credible intervals & medians
    mod <- as.matrix(mod)
    
    ci95 <- posterior_interval(mod, prob = prob_outer)
    ci50 <- posterior_interval(mod, prob = prob_inner)
    
    if (length(pars) > 1){
        medians <- apply(mod[,pars], 2, median)
        }
    else{
        medians <- median(mod[,pars])
        names(medians) <- pars}

    df <- cbind(ci95, ci50)
    colnames(df) <- c("cil", "ciu", "ciml", "cimu") # lower, upper, middle lower, middle upper credible intervals
    
    if (length(pars) == 1){
        df <- cbind(t(df[pars,]), medians)
    }else{
        df <- cbind(df[pars,], medians)
    }
    
    df <- cbind(df, data.frame("loc" = seq(dim(df)[1], 1, -1)))

    df <- df[order(nrow(df):1),]
    

    ## Plot the figure
    p <- ggplot(df, aes(medians, rownames(df))) +
    geom_vline(xintercept = 0, color = "darkgray") +
    geom_errorbar(aes(xmin = cil, xmax = ciu),
                  color = color_scheme_get()$mid,
                 width = 0, size = 0.5) +
    geom_errorbar(aes(xmin = ciml, xmax = cimu), #aes(xmin = ciml, xmax = cimu)
              color = color_scheme_get()$dark,
             size = 2, width = 0) +
    geom_point(shape = 21, size = 4,
               colour = color_scheme_get()$dark_highlight,
               fill = color_scheme_get()$light) +
    xlab("") +
    ylab("") +
    scale_y_discrete(limits = rownames(df))
    
    p
}

site_figs <- function(list_of_mods, plotfun, params, as_list = FALSE, is_model1 = FALSE, no_xlim = FALSE, ...){
    #' Call plotfun on each site-specific model.
    #'
    #' @param list_of_mods List of the fitted models of interest.
    #' @param plotfun Plotting function to be used, such as the mcmc_intervals function from the bayesplot package.
    #' @param params Covariates of interest in a vector, such as c("Mupirocin", "Chlorhexidine").
    #' @param as_list If true, returns a list of the ggplot2 figure objects instead of returning the figures as a 2x2 grid object (see package cowplot, plot_grid).
    #' @param is_model1 If true, the models passed to this function where separate models for each treatments. 
    #' @param no_xlim If true, do not constrain x-axis.

    # Assumptions: site mods are named mod_n, mod_t, mod_s and mod_w
    # The title will always be just the corresponding site
    # Passes extra arguments to the plotting function.

    if (is_model1){
        modn <- extract_params(modn)
        modt <- extract_params(modt)
        mods <- extract_params(mods)
        modw <- extract_params(modw)
        
        list_of_mods <- list("mod_n" = modn, "mod_t" = modt, "mod_s" = mods, "mod_w" = modw)
    }
    
    # Find a common x-axis: 
    
    xlims <- get_mods_quantile_xlims(list_of_mods, params = params)
    maxx <- xlims$maxx
    minx <- xlims$minx
    
    if(no_xlim){
        pn <- plotfun(list_of_mods$mod_n, params = params, ...) + ggtitle("Nares")
        pt <- plotfun(list_of_mods$mod_t, params = params, ...) + ggtitle("Throat")
        ps <- plotfun(list_of_mods$mod_s, params = params, ...) + ggtitle("Skin")
        pw <- plotfun(list_of_mods$mod_w, params = params, ...) + ggtitle("Wound")  
        }
    else{
        pn <- plotfun(list_of_mods$mod_n, params = params, ...) + ggtitle("Nares") + coord_cartesian(xlim = c(minx, maxx))
        pt <- plotfun(list_of_mods$mod_t, params = params, ...) + ggtitle("Throat") + coord_cartesian(xlim = c(minx, maxx))
        ps <- plotfun(list_of_mods$mod_s, params = params, ...) + ggtitle("Skin") + coord_cartesian(xlim = c(minx, maxx))
        pw <- plotfun(list_of_mods$mod_w, params = params, ...) + ggtitle("Wound") + coord_cartesian(xlim = c(minx, maxx))
    }
    # Return a list of figures
    
    if (as_list){
        return(list("pn" = pn, "pt" = pt, "ps" = ps, "pw" = pw))
    }
    
    return(plot_grid(pn, pt, ps, pw))
}

arm_figs <- function(modd, mode, plotfun, params, no_xlim = FALSE, is_model1 = FALSE, as_list = FALSE, title = "Impact of resistance on MRSA clearance", vjust = 1, ...){
    #' Plot a figure for fitted decolonization and education models as specified by plotfun.
    #'
    #' @param modd Fitted model for the decolonization arm.
    #' @param mode Fitted model for the education arm.
    #' @param plotfun Plotting function, such as the mcmc_intervals function from the bayesplot package.
    #' @param params Parameters to be plotted, such as the beta parameter for mupirocin ("Mupirocin")
    #'
    #' For the rest of the function parameters, see site_figs.
    
    if (is_model1){
        modd_ext <- extract_params(modd, params = params)
        mode_ext <- extract_params(mode, params = params)
        xlims <- get_mods_quantile_xlims(list("modd" = modd_ext, "mode" = mode_ext), params = params)
    }else{
         xlims <- get_mods_quantile_xlims(list("modd" = modd, "mode" = mode), params = params)
    }
    
    if (!no_xlim){
        color_scheme_set("blue")
        pd <- plotfun(modd, params = params, col_scheme = "blue", is_model1 = is_model1, ...) + ggtitle("Decolonization") + coord_cartesian(xlim = c(xlims$minx, xlims$maxx))
        color_scheme_set("pink") 
        pe <- plotfun(mode,params = params,col_scheme = "pink", is_model1 = is_model1,...) + ggtitle("Education") + coord_cartesian(xlim = c(xlims$minx, xlims$maxx))
        }
    else{
        color_scheme_set("blue")
        pd <- plotfun(modd, params = params, col_scheme = "blue", is_model1 = is_model1,...) + ggtitle("Decolonization")
        color_scheme_set("pink")
        pe <- plotfun(mode,params = params, col_scheme = "pink", is_model1 = is_model1,...) + ggtitle("Education")
        }
    
    if (as_list){
        return(list("pd" = pd, "pe"= pe))
     }
    
     title <- get_grid_title(title, vjust = vjust)

    return(plot_grid(title, pd, pe, ncol = 1, rel_heights = c(0.1, 1, 1)))
    
    #return(plot_grid(pd, pe))
}

# Extract parameters from model1 (each antibiotic in a separate model), used by get_mcmc_intervals

extract_params <- function(mod1, params = treatments){
    #' Add the treatments from separate models in model1 to a single data frame.
    #'
    #' Enables plotting each antibiotic in the same figure using mcmc_intervals, instead of plotting one antibiotic/antimicrobial per figure.
    #'
    #' @param mod1 Separately fitted models for each antimicrobial.
    #' @param params Antimicrobials of interest.
    #' @return df All estimated parameters of interest from separate models in one data frame.
    
    df <- array()
    for (m in mod1){
        coef_name <- names(m$coefficients)[names(m$coefficients) %in% params]
        df <- cbind(df, as.matrix(m, pars = coef_name))
    }
    df <- as.data.frame(df[,-1])
    
    return(df)
    
}

# MCMC intervals, not site-specific



get_mcmc_intervals <- function(mod, params, title = "Posterior Intervals", is_model1 = FALSE, col_scheme = "blue"){
    #' Plot a figure using a custom version of the mcmc_intervals function from the bayesplot figure.
    #'
    #' @param mod Fitted model of interest.
    #' @param title Title for the figure.
    #' @param is_model1 True indicates separate models for all antibiotics.
    #' @param col_scheme Color scheme used in the figure, see bayesplot mcmc_intervals for details.
    #' @return pd Ggplot2 object.
    
    if (is_model1){
        
        mod <- extract_params(mod, params = params)
        pars <- params
    }
    else{
        pars <- names(mod$coefficients)[names(mod$coefficients) %in% params]
    }
    
    color_scheme_set(col_scheme)

    pd <- custom_mcmc_intervals(mod, pars = pars, prob_outer = 0.95) +
    #vline_0(color = "darkgray", size = 0.75) +
    #geom_vline(xintercept = 0, color = "darkgray", size = 0.75) +
    labs(title = title)
    
    return(pd)
}

# Summary statistics: posterior draw means, quantiles and SD

get_summary_statistics <- function(mod, params = treatments, re_param = "(Intercept)", raneff = "host_id"){
    #' Get posterior mean, SD, median and 95% credible interval.
    #'
    #' Code modified from: https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html
    #'
    #' @param mod Fitted model of interest.
    #' @param params Antimicrobials of interest in a character vector.
    #' @param re_param Random effect param - here, we use Intercept.
    #' @param raneff Name of the random effect of interest. Either host_id or strain.
    #' @return summary_stats Return the aforementioned statistics in a data frame.

    # Extract random effects (host level errors)
    if (re_param == "(Intercept)"){
        re_param = "\\(Intercept\\)"
    }
    
    draws <- as.matrix(mod, regex_pars = paste(paste("b\\[", re_param, sep = ""), paste(paste(" ", raneff, sep = ""), "\\:", sep = ""), sep = ""))
    
    draw_means <- apply(draws, 2, mean)
    draw_sd <- apply(draws, 2, sd)

    # Posterior 95% credible interval
    draw_q <- data.frame(t(apply(draws, 2, quantile, probs = c(0.025, 0.50, 0.975))))
    names(draw_q) <- c("Q2.5", "Q50", "Q97.5")

    # Combine summary statistics of posterior simulation draws
    summary_stats <- data.frame(draw_means, draw_sd, draw_q)
    
    return(summary_stats)    
}

get_posterior_summary_arms <- function(modD, modE, mod_id1 = "D", mod_id2 = "E", raneff){
    #' Call get_summary_statistics for fitted models modD (decolonization) and modE (education) and combine these in a data frame.
    
    postsumD <- get_summary_statistics(modD, raneff = raneff)
    postsumE <- get_summary_statistics(modE, raneff = raneff)

    postsum <- rbind(postsumD, postsumE)
    postsum$arm <- c(rep(mod_id1, dim(postsumD)[1]), rep(mod_id2, dim(postsumE)[1]))
    postsum$raneff <- c(rep(raneff, dim(postsumD)[1]), rep(raneff,  dim(postsumE)[1]))
    
    return(postsum)
}


## Utilities

get_grid_title <- function(title, angle = 0, vjust = 1.5){
    #' Return a plot_grid (see package cowplot) compatible title object.
    
    hjust = 0
    if(angle != 0){
        hjust = 0.5
    }
    
      title <- ggdraw() + 
      draw_label(
        title,
        fontface = 'bold',
        fontfamily = "Times",
        x = 0,
        hjust = hjust,
          vjust = vjust,
        angle = angle
      )
    
    return(title)
    
}  
                         
save_figs_by_site <- function(sitefigs, filename, figdir = fig_savedir){
    #' Save a list of site-specific figures as pdf.
    #'
    #' @param sitefigs A list of site-specific figure objects. Assumes the following naming convention for sitefigs list: nares = pn, throat = pt, skin = ps, wound = pw.
    #' @param filename Name of the pdf file. It will be combined with a site name (nares, skin, throat or wound).
    #' @param figdir The directory where the pdf file will be saved.
    
    #print(paste(figdir, paste("nares", filename, sep = "_"), sep = ""))
          
    pdf(paste(figdir, paste("nares", filename, sep = "_"), sep = ""))
    print(sitefigs$pn)
    dev.off()

    pdf(paste(figdir, paste("throat", filename, sep = "_"), sep = ""))
    print(sitefigs$pt)
    dev.off()

    pdf(paste(figdir, paste("skin", filename, sep = "_"), sep = ""))
    print(sitefigs$ps)
    dev.off()

    pdf(paste(figdir, paste("wound", filename, sep = "_"), sep = ""))
    print(sitefigs$pw)
    dev.off()
    
}
