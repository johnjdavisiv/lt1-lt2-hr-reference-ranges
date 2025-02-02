# Meta-analytic reference ranges for %HRmax, %HRR, and %VO2max at LT1 and LT2
# John J Davis
# RunningWritings.com

library(metafor)
library(tidyverse)
library(cowplot)


df <- read_csv("VO2 and heart rate at LT1 and LT2 meta-analysis.csv") %>%
  mutate(study = factor(study),
         study_cluster = factor(study_cluster)) %>%
  mutate(vi = (sd**2)/n_subjects) %>%
  mutate(pct_male = n_male/n_subjects) %>%
  # Use modified Siegel methods for logmean and log sd
  mutate(ybar = 100 - mean,
         ybarstar = log(
           ybar / sqrt(
             1 + (
               ((n_subjects - 1)/n_subjects)*(sd^2)
             ) / (ybar^2)
           )
           ),
         #Compare with siegel function, ensure they are the same
         ybarstar_siegel = log(ybar/sqrt(1 + ((n_subjects -1)/(n_subjects)) * sd^2/ybar^2)),
         si_star = sqrt(log(
           1 + (((n_subjects-1)/n_subjects)*(sd^2))/(ybar^2)
         )),
         si_star_siegel = sqrt(log(1 + ((n_subjects - 1)/(n_subjects)) *sd^2/ybar^2))
         ) %>%
  #Now try inverting transform - easy, you just do 100 minus exp(value) *but* need to do limits, densities, etc on log scale
  mutate(ybarstar_hi = ybarstar + 2*si_star,
         ybarstar_lo = ybarstar - 2*si_star,
         ybar_invert = 100 - exp(ybarstar),
         yhi_invert = 100 - exp(ybarstar_hi),
         ylo_invert = 100 - exp(ybarstar_lo)) %>%
  mutate(vi_star = (si_star**2)/n_subjects)

df %>% glimpse()


# Wow this actually looks so good, log(100-x) is so clutch
df %>% 
  filter(metric %in% c("LT1","LT2")) %>%
  filter(percent_metric %in% c("pct_vo2max","pct_hrmax")) %>%
  ggplot(aes(x=ybar_invert, y=study)) + 
  geom_linerange(aes(xmin = mean - 2*sd, xmax = mean  + 2*sd), color = "red", alpha = 0.3, linewidth = 4) + 
  geom_point(aes(x=mean), color = "red", size=5, alpha = 0.95)+
  geom_point(size=3) + 
  geom_linerange(aes(xmin=ylo_invert, xmax=yhi_invert), linewidth = 1.25) + 
  geom_vline(xintercept = 100, linewidth = 1, linetype = "dotted") + 
  facet_grid(percent_metric ~ metric) + 
  scale_x_continuous(limits = c(50,107),
                     breaks = seq(50,100,10))

#QQ mostly fine
df %>%
  filter(metric %in% c("LT1","LT2")) %>%
  filter(percent_metric %in% c("pct_vo2max","pct_hrmax")) %>%
  ggplot(aes(sample = ybarstar)) + 
  geom_qq_line() + 
  geom_qq(pch=21, fill="white", size=3) + 
  facet_grid(percent_metric ~ metric, scales = "free")


# --- also make plots decomposing stuff intov araince compoemnts 
# -- Setup functions: one for standard meta analysis, one for clusters

fit_log_meta_model <- function(input_df, training_metric, threshold_metric, conf_level = 0.9){
  #Input df is JUST DF, we filter inside the function
  # conf_level is confidence level of the reference range, not the mean estimate CIs!
  # training_metric - Heart rate, VO2max, etc OPTIONS: pct_hrmax, pct_vo2max, pct_hrreserve <-- care, one study for HRR
  use_df <- input_df %>% filter(percent_metric == training_metric,
                                metric == threshold_metric)
  
  print(sprintf("Analyzing %s as a function of %s in %d studies with a total of %d subjects",
        threshold_metric, training_metric, dim(use_df)[1], sum(use_df$n_subjects)))
  #Throw error if we have clusters within each study
  if (dim(use_df)[1] != length(unique(use_df$study_cluster))) stop("Clustered data detected! use multilevel function")
  
  #Fit the REML model
  meta_mod <- rma(yi = ybarstar, vi = vi_star, data = use_df, method = "REML")
  print(summary(meta_mod))
  
  #Mean of model (true mean)
  mu_hat <- meta_mod$beta[1]
  
  se_hat <- meta_mod$se
  
  #Tau squared, variation across studies
  tau_hat_sq <- meta_mod$tau2
  
  #Pooled variance - subject-level variance within studies (ref: Siegel 2021)
  sigma_hat_sq <- sum((use_df$n_subjects - 1) * use_df$si_star^2)/sum(use_df$n_subjects - 1)

  #Total variance - accounts for variation across studies AND variationa cross subjects within studies  
  total_var <- sigma_hat_sq + tau_hat_sq
  upper_tail <- (1+conf_level)/2
  lower_tail <- 1 - upper_tail
  
  #Reference range is ~ N(mean = mu_hat, variance = total_var), just use normal quantiles to get it 
  lower_ref <- qnorm(lower_tail, mean = mu_hat, sd = sqrt(total_var))
  upper_ref <- qnorm(upper_tail, mean = mu_hat, sd = sqrt(total_var))
  
  #Transform back
  orig_mu <- 100 - exp(mu_hat)
  orig_lower <- 100 - exp(upper_ref) #note flipped
  orig_upper <- 100 - exp(lower_ref)
  
  result_string <- sprintf("%s Estimate: %.1f %s, %.0f%% reference range: %.1f - %.1f",
          threshold_metric, 
          orig_mu,
          training_metric,
          conf_level*100,
          orig_lower,
          orig_upper)
  print(result_string)
  #Return basically everything
  result_list <- list(mean_estimate = mu_hat,
                      se_of_mean = se_hat,
                      tau_sq = tau_hat_sq,
                      sigma_sq = sigma_hat_sq,
                      total_var = total_var,
                      individual_fraction = sigma_hat_sq/total_var,
                      lab_fraction = tau_hat_sq/total_var,
                      conf_level = conf_level, 
                      lower_ref = lower_ref,
                      upper_ref = upper_ref,
                      orig_mean = orig_mu,
                      orig_upper = orig_upper,
                      orig_lower = orig_lower,
                      meta_model = meta_mod,
                      model_type = "single",
                      analysis_df = use_df,
                      result_string = result_string)
  return(result_list)
}



#Multilevel version - for studies with sub-group data (men, women, elite, recreational, etc.)
fit_log_ml_meta_model <- function(input_df, training_metric, threshold_metric, conf_level = 0.9){
  # Clamp range: truncate above 100?
  # Percent metric is the training metric, e.g. pct_hrmax or pct_vo2max or pct_hrreserve
  use_df <- input_df %>% filter(percent_metric == training_metric,
                                metric == threshold_metric)
  
  print(sprintf("MULTILEVEL: Analyzing %s as a function of %s in %d groups from %d studies with a total of %d subjects",
                threshold_metric, training_metric, 
                dim(use_df)[1],
                length(unique(use_df$study_cluster)),
                sum(use_df$n_subjects)))
  
  if (dim(use_df)[1] == length(unique(use_df$study_cluster))) stop("Should not have equal number of clusters and groups! use normal REML function")
  
  #Fit model, using mv / multilevel structure
  meta_mod <- rma.mv(yi = ybarstar, V = vi_star, random = ~ 1 | study_cluster/study, 
                   data = use_df, 
                   method = "REML")
  
  print(summary(meta_mod))
  
  meta_mod$sigma2 #2 element vector!
  
  #Mean of model (true mean)
  mu_hat <- meta_mod$beta[1]
  se_hat <- meta_mod$se
  
  #Tau squared, variation across studies
  tau_hat_1_sq <- meta_mod$sigma2[1]
  tau_hat_2_sq <- meta_mod$sigma2[2]
  
  #Pooled variance - subject-level variance within studies
  sigma_hat_sq <- sum((use_df$n_subjects - 1) * use_df$si_star^2)/sum(use_df$n_subjects - 1)
  
  # -- DANGER -- experimental unpublished statistics ahead! (assume uncorrelated level 1 & 2)
  #Total variance - accounts for variation across studies AND variation across subjects within studies  
  total_var <- sigma_hat_sq + tau_hat_1_sq + tau_hat_2_sq
  
  upper_tail <- (1+conf_level)/2
  lower_tail <- 1 - upper_tail
  
  #Transform back
  
  #Reference range is ~ N(mean = mu_hat, variance = total_var), just use normal quantiles to get it 
  lower_ref <- qnorm(lower_tail, mean = mu_hat, sd = sqrt(total_var))
  upper_ref <- qnorm(upper_tail, mean = mu_hat, sd = sqrt(total_var))
  
  #Transform back
  orig_mu <- 100 - exp(mu_hat)
  orig_lower <- 100 - exp(upper_ref) #note flipped bw/ 100-x transform
  orig_upper <- 100 - exp(lower_ref)

  
  result_string <- sprintf("%s Estimate: %.1f %s, %.0f%% reference range: %.1f - %.1f",
                           threshold_metric, 
                           orig_mu,
                           training_metric,
                           conf_level*100,
                           orig_lower,
                           orig_upper)
  print(result_string)
  
  result_list <- list(mean_estimate = mu_hat,
                      se_of_mean = se_hat,
                      tau_1_sq = tau_hat_1_sq,
                      tau_2_sq = tau_hat_2_sq,
                      sigma_sq = sigma_hat_sq,
                      total_var = total_var,
                      #DEBATABLE! Two reasonable ways to do it
                      individual_fraction = (sigma_hat_sq + tau_hat_2_sq)/total_var,
                      alt_individual_fraction = sigma_hat_sq/total_var,
                      #DEBATABLE, but since variance between groups in the same lab...should add!
                      lab_fraction = tau_hat_1_sq/total_var,
                      alt_lab_fraction = (tau_hat_1_sq + tau_hat_2_sq)/total_var,
                      conf_level = conf_level, 
                      lower_ref = lower_ref,
                      upper_ref = upper_ref,
                      orig_mean = orig_mu,
                      orig_lower = orig_lower,
                      orig_upper = orig_upper,
                      meta_model = meta_mod,
                      model_type = "cluster",
                      analysis_df = use_df,
                      result_string = result_string)
  
  return(result_list)
}

#Plotting / evaluation functions
make_forest <- function(meta_results){
  par(mfrow=c(1,1))
  forest(meta_results$meta_model, slab = meta_results$analysis_df$study)
}

#Need to check to make sure we actually have a maximum in the likelihood function (i.e. not overparameterized)
make_profile_plot <- function(meta_ml_results){
  #trying to be adaptive here
  if (any(class(meta_ml_results$meta_model) == "rma.uni")) {
    par(mfrow=c(1,1))
    profile(meta_ml_results$meta_model)
  } else {
    par(mfrow=c(2,1))
    profile(meta_ml_results$meta_model, sigma2=1)
    profile(meta_ml_results$meta_model, sigma2=2)
  }
}

# -- Generalized I**2 from URL: https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
get_generalized_I2 <- function(model) {
  # Construct the weight matrix
  W <- diag(1 / model$vi)
  
  # Model matrix (fixed-effects design matrix)
  X <- model.matrix(model)
  
  # Projection matrix
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  
  # Sum of all variance components
  total_heterogeneity <- sum(model$sigma2)
  
  # Effective sampling variance
  sampling_var <- (model$k - model$p) / sum(diag(P))
  
  # Denominator for generalized I^2 calculations
  total_variance <- total_heterogeneity + sampling_var
  
  # Overall I^2 (in %)
  I2_total <- 100 * total_heterogeneity / total_variance
  
  # Breakdown by each variance component (in %)
  I2_components <- 100 * model$sigma2 / total_variance
  
  # If sigma2 has no names, assign defaults (e.g., "sigma2.1", "sigma2.2", ...)
  if (is.null(names(I2_components))) {
    names(I2_components) <- paste0("sigma2.", seq_along(I2_components))
  }
  
  # Fraction due to sampling variance (in %)
  sampling_fraction <- 100 * sampling_var / total_variance
  
  # Return a list with results
  list(
    I2_total      = I2_total,
    I2_components = I2_components,
    Sampling      = sampling_fraction
  )
}

print_generalized_I2 <- function(x) {
  # x is the list returned by get_generalized_I2()
  
  # Overall I^2
  print(sprintf("Overall I^2: %.1f%%", x$I2_total))
  
  # Sampling variance fraction
  print(sprintf("Sampling variance: %.1f%%", x$Sampling))
  
  # Breakdown of variance components
  print("Breakdown by variance components:")
  for (comp_name in names(x$I2_components)) {
    print(sprintf("  %s: %.1f%%", comp_name, x$I2_components[[comp_name]]))
  }
}


make_plot_df <- function(fitted_results, ng = 501, grid_lower = 50){
  #Grid of values to evaluate dnorm across (don't go to 100 obviously)
  pct_metric_grid <- seq(grid_lower,99.9, length.out=ng)
  
  #Individual variance
  sd_individual <- sqrt(fitted_results$sigma_sq)
  
  if (fitted_results$model_type == "single"){
    sd_lab <- sqrt(fitted_results$tau_sq)
  } else {
    # For multilevel meta-analysis, groups within studies 
    #DEBATABLE: could also do tau_1_sq + tau_2_sq, depends on whether you think lab variance "leaks" into tau_2
    # In practice it ends up not being a huge difference in most cases
    sd_lab <- sqrt(fitted_results$tau_1_sq)
    #sd_lab <- sqrt(fitted_results$tau_1_sq + fitted_results$tau_2_sq)
  }
  
  log_x_grid <- log(100 - pct_metric_grid)
  
  # --- IMPORTANT! Do not forget the Jacobian correction ----
  #Individual variance
  dist_individual_star = dnorm(log_x_grid, mean = fitted_results$mean_estimate, sd = sd_individual)
  dist_individual_orig = dist_individual_star * (1 / (100 - pct_metric_grid)) #<-- Jacobian correction
  
  #Lab variance (and also different groups or  subject pools)
  dist_lab_star = dnorm(log_x_grid, mean = fitted_results$mean_estimate, sd = sd_lab)
  dist_lab_orig = dist_lab_star * (1 / (100 - pct_metric_grid)) #<-- Jacobian correction
  
  
  # Total variance (indivdiual + lab)
  dist_total_star = dnorm(log_x_grid, mean = fitted_results$mean_estimate, sd = sqrt(fitted_results$total_var))
  dist_total_orig = dist_total_star * (1 / (100 - pct_metric_grid)) #<-- Jacobian correction! 
  
  plot_df <- data.frame(pct_metric = pct_metric_grid,
                        log_pct_metric = log_x_grid,
                        dist_individual_star = dist_individual_star,
                        dist_individual_orig = dist_individual_orig,
                        dist_lab_star = dist_lab_star,
                        dist_lab_orig = dist_lab_orig,
                        dist_total_star = dist_total_star,
                        dist_total_orig = dist_total_orig)
  
  below_df <- plot_df %>%
    filter(pct_metric < fitted_results$orig_lower)
  
  above_df <- plot_df %>%
    filter(pct_metric > fitted_results$orig_upper)
  
  res_list <- list(plot_df = plot_df,
                   below_df = below_df,
                   above_df = above_df)
  return(res_list)
}



# ---- Plot setup, general ----
fnt <- 12
fnt_ax <- 12
fnt_title <- 12

#Axis ticks
tick_wd <- 0.08
arrow_lwd <- 0.4
ant_fnt <- 6
ax_lwd <- 0.25
ax_x_fnt <- 12

# Unify theme
theme_mee <-   theme(plot.title = element_text(hjust=0.5, size = fnt_title),
                     axis.line = element_line(color="black", linewidth=ax_lwd,
                                              lineend = "square"),
                     axis.ticks = element_line(color="black", linewidth=ax_lwd,
                                               lineend = "square"),
                     axis.ticks.length = unit(tick_wd, "cm"),
                     axis.text = element_text(color = "black", size = fnt),
                     axis.title = element_text(color = "black", size = fnt_ax),
                     plot.background = element_rect(fill="white", color="white"),
                     plot.caption = element_text(size = fnt, face = "bold"))


# --- Analyze! ---
df %>% glimpse()


#Plot on log(ybar*) range
df %>% 
  filter(metric %in% c("LT1","LT2")) %>%
  filter(percent_metric %in% c("pct_vo2max","pct_hrmax")) %>%
  ggplot(aes(x=ybarstar, y=study)) + 
  geom_point() + 
  geom_linerange(aes(xmin=ybarstar-2*si_star, xmax=ybarstar+2*si_star)) + 
  facet_grid(percent_metric ~ metric)


# ----- LT1, % HRmax  ------------ 

# analyze ybarstar, si_star as mean, sd
LT1_hr_results <- fit_log_meta_model(df, training_metric = "pct_hrmax", threshold_metric = "LT1", conf_level = 0.9)

make_profile_plot(LT1_hr_results)
make_forest(LT1_hr_results)


LT1_hr_plot_res <- make_plot_df(LT1_hr_results)


plot_df <- LT1_hr_plot_res$plot_df
below_df <- LT1_hr_plot_res$below_df
above_df <- LT1_hr_plot_res$above_df


plot_df %>% glimpse()

plot_df %>%
  ggplot(aes(x=pct_metric, y=dist_total_orig)) + 
  geom_line() + 
  geom_line(aes(y=dist_individual_orig), color = "blue") + 
  geom_line(aes(y = dist_lab_orig), color = "green")



sum(plot_df$dist_individual_orig)
sum(plot_df$dist_total_orig)


# --- Plot -----

ymax <- max(plot_df$dist_total_orig)*1.08
ymax_plot <- ymax*1.1

ant_col <- "#94a3b8"
ant_lwd <- 0.75
ant_type <- "solid"



within_alf <- 0.6
outer_alf <- 0.5

hr_fill <- "#4ade80"
LT1_col <- hr_fill

ribbon_fill <- hr_fill
ribbon_line_col <- "black"
ribbon_lwd <- 0.5

#Axis setup
x_limits <- c(50,100)
x_breaks <- seq(50,100,by=10)
x_labels <- sprintf("%d%%", x_breaks)
x_name <- "Percent of HRmax"

range_text <- sprintf("90%% range: %.0f - %.0f%% HRmax",
                      LT1_hr_results$orig_lower,
                      LT1_hr_results$orig_upper)

title_string <- "LT1 as a percentage of HRmax in runners"


subtitle_string <- sprintf("Pooled from %d studies on a total of %d runners",
                           length(LT1_hr_results$analysis_df$study_cluster),
                           sum(LT1_hr_results$analysis_df$n_subjects))

x_middle <- (LT1_hr_results$orig_lower + LT1_hr_results$orig_upper)/2


single_title_fnt <- 18
range_fnt <- 5

sub_fnt <- 13

#Now try with fill
total_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  #Lower and upper shade
  #Upper/lower
  annotate(geom="segment", x = LT1_hr_results$orig_lower, xend = LT1_hr_results$orig_lower,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom="segment", x = LT1_hr_results$orig_upper, xend = LT1_hr_results$orig_upper,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom = "text", x = x_middle, y = ymax, label = range_text, hjust=0.5, vjust=0.5, size = range_fnt) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = below_df,
              fill = "white", alpha = outer_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = above_df,
              fill = "white", alpha = outer_alf) + 
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0,ymax_plot),
                     expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  labs(caption = "RunningWritings.com") + 
  ggtitle(title_string,
          subtitle = subtitle_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        plot.subtitle = element_text(size=sub_fnt, hjust=0.5),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

total_plot


ggsave("01 - Heart rate at LT1 in runners - individual variation.png",
       total_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")


# -- Now generate for individul vs. lab vs. total
title_string <- sprintf("True individual variation: %.0f%% of the total variation", 100*LT1_hr_results$individual_fraction)
ribbon_fill <- "#38bdf8"

# -- Individual component
ind_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_individual_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_individual_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1))
ind_plot


# -- Lab component
title_string <- sprintf("Lab-to-lab variation: %.0f%% of the total variation", 100*LT1_hr_results$lab_fraction)
ribbon_fill <- "#818cf8"

lab_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_lab_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_lab_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
        )

lab_plot



# -- Total (again)
title_string <- "Total variation (individual + lab-to-lab)"
ribbon_fill <- hr_fill

total_plot_2 <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
        )

total_plot_2

hr_stack_plot <- plot_grid(ind_plot, lab_plot, total_plot_2, ncol=1, rel_heights = c(43,57,100), scale = 0.95)
hr_stack_plot



# Ughh...HOW do you do title + caption again?


sup_fnt <- 20
cap_fnt <- 14


make_title <- function(title_string){
  title_1 <- ggdraw() + 
    draw_label(
      title_string,
      fontface = 'bold', x = 0.5, y=0.6, hjust = 0.5, size=sup_fnt) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 2, 0, 2)
    )
  return(title_1)
  
}

#plot margin is TOP, RIGHT, BOTTOM, LEFT
caption <- ggdraw() + 
  draw_label(
    "RunningWritings.com",
    fontface = 'bold', x = 0.97, hjust = 1, size=cap_fnt) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 2, 0, 2)
  )


super_combo_hr <- plot_grid(make_title("Heart rate at LT1: Sources of variation"), 
                            hr_stack_plot, 
                            caption, 
                            ncol=1, rel_heights = c(1,10,1), scale=1.0)
super_combo_hr


ggsave("02 - Sources of variation in heart rate at LT1 in runners.png",
       super_combo_hr, width = 2*1024, height = 1024*2, units="px",
       bg = "white")


LT1_plot_df <- plot_df


# ---- LT2 and HRmax

LT2_hr_results <- fit_log_ml_meta_model(df, training_metric = "pct_hrmax", threshold_metric = "LT2", conf_level = 0.9)

I2_res <- get_generalized_I2(LT2_hr_results$meta_model)
print_generalized_I2(I2_res)

make_profile_plot(LT2_hr_results)
make_forest(LT2_hr_results)


LT2_hr_plot_res <- make_plot_df(LT2_hr_results)


plot_df <- LT2_hr_plot_res$plot_df
below_df <- LT2_hr_plot_res$below_df
above_df <- LT2_hr_plot_res$above_df


plot_df %>% glimpse()

plot_df %>%
  ggplot(aes(x=pct_metric, y=dist_total_orig)) + 
  geom_line() + 
  geom_line(aes(y=dist_individual_orig), color = "blue") + 
  geom_line(aes(y = dist_lab_orig), color = "green")



sum(plot_df$dist_individual_orig)
sum(plot_df$dist_total_orig)


# --- Plot -----

ymax <- max(plot_df$dist_total_orig)*1.08
ymax_plot <- ymax*1.1

ant_col <- "#94a3b8"
ant_lwd <- 0.75
ant_type <- "solid"


hr_fill <- "#d946ef"
LT2_col <- hr_fill

ribbon_fill <- hr_fill
ribbon_line_col <- "black"
ribbon_lwd <- 0.5

#Axis setup
x_limits <- c(50,100)
x_breaks <- seq(50,100,by=10)
x_labels <- sprintf("%d%%", x_breaks)
x_name <- "Percent of HRmax"

range_text <- sprintf("90%% range: %.0f - %.0f%% HRmax",
                      LT2_hr_results$orig_lower,
                      LT2_hr_results$orig_upper)

title_string <- "LT2 as a percentage of HRmax in runners"


subtitle_string <- sprintf("Pooled from %d studies on a total of %d runners",
                           length(LT2_hr_results$analysis_df$study_cluster),
                           sum(LT2_hr_results$analysis_df$n_subjects))


x_middle <- (LT2_hr_results$orig_lower + LT2_hr_results$orig_upper)/2


single_title_fnt <- 18
range_fnt <- 5

#Now try with fill
total_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  #Lower and upper shade
  #Upper/lower
  annotate(geom="segment", x = LT2_hr_results$orig_lower, xend = LT2_hr_results$orig_lower,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom="segment", x = LT2_hr_results$orig_upper, xend = LT2_hr_results$orig_upper,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom = "text", x = x_middle, y = ymax*1.04, label = range_text, hjust=0.5, vjust=0.5, size = range_fnt) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = below_df,
              fill = "white", alpha = outer_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = above_df,
              fill = "white", alpha = outer_alf) + 
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0,ymax_plot),
                     expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  labs(caption = "RunningWritings.com") + 
  ggtitle(title_string,
          subtitle = subtitle_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size=sub_fnt, hjust=0.5))

total_plot


ggsave("03 - Heart rate at LT2 in runners - individual variation.png",
       total_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")




# -- Now generate for individul vs. lab vs. total
title_string <- sprintf("True individual variation: %.0f%% of the total variation", 100*LT2_hr_results$individual_fraction)
ribbon_fill <- "#38bdf8"

# -- Individual component
ind_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_individual_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_individual_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1))
ind_plot


# -- Lab component
title_string <- sprintf("Lab-to-lab variation: %.0f%% of the total variation", 100*LT2_hr_results$lab_fraction)
ribbon_fill <- "#818cf8"

lab_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_lab_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_lab_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
  )

lab_plot



# -- Total (again)
title_string <- "Total variation (individual + lab-to-lab)"
ribbon_fill <- hr_fill

total_plot_2 <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
  )

total_plot_2

hr_stack_plot <- plot_grid(ind_plot, lab_plot, total_plot_2, ncol=1, rel_heights = c(70,30,100), scale = 0.95)
hr_stack_plot



# Ughh...HOW do you do title + caption again?


make_title <- function(title_string){
  title_1 <- ggdraw() + 
    draw_label(
      title_string,
      fontface = 'bold', x = 0.5, y=0.6, hjust = 0.5, size=sup_fnt) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 2, 0, 2)
    )
  return(title_1)
  
}

#plot margin is TOP, RIGHT, BOTTOM, LEFT
caption <- ggdraw() + 
  draw_label(
    "RunningWritings.com",
    fontface = 'bold', x = 0.97, hjust = 1, size=cap_fnt) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 2, 0, 2)
  )


super_combo_hr <- plot_grid(make_title("Heart rate at LT2: Sources of variation"), 
                            hr_stack_plot, 
                            caption, 
                            ncol=1, rel_heights = c(1,10,1), scale=1.0)
super_combo_hr


ggsave("04 - Sources of variation in heart rate at LT2 in runners.png",
       super_combo_hr, width = 2*1024, height = 1024*2, units="px",
       bg = "white")




# ------- LT1 / LT2 overlap ------------------

single_title_fnt <- 16
title_string <- "Substantial overlap in LT1 and LT2 heart rate ranges"


new_alf <- 0.5

leg_fnt <- 14

overlap_plot <- ggplot() +
  # LT1 data with mapped fill label "LT1"
  geom_ribbon(data = LT1_plot_df, 
              aes(x = pct_metric, ymin = 0, ymax = dist_total_orig, fill = "LT1"), 
              alpha = new_alf) + 
  geom_line(data = LT1_plot_df,
            aes(x = pct_metric, y = dist_total_orig),
            linewidth = ribbon_lwd, color = ribbon_line_col) + 
  # LT2 data with mapped fill label "LT2"
  geom_ribbon(data = plot_df, 
              aes(x = pct_metric, ymin = 0, ymax = dist_total_orig, fill = "LT2"), 
              alpha = new_alf) + 
  geom_line(data = plot_df,
            aes(x = pct_metric, y = dist_total_orig),
            linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0, ymax_plot),
                     expand = c(0, 0)) + 
  coord_cartesian(clip = 'off') + 
  guides(x = guide_axis(cap = "both"),
         fill = guide_legend(title.position = "top", title.hjust=0.5)) +  
  labs(caption = "RunningWritings.com",
       title = title_string,
       fill = "Distribution across runners") +  # Label for the fill legend
  scale_fill_manual(values = c("LT1" = LT1_col, "LT2" = ribbon_fill)) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        plot.subtitle = element_text(size=sub_fnt, hjust = 0.5),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = leg_fnt),
        legend.title = element_text(size=leg_fnt, hjust=0.5),
        #legend.position = c(0.1, 0.9),
        #legend.justification = c(0, 1),
        legend.position = "inside", 
        legend.position.inside =  c(.2, .85),
        legend.direction = "horizontal",
        legend.box = "vertical",  # Stack title above keys
        ) 
  # Manually set the fill colors for each label

overlap_plot


ggsave("04 - Heart rate at LT1 vs. LT2 in runners - extensive overlap.png",
       overlap_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")


#Sums to ~ the same!
sum(LT1_plot_df$dist_total_orig)
sum(plot_df$dist_total_orig)




# TODO: HRR and VO2max plots

# --------- LT1 and HRR -----------------

#It's really just Weltman data so manually implement single study analysis 


#df_hrr

welt_LT1 <- df %>% filter(percent_metric == "pct_hrreserve",
              metric == "LT1")


hrr_z <- 1

#REcall thsi is on transformed scale
hrr_mean <- welt_LT1$ybarstar
hrr_sd <- welt_LT1$si_star


conf_level <- 0.9
upper_tail <- (1+conf_level)/2
lower_tail <- 1 - upper_tail

#Weltman range only (transformeds cale)
lower_ref <- qnorm(lower_tail, mean = hrr_mean, sd = hrr_sd)
upper_ref <- qnorm(upper_tail, mean = hrr_mean, sd = hrr_sd)

#INvert to PCT HRresrve scale
orig_mu <- 100 - exp(hrr_mean)
orig_lower <- 100 - exp(upper_ref) #note flipped
orig_upper <- 100 - exp(lower_ref)


#Grid up the dnorm then invert w/ jacobian correction
ng <- 501
pct_metric_grid <- seq(50,99.9, length.out=ng)
log_x_grid <- log(100 - pct_metric_grid)
dist_star = dnorm(log_x_grid, mean = hrr_mean, sd = hrr_sd)
dist_orig = dist_star * (1 / (100 - pct_metric_grid)) #<-- Jacobian correction


plot_df <- data.frame(pct_metric = pct_metric_grid,
                      log_pct_metric = log_x_grid,
                      dist_star = dist_star,
                      dist_orig = dist_orig)


below_df <- plot_df %>%
  filter(pct_metric < orig_lower)

above_df <- plot_df %>%
  filter(pct_metric > orig_upper)

# --- Now plot LT1 as HRR




ymax <- max(dist_orig)*1.08
ymax_plot <- ymax*1.1

ant_col <- "#94a3b8"
ant_lwd <- 0.75
ant_type <- "solid"


hr_fill <- "#1e40af"
LT1_col <- hr_fill

ribbon_fill <- hr_fill
ribbon_line_col <- "black"
ribbon_lwd <- 0.5

#Axis setup
x_limits <- c(50,100)
x_breaks <- seq(50,100,by=10)
x_labels <- sprintf("%d%%", x_breaks)
x_name <- "Percent of HR reserve"

range_text <- sprintf("90%% range: %.0f - %.0f%% HR reserve",
                      orig_lower,
                      orig_upper)

title_string <- "LT1 as a percentage of HR reserve in runners (version 1)"


subtitle_string <- "Data from one study on 31 runners"


x_middle <- (orig_lower + orig_upper)/2


single_title_fnt <- 18
range_fnt <- 5

#Now try with fill
total_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_orig)) +
  #Lower and upper shade
  #Upper/lower
  annotate(geom="segment", x = orig_lower, xend = orig_lower,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom="segment", x = orig_upper, xend = orig_upper,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom = "text", x = x_middle, y = ymax*1.04, label = range_text, hjust=0.5, vjust=0.5, size = range_fnt) + 
  geom_ribbon(aes(ymin=0, ymax = dist_orig), fill = ribbon_fill, alpha = within_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_orig), data = below_df,
              fill = "white", alpha = outer_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_orig), data = above_df,
              fill = "white", alpha = outer_alf) + 
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0,ymax_plot),
                     expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  labs(caption = "RunningWritings.com") + 
  ggtitle(title_string,
          subtitle = subtitle_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size=sub_fnt, hjust=0.5))

total_plot


ggsave("05 - Heart rate reserve at LT1 in runners - individual variation.png",
       total_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")

welt_LT1_plot_df <- plot_df


# --------- LT2 and HRR -----------------

#It's really just Weltman data so manually implement single study analysis 


#df_hrr

welt_LT2 <- df %>% filter(percent_metric == "pct_hrreserve",
                          metric == "LT2")


hrr_z <- 1

#REcall thsi is on transformed scale
hrr_mean <- welt_LT2$ybarstar
hrr_sd <- welt_LT2$si_star


conf_level <- 0.9
upper_tail <- (1+conf_level)/2
lower_tail <- 1 - upper_tail

#Weltman range only (transformeds cale)
lower_ref <- qnorm(lower_tail, mean = hrr_mean, sd = hrr_sd)
upper_ref <- qnorm(upper_tail, mean = hrr_mean, sd = hrr_sd)

#INvert to PCT HRresrve scale
orig_mu <- 100 - exp(hrr_mean)
orig_lower <- 100 - exp(upper_ref) #note flipped
orig_upper <- 100 - exp(lower_ref)


#Grid up the dnorm then invert w/ jacobian correction
ng <- 501
pct_metric_grid <- seq(50,99.9, length.out=ng)
log_x_grid <- log(100 - pct_metric_grid)
dist_star = dnorm(log_x_grid, mean = hrr_mean, sd = hrr_sd)
dist_orig = dist_star * (1 / (100 - pct_metric_grid)) #<-- Jacobian correction


plot_df <- data.frame(pct_metric = pct_metric_grid,
                      log_pct_metric = log_x_grid,
                      dist_star = dist_star,
                      dist_orig = dist_orig)


below_df <- plot_df %>%
  filter(pct_metric < orig_lower)

above_df <- plot_df %>%
  filter(pct_metric > orig_upper)

# --- Now plot LT2 as HRR


ymax <- max(dist_orig)*1.08
ymax_plot <- ymax*1.1

ant_col <- "#94a3b8"
ant_lwd <- 0.75
ant_type <- "solid"


hr_fill <- "#b91c1c"
LT2_col <- hr_fill

ribbon_fill <- hr_fill
ribbon_line_col <- "black"
ribbon_lwd <- 0.5

#Axis setup
x_limits <- c(50,100)
x_breaks <- seq(50,100,by=10)
x_labels <- sprintf("%d%%", x_breaks)
x_name <- "Percent of HR reserve"

range_text <- sprintf("90%% range: %.0f - %.0f%% HR reserve",
                      orig_lower,
                      orig_upper)

title_string <- "LT2 as a percentage of HR reserve in runners (version 1)"


subtitle_string <- "Data from one study on 31 runners"


x_middle <- (orig_lower + orig_upper)/2


single_title_fnt <- 18
range_fnt <- 5

#Now try with fill
total_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_orig)) +
  #Lower and upper shade
  #Upper/lower
  annotate(geom="segment", x = orig_lower, xend = orig_lower,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom="segment", x = orig_upper, xend = orig_upper,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom = "text", x = x_middle, y = ymax*1.04, label = range_text, hjust=0.5, vjust=0.5, size = range_fnt) + 
  geom_ribbon(aes(ymin=0, ymax = dist_orig), fill = ribbon_fill, alpha = within_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_orig), data = below_df,
              fill = "white", alpha = outer_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_orig), data = above_df,
              fill = "white", alpha = outer_alf) + 
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0,ymax_plot),
                     expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  labs(caption = "RunningWritings.com") + 
  ggtitle(title_string,
          subtitle = subtitle_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size=sub_fnt, hjust=0.5))

total_plot


ggsave("06 - Heart rate reserve at LT2 in runners - individual variation.png",
       total_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")

welt_LT2_plot_df <- plot_df


# --- Overlay plot

# ------- LT1 / LT2 overlap ------------------

single_title_fnt <- 16
title_string <- "Substantial overlap in LT1 and LT2 heart rate reserve values"
subtitle_string <- "Caution: Data from only one study on 31 runners"


new_alf <- 0.5

leg_fnt <- 14

overlap_plot <- ggplot() +
  # LT1 data with mapped fill label "LT1"
  geom_ribbon(data = welt_LT1_plot_df, 
              aes(x = pct_metric, ymin = 0, ymax = dist_orig, fill = "LT1"), 
              alpha = new_alf) + 
  geom_line(data = welt_LT1_plot_df,
            aes(x = pct_metric, y = dist_orig),
            linewidth = ribbon_lwd, color = ribbon_line_col) + 
  # LT2 data with mapped fill label "LT2"
  geom_ribbon(data = welt_LT2_plot_df, 
              aes(x = pct_metric, ymin = 0, ymax = dist_orig, fill = "LT2"), 
              alpha = new_alf) + 
  geom_line(data = welt_LT2_plot_df,
            aes(x = pct_metric, y = dist_orig),
            linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0, ymax_plot),
                     expand = c(0, 0)) + 
  coord_cartesian(clip = 'off') + 
  guides(x = guide_axis(cap = "both"),
         fill = guide_legend(title.position = "top", title.hjust=0.5)) +  
  labs(caption = "RunningWritings.com",
       title = title_string,
       subtitle = subtitle_string,
       fill = "Distribution across runners") +  # Label for the fill legend
  scale_fill_manual(values = c("LT1" = LT1_col, "LT2" = LT2_col)) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        plot.subtitle = element_text(size=sub_fnt, hjust = 0.5),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = leg_fnt),
        legend.title = element_text(size=leg_fnt, hjust=0.5),
        #legend.position = c(0.1, 0.9),
        #legend.justification = c(0, 1),
        legend.position = "inside", 
        legend.position.inside =  c(.2, .85),
        legend.direction = "horizontal",
        legend.box = "vertical",  # Stack title above keys
  ) 
# Manually set the fill colors for each label

overlap_plot


ggsave("07 - Heart rate reserve overlap at LT1 and LT2 in runners.png",
       overlap_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")



# ----------- VO2max ----------------------


# analyze ybarstar, si_star as mean, sd
LT1_vo2_results <- fit_log_ml_meta_model(df, training_metric = "pct_vo2max", threshold_metric = "LT1", conf_level = 0.9)

#fit_log_ml_meta_model

make_profile_plot(LT1_vo2_results)
make_forest(LT1_vo2_results)

I2_res <- get_generalized_I2(LT1_vo2_results$meta_model)
print_generalized_I2(I2_res)

LT1_vo2_plot_res <- make_plot_df(LT1_vo2_results)


plot_df <- LT1_vo2_plot_res$plot_df
below_df <- LT1_vo2_plot_res$below_df
above_df <- LT1_vo2_plot_res$above_df



sum(plot_df$dist_individual_orig)
sum(plot_df$dist_total_orig)


plot_df %>% glimpse()

plot_df %>%
  ggplot(aes(x=pct_metric, y=dist_total_orig)) + 
  geom_line() + 
  geom_line(aes(y=dist_individual_orig), color = "blue") + 
  geom_line(aes(y = dist_lab_orig), color = "green")




ymax <- max(plot_df$dist_total_orig)*1.08
ymax_plot <- ymax*1.1

ant_col <- "#94a3b8"
ant_lwd <- 0.75
ant_type <- "solid"


hr_fill <- "#14b8a6"
LT1_col <- hr_fill

ribbon_fill <- hr_fill
ribbon_line_col <- "black"
ribbon_lwd <- 0.5

#Axis setup
x_limits <- c(50,100)
x_breaks <- seq(50,100,by=10)
x_labels <- sprintf("%d%%", x_breaks)
x_name <- "Percent of VO2max (and arguably, heart rate reserve)"

range_text <- sprintf("90%% range: %.0f - %.0f%% VO2max / %%HRreserve",
                      LT1_vo2_results$orig_lower,
                      LT1_vo2_results$orig_upper)

title_string <- "LT1 as a percentage of HR reserve in runners (version 2)"


subtitle_string <- sprintf("Pooled from %d studies on a total of %d runners, using VO2max-HRR correspondence",
                           length(LT1_vo2_results$analysis_df$study_cluster),
                           sum(LT1_vo2_results$analysis_df$n_subjects))


x_middle <- (LT1_vo2_results$orig_lower + LT1_vo2_results$orig_upper)/2


single_title_fnt <- 18
range_fnt <- 5
sub_fnt <- 12

#Now try with fill
total_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  #Lower and upper shade
  #Upper/lower
  annotate(geom="segment", x = LT1_vo2_results$orig_lower, xend = LT1_vo2_results$orig_lower,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom="segment", x = LT1_vo2_results$orig_upper, xend = LT1_vo2_results$orig_upper,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom = "text", x = x_middle, y = ymax*1.04, label = range_text, hjust=0.5, vjust=0.5, size = range_fnt) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = below_df,
              fill = "white", alpha = outer_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = above_df,
              fill = "white", alpha = outer_alf) + 
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0,ymax_plot),
                     expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  labs(caption = "RunningWritings.com") + 
  ggtitle(title_string,
          subtitle = subtitle_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size=sub_fnt, hjust=0.5))

total_plot


ggsave("08 - Heart rate reserve at LT1 in runners - VO2max version - individual variation.png",
       total_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")

LT1_vo2_plot_df <- plot_df


# --- Now let's break it down


# -- Now generate for individul vs. lab vs. total
title_string <- sprintf("True individual variation: %.0f%% of the total variation", 100*LT1_vo2_results$individual_fraction)
ribbon_fill <- "#38bdf8"

# -- Individual component
ind_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_individual_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_individual_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1))
ind_plot


# -- Lab component
title_string <- sprintf("Lab-to-lab variation: %.0f%% of the total variation", 100*LT1_vo2_results$lab_fraction)
ribbon_fill <- "#818cf8"

lab_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_lab_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_lab_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
  )

lab_plot



# -- Total (again)
title_string <- "Total variation (individual + lab-to-lab)"
ribbon_fill <- hr_fill

total_plot_2 <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
  )

total_plot_2


vo2_LT1_stack_plot <- plot_grid(ind_plot, lab_plot, total_plot_2, ncol=1, rel_heights = c(80,20,100), scale = 0.95)
vo2_LT1_stack_plot

super_combo_vo2 <- plot_grid(make_title("VO2max (and HRR) at LT1: Sources of variation"), 
                            vo2_LT1_stack_plot, 
                            caption, 
                            ncol=1, rel_heights = c(1,10,1), scale=1.0)
super_combo_vo2


ggsave("09 - Sources of variation in heart rate reserve and vo2max at LT1 in runners.png",
       super_combo_vo2, width = 2*1024, height = 1024*2, units="px",
       bg = "white")


# ---- Samesies for LT2 at VO2max pct ------------


# analyze ybarstar, si_star as mean, sd
LT2_vo2_results <- fit_log_ml_meta_model(df, training_metric = "pct_vo2max", threshold_metric = "LT2", conf_level = 0.9)

#fit_log_ml_meta_model

make_profile_plot(LT2_vo2_results)
make_forest(LT2_vo2_results)

I2_res <- get_generalized_I2(LT2_vo2_results$meta_model)
print_generalized_I2(I2_res)

LT2_vo2_plot_res <- make_plot_df(LT2_vo2_results)


plot_df <- LT2_vo2_plot_res$plot_df
below_df <- LT2_vo2_plot_res$below_df
above_df <- LT2_vo2_plot_res$above_df



sum(plot_df$dist_individual_orig)
sum(plot_df$dist_total_orig)


plot_df %>% glimpse()

plot_df %>%
  ggplot(aes(x=pct_metric, y=dist_total_orig)) + 
  geom_line() + 
  geom_line(aes(y=dist_individual_orig), color = "blue") + 
  geom_line(aes(y = dist_lab_orig), color = "green")




ymax <- max(plot_df$dist_total_orig)*1.08
ymax_plot <- ymax*1.1

ant_col <- "#94a3b8"
ant_lwd <- 0.75
ant_type <- "solid"


hr_fill <- "#6d28d9"
LT2_col <- hr_fill

ribbon_fill <- hr_fill
ribbon_line_col <- "black"
ribbon_lwd <- 0.5

#Axis setup
x_limits <- c(50,100)
x_breaks <- seq(50,100,by=10)
x_labels <- sprintf("%d%%", x_breaks)
x_name <- "Percent of VO2max (and arguably, heart rate reserve)"

range_text <- sprintf("90%% range: %.0f - %.0f%% VO2max / %%HRreserve",
                      LT2_vo2_results$orig_lower,
                      LT2_vo2_results$orig_upper)

title_string <- "LT2 as a percentage of HR reserve in runners (version 2)"


subtitle_string <- sprintf("Pooled from %d studies on a total of %d runners, using VO2max-HRR correspondence",
                           length(LT2_vo2_results$analysis_df$study_cluster),
                           sum(LT2_vo2_results$analysis_df$n_subjects))


x_middle <- (LT2_vo2_results$orig_lower + LT2_vo2_results$orig_upper)/2


single_title_fnt <- 18
range_fnt <- 5
sub_fnt <- 12

#Now try with fill
total_plot <- plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  #Lower and upper shade
  #Upper/lower
  annotate(geom="segment", x = LT2_vo2_results$orig_lower, xend = LT2_vo2_results$orig_lower,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom="segment", x = LT2_vo2_results$orig_upper, xend = LT2_vo2_results$orig_upper,
           y = 0, yend = ymax, linewidth = ant_lwd, linetype = ant_type, color = ant_col) + 
  annotate(geom = "text", x = x_middle, y = ymax*1.04, label = range_text, hjust=0.5, vjust=0.5, size = range_fnt) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = below_df,
              fill = "white", alpha = outer_alf) + 
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), data = above_df,
              fill = "white", alpha = outer_alf) + 
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0,ymax_plot),
                     expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  labs(caption = "RunningWritings.com") + 
  ggtitle(title_string,
          subtitle = subtitle_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size=sub_fnt, hjust=0.5))

total_plot
LT2_vo2_plot_df <- plot_df


ggsave("10 - Heart rate reserve at LT2 in runners - VO2max version - individual variation.png",
       total_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")



# ------- LT1 / LT2 overlap ------------------

single_title_fnt <- 16
title_string <- "Substantial overlap in LT1 and LT2 VO2 (and HRR) values"
subtitle_string <- "Caution: Data from only one study on 31 runners"


new_alf <- 0.5

leg_fnt <- 14

overlap_plot <- ggplot() +
  # LT1 data with mapped fill label "LT1"
  geom_ribbon(data = LT1_vo2_plot_df, 
              aes(x = pct_metric, ymin = 0, ymax = dist_total_orig, fill = "LT1"), 
              alpha = new_alf) + 
  geom_line(data = LT1_vo2_plot_df,
            aes(x = pct_metric, y = dist_total_orig),
            linewidth = ribbon_lwd, color = ribbon_line_col) + 
  # LT2 data with mapped fill label "LT2"
  geom_ribbon(data = LT2_vo2_plot_df, 
              aes(x = pct_metric, ymin = 0, ymax = dist_total_orig, fill = "LT2"), 
              alpha = new_alf) + 
  geom_line(data = LT2_vo2_plot_df,
            aes(x = pct_metric, y = dist_total_orig),
            linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(limits = c(0, ymax_plot),
                     expand = c(0, 0)) + 
  coord_cartesian(clip = 'off') + 
  guides(x = guide_axis(cap = "both"),
         fill = guide_legend(title.position = "top", title.hjust=0.5)) +  
  labs(caption = "RunningWritings.com",
       title = title_string,
       fill = "Distribution across runners") +  # Label for the fill legend
  scale_fill_manual(values = c("LT1" = LT1_col, "LT2" = LT2_col)) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        plot.title = element_text(size = single_title_fnt, face="bold"),
        plot.subtitle = element_text(size=sub_fnt, hjust = 0.5),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = leg_fnt),
        legend.title = element_text(size=leg_fnt, hjust=0.5),
        #legend.position = c(0.1, 0.9),
        #legend.justification = c(0, 1),
        legend.position = "inside", 
        legend.position.inside =  c(.2, .85),
        legend.direction = "horizontal",
        legend.box = "vertical",  # Stack title above keys
  ) 
# Manually set the fill colors for each label

overlap_plot


ggsave("11 - VO2max percentage overlap at LT1 and LT2 in runners.png",
       overlap_plot, width = 2*1024, height = 1536, units="px",
       bg = "white")


# --- Finally, ONE more breakdown -----



# -- Now generate for individul vs. lab vs. total
title_string <- sprintf("True individual variation: %.0f%% of the total variation", 100*LT2_vo2_results$individual_fraction)
ribbon_fill <- "#38bdf8"

# -- Individual component
ind_plot <- LT2_vo2_plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_individual_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_individual_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1))
ind_plot


# -- Lab component
title_string <- sprintf("Lab-to-lab variation: %.0f%% of the total variation", 100*LT2_vo2_results$lab_fraction)
ribbon_fill <- "#818cf8"

lab_plot <- LT2_vo2_plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_lab_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_lab_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
  )

lab_plot



# -- Total (again)
title_string <- "Total variation (individual + lab-to-lab)"
ribbon_fill <- hr_fill

total_plot_2 <- LT2_vo2_plot_df %>% 
  ggplot(aes(x=pct_metric, y=dist_total_orig)) +
  geom_ribbon(aes(ymin=0, ymax = dist_total_orig), fill = ribbon_fill, alpha = within_alf) + 
  #Lower and upper shade
  geom_line(linewidth = ribbon_lwd, color = ribbon_line_col) + 
  scale_x_continuous(limits = x_limits,
                     breaks = x_breaks,
                     labels = x_labels,
                     name = x_name) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(clip='off') + 
  guides(x = guide_axis(cap = "both")) +
  ggtitle(title_string) + 
  theme_bw() + 
  theme_mee + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust=1)
  )

total_plot_2


vo2_LT2_stack_plot <- plot_grid(ind_plot, lab_plot, total_plot_2, ncol=1, rel_heights = c(60,40,100), scale = 0.95)
vo2_LT2_stack_plot

super_combo_vo2_2 <- plot_grid(make_title("VO2max (and HRR) at LT2: Sources of variation"), 
                               vo2_LT2_stack_plot, 
                             caption, 
                             ncol=1, rel_heights = c(1,10,1), scale=1.0)
super_combo_vo2_2


ggsave("12 - Sources of variation in heart rate reserve and vo2max at LT2 in runners.png",
       super_combo_vo2_2, width = 2*1024, height = 1024*2, units="px",
       bg = "white")



#N studies, N subjects

df %>%
  group_by(study) %>%
  slice(1) %>%
  pull(n_subjects) %>% sum()


unique(df$study) %>% length()