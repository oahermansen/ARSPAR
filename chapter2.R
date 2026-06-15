## Initial Investigation Script
## Using master_functions.R

# Implementation notes
# Load data
source("SPAR-rcode/master_functions.R")
load("data/ekofisk_wave_surge.rda")
library(dplyr)
library(ggplot2)
# Plot timeseries
# Select parameters
theme_set(theme_gray(ink="blue"))
ekofisk_tm2_hs = ekofisk_wave_surge[c("tm2", "hs", "doy")]
zero_t <- data.frame(tm2=mean(ekofisk_tm2_hs$tm2), hs=mean(ekofisk_tm2_hs$hs))
inv_s <- data.frame(tm2 = sd(ekofisk_tm2_hs$tm2), hs=sd(ekofisk_tm2_hs$hs))

ekofisk_tscaled <- ekofisk_tm2_hs |> mutate(tm2=(tm2-mean(tm2))/sd(tm2), hs = (hs-mean(hs))/sd(hs))
# SPAR transforms an
k <- 25
k_shape <- 35
pred_q <- seq_between(c(-2, 2), length=2001)
nna_ind <- pmax(abs(ekofisk_tscaled$tm2), abs(ekofisk_tscaled$hs)) > 0
ekofisk_tscaled <- ekofisk_tscaled[nna_ind, ]

l1_radial <- function(data) {
  abs(data[, "tm2"]) + abs(data[, "hs"])
}
l1_angle <- function(data) {
  s <- data[, "R_l1"]
  a <- sign(data[, "hs"])
  a[a == 0] <- 1
  out <- a * (1 - data[, "tm2"] / s)
  out[s == 0] <- NA_real_
  out
}
l2_radial <- function(data){
  return(sqrt(data[,"tm2"]^2+data[,"hs"]^2))
}
l2_angle <- function(data){
  return(sign(data[,"hs"])*(2/pi)*acos(data[,"tm2"]/data[,"R_l2"]))
}

ekofisk_tscaled$R_l1 <- l1_radial(ekofisk_tscaled)
ekofisk_tscaled$phi_l1 <- l1_angle(ekofisk_tscaled)
ekofisk_tscaled$R_l2 <- l2_radial(ekofisk_tscaled)
ekofisk_tscaled$phi_l2 <- l2_angle(ekofisk_tscaled)

# Threshold fits
ald_fits <- list(l1 = list(), l2 = list())
ekofisk_tscaled$logR_l1 <- log(ekofisk_tscaled$R_l1)
ekofisk_tscaled$logR_l2 <- log(ekofisk_tscaled$R_l2)
levels <- c(0.8, 0.85, 0.9, 0.95)
names(levels) <- c("Tau = 0.80", "Tau = 0.85", "Tau = 0.90", "Tau = 0.95")

excess <- list(l1 = list(), l2 = list())

threshold_function <- function(data, gauge, level){
  name <- which(levels == level)
  fit <- ald_fits[[gauge]][[name]]
  return(exp(predict(fit, newdata=data)$location))
}

k_ald_l1 <- c(25, 25, 25, 35)
k_ald_l2 <- c(35, 35, 35, 35)
ald_fmla <- function(gauge, k) {
  return(list(as.formula(paste0("logR_", gauge, " ~ s(phi_", gauge, ", bs='cc', k=", k, ")"))))
}
ald_fmla_l1 <- list(as.formula(paste0("logR_l1 ~ s(phi_l1, bs='cc', k=", k, ")")))
ald_fmla_l2 <- list(as.formula(paste0("logR_l2 ~ s(phi_l2, bs='cc', k=", k, ")")))
threshold_grid <- list(l1 = list(), l2 = list())
for(i in 1:length(levels)){
  level <- names(levels)[i]
  print(paste("Fitting L1 threshold distribution", level))
  ald_fits$l1[[i]] <- evgam(ald_fmla("l1", k_ald_l1[i]), data = ekofisk_tscaled, family = "ald", ald.args = list(tau = levels[i]), trace = 2)
  names(ald_fits$l1)[i] <- level
  
  print(paste("Predicting L1 excess", level))
  excess$l1[[i]] <- ekofisk_tscaled$R_l1 - threshold_function(ekofisk_tscaled, "l1", levels[i])
  is.na(excess$l1[[i]]) <- excess$l1[[i]] < 0
  names(excess$l1)[i] <- level
  threshold_grid$l1[[i]] <- threshold_function(data.frame(phi_l1 = pred_q), "l1", levels[i])
  names(threshold_grid$l1)[i] <- level
}

for(i in 1:length(levels)){
  level <- names(levels)[i]
  print(paste("Fitting L2 threshold distribution", level))
  ald_fits$l2[[i]] <- evgam(ald_fmla("l2", k_ald_l2[i]), data = ekofisk_tscaled, family = "ald", ald.args = list(tau = levels[i]), trace = 2)
  names(ald_fits$l2)[i] <- level
  
  print(paste("Predicting L2 excess", level))
  excess$l2[[i]] <- ekofisk_tscaled$R_l2 - threshold_function(ekofisk_tscaled, "l2", levels[i])
  is.na(excess$l2[[i]]) <- excess$l2[[i]] < 0
  names(excess$l2)[i] <- level
  threshold_grid$l2[[i]] <- threshold_function(data.frame(phi_l2 = pred_q), "l2", levels[i])
  names(threshold_grid$l2)[i] <- level
}

# Exceedances
for(i in 1:length(levels)){
  ekofisk_tscaled[, paste0("exc_", levels[i]*100, "_l1")] <- excess$l1[[i]]
  ekofisk_tscaled[, paste0("exc_", levels[i]*100, "_l2")] <- excess$l2[[i]]
}


# Gauge vectors from https://github.com/callumbarltrop/SPAR/
# x = R * u, y = R * v
# u and v are the x and y value corresponding to each angular value in the grid projected onto the norm unit circle
gauge_vecs <- list(l1 = list(u = ifelse(pred_q>=0,(1-pred_q),(pred_q+1))), 
                   l2 = list(u = ifelse(pred_q>=0,cos(pi*pred_q/2),cos(-pi*pred_q/2))))
gauge_vecs$l1$v <- ifelse(pred_q>=0, 1-abs(gauge_vecs$l1$u),-1+abs(gauge_vecs$l1$u))
gauge_vecs$l2$v <- ifelse(pred_q>=0, sqrt(1-gauge_vecs$l2$u^2),-sqrt(1-gauge_vecs$l2$u^2))

plot_threshold <- function(gauge, level){
  plot <- ggplot(data = ekofisk_tscaled) + geom_point(mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs, color = .data[[paste0("exc_",levels[level]*100,"_",gauge)]] > 0), size = 0.1) +
    geom_path(data=data.frame(tm2 = gauge_vecs[[gauge]]$u*threshold_grid[[gauge]][[level]], hs = gauge_vecs[[gauge]]$v*threshold_grid[[gauge]][[level]]), mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs), color="red") +
    labs(x = "tm2", y = "hs", title = paste("Threshold level contour,", gauge_names[gauge], "gauge,", level), color="Exceedance")
  ggsave(paste0("plots/chapter 2/threshold/", format(levels[level], nsmall = 2), "_", gauge, ".png"), plot=plot)
}

gauges <- c("l1", "l2")
gauge_names <- c("L1", "L2")
names(gauge_names) <- gauges

for(gauge in gauges){
  for(level in names(levels)){
    plot_threshold(gauge, level)
  }
}

# GPD fits
gpd_fmla <- function(gauge, level) {
  return(list(as.formula(paste0("exc_", level*100, "_", gauge, " ~ s(phi_", gauge, ", bs='cc', k=", k, ")")),
              as.formula(paste0("exc_", level*100, "_", gauge, " ~ s(phi_", gauge, ", bs='cc', k=", k_shape, ")"))))
}
gpd_fits <- list(l1 = list(), l2 = list(), c_scale = list())

for(i in 1:length(levels)){
  level <- names(levels)[i]
  print(paste("Fitting L1 GPD distribution at threshold", level))
  gpd_fits$l1[[i]] <- evgam(gpd_fmla("l1", levels[i]), data = ekofisk_tscaled[!is.na(excess$l1[[i]]),], family = "gpd", trace = 2)
  names(gpd_fits$l1)[i] <- level
  print(paste("Fitting L2 GPD distribution at threshold", level))
  gpd_fits$l2[[i]] <- evgam(gpd_fmla("l2", levels[i]), data = ekofisk_tscaled[!is.na(excess$l2[[i]]),], family = "gpd", trace = 2)
  names(gpd_fits$l2)[i] <- level
}

gpd_estimates <- list(l1 = list(), l2 = list(), c_scale = list())

for(i in 1:length(levels)){
  level <- names(levels)[i]
  print(paste("Estimating GPD parameter estimates at threshold", level))
  for(gauge in gauges){
    tempdf <- data.frame(pred_q)
    ang_name <- paste0("phi_", gauge)
    exc_name <- paste0("exc_", levels[i]*100, "_", gauge)
    colnames(tempdf) <- c(ang_name)
    gpd_estimates[[gauge]][[level]]$phi_grid <- predict(gpd_fits[[gauge]][[level]], newdata = tempdf, type = "response", se.fit = TRUE)
    gpd_estimates[[gauge]][[level]]$phi_grid$phi <- as.numeric(tempdf[, ang_name])
    
    exc_inds <- !is.na(ekofisk_tscaled[, exc_name])
    tempdf <- data.frame(ekofisk_tscaled[exc_inds, ang_name])
    colnames(tempdf) <- c(ang_name)
    gpd_estimates[[gauge]][[level]]$exceedances <- predict(gpd_fits[[gauge]][[level]], newdata = tempdf, type = "response", se.fit = TRUE)
    gpd_estimates[[gauge]][[level]]$exceedances$phi <- as.numeric(tempdf[, ang_name])
  }
}

gpd_parameters_table <- function(){
  table <- data.frame()
  for(gauge in gauges){
    for(level in names(levels)){
      for(param in c("scale", "shape")){
        table <- data.table::rbindlist(list(table, data.frame(response = gpd_estimates[[gauge]][[level]]$phi_grid$fitted[, param],
                                                              phi = pred_q,
                                                              param = rep(param, length(pred_q)),
                                                              gauge = rep(toupper(gauge), length(pred_q)),
                                                              level = rep(level, length(pred_q)))))
      }
    }
  }
  return(table)
}

print("Plotting GPD parameter estimates")
plot <- ggplot(data = gpd_parameters_table()) + geom_path(mapping = aes(x = phi, y = response, color = level)) + facet_grid(param ~ gauge, scales="free") +
  labs(title = "GPD Parameter Estimates", y = "estimate", color = "Threshold level")
ggsave("plots/chapter 2/GPD Parameters/GPD_estimates.png", plot)

# Computing angular densities
print("Estimating angular densitites")

angular_density_grid <- list(l1 = list(), l2 = list())
bw <- 50
for(gauge in gauges){
  phi <- ekofisk_tscaled[, paste0("phi_", gauge)]
  scaled_phi <- (phi+2)*pi/2
  dens_est = density.circular(as.circular(scaled_phi,type="angles",units="radians",template="none",modulo="2pi",zero=0,rotation="counter"), bw=bw,kernel="vonmises")
  dens_est$x[length(dens_est$x)] = 2*pi
  f_phi = approxfun(x = (2*dens_est$x)/pi - 2,y = (pi/2)*dens_est$y) 
  angular_density_grid[[gauge]] <- f_phi(pred_q)
}

qq_plot_regions <- list(list(name = "All Exceedances", min = -2, max = 2),
                        list(name = "First Quadrant", min = 0, max = 1),
                        list(name = "Second Quadrant", min = -1, max = 0),
                        list(name = "Third Quadrant", min = -2, max = -1),
                        list(name = "Fourth Quadrant", min = 1, max = 2),
                        list(name = "Growth Region", min = -0.5, max = 0.7),
                        list(name = "Growth Asymptotic Region", min = 0.625, max = 0.675))
plot_qq_regions <- function(gpd_fit, level, gauge){
  gauge_names <- c("L1", "L2")
  names(gauge_names) <- gauges
  for(region in qq_plot_regions){
    plot.new()
    png(paste(paste0("plots/chapter 2/",region$name, "/", format(levels[level], nsmall=2)), region$name, gauge, "qqplot.png"))
    ang <- ekofisk_tscaled[, paste0("phi_", gauge)]
    indexset <- !is.na(ekofisk_tscaled[, paste0("exc_", levels[level]*100, "_", gauge)]) & ang >= region$min & ang <= region$max
    predict(gpd_fit, newdata = ekofisk_tscaled[indexset,], type = "qqplot")
    title(main = paste(region$name, gauge_names[gauge], "Representation,", level), line = 0.5, cex.main = 0.75)
    dev.off()
  }
}

print(paste("Generating exceedance Q-Q plots"))
for(i in 1:length(levels)){
  name <- names(levels)[i]
  level <- levels[i]
  print(paste("Plotting for", name))
  for(gauge in gauges){
    plot_qq_regions(gpd_fits[[gauge]][[name]], name, gauge)
  }
}

# Isodensity contours (code from https://github.com/callumbarltrop/SPAR/)
density_contours <- list(l1 = list(), l2 = list())
total_return_contours <- list(l1 = list(), l2 = list())
density_levels <- 10^-(3:8)
ret_period_years <- 10
obs_year <- 365.25*24
pred_Q <- pred_q
for(gauge in gauges){
  fits <- gpd_fits[[gauge]]
  fit_levels <- names(fits)
  for(level in fit_levels){
    spar_obj <- list(pred_thresh = threshold_grid[[gauge]][[level]], 
                     pred_para = gpd_estimates[[gauge]][[level]]$phi_grid$fitted,
                     thresh_prob = levels[level],
                     norm_choice = toupper(gauge),
                     pred_Q = pred_q)
    density_contours[[gauge]][[level]] <- SPAR_equidensity_contours(density_levels, spar_obj, angular_density_grid[[gauge]])
    total_return_contours[[gauge]][[level]] <- SPAR_ret_level_sets(ret_period_years, obs_year, spar_obj)
  }
}


plot_density_contours <- function(data, level, gauge){
  contours <- data.table::rbindlist(lapply(data[[gauge]][[level]], as.data.frame), idcol = "Level")
  contours$Level <- format(density_levels[contours$Level])
  plot <- ggplot() + geom_point(data = ekofisk_tscaled, mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs), size = 0.1) +
    geom_path(data = contours, mapping = aes(x = V1 * inv_s$tm2 + zero_t$tm2, y = V2 * inv_s$hs + zero_t$hs, color = Level)) +
    labs(x = "tm2", y = "hs", title = paste("Equidensity contours", gauge_names[gauge], "gauge,", level), color = "Level")
  ggsave(paste0("plots/chapter 2/equidensity contours/", format(levels[level], nsmall = 2), "_", gauge, ".png"), plot=plot)
}

plot_total_return_level <- function(data, level, gauge, ret_period_years){
  plot <- ggplot() + geom_point(data = ekofisk_tscaled, mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs), size = 0.1) +
    geom_path(data = as.data.frame(data[[gauge]][[level]]), mapping = aes(x = V1 * inv_s$tm2 + zero_t$tm2, y = V2 * inv_s$hs + zero_t$hs), color = "red") +
    labs(x = "tm2", y = "hs", title = paste(paste0("Total ", ret_period_years, "-year Return-Level Contour"), gauge_names[gauge], "gauge,", level))
  ggsave(paste0("plots/chapter 2/total return level contours/", format(levels[level], nsmall = 2), "_", gauge, ".png"), plot=plot)
}

print("Generating total return-level and equidensity contour plots.")
for(gauge in gauges){
  for(level in names(levels)){
    plot_density_contours(density_contours, level, gauge)
    plot_total_return_level(total_return_contours, level, gauge, ret_period_years)
  }
}


# Joint return level
joint_return_level <- function(phi_grid, gauge, level, ret_period_years, obs_year, grid_width = 4/(length(phi_grid) - 1), theta = 1){
  params <- gpd_estimates[[gauge]][[level]]$phi_grid$fitted
  thresh <- threshold_grid[[gauge]][[level]]
  phi_dens <- angular_density_grid[[gauge]]
  ret_period_obs <- ret_period_years * obs_year
  return(lapply(ret_period_obs, 
                FUN = function(ret){ 
                  data.frame(phi = pred_q, 
                             R = thresh + (params$scale/params$shape)*((grid_width*phi_dens*ret*(1 - levels[level]))^(params$shape) - 1)) |>
                             mutate(tm2 = R * gauge_vecs[[gauge]]$u, hs = R * gauge_vecs[[gauge]]$v)
                }))
}

plot_joint_return_level <- function(data, level, gauge){
  contour <- data[[gauge]][[level]]
  plot <- ggplot() + geom_point(data = ekofisk_tscaled, mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs, color = as.factor(.data[[ret_level_compare_colname(gauge, levels[level])]])), size = 0.1) +
    geom_path(data = contour, mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs, color = Period)) +
    labs(x = "tm2", y = "hs", title = paste("Joint Return-Level Contours", gauge_names[gauge], "gauge,", paste0(level, ", Theta = 1")), color = "Period")
  ggsave(paste0("plots/chapter 2/joint return level contours/", format(levels[level], nsmall = 2), "_", gauge, ".png"), plot=plot)
}

joint_retlevel_contours <- list(l1 = list(), l2 = list())
ret_levels <- c(10, 20, 40, 60, 100)
print("Estimating joint return-level contours")
for(gauge in gauges){
  for(level in names(levels)){
    contours <- data.table::rbindlist(joint_return_level(pred_q, gauge, level, ret_levels, obs_year), idcol = "Period")
    contours$Period <- as.character(ret_levels[contours$Period])
    joint_retlevel_contours[[gauge]][[level]] <- contours
  }
}

# Identifying observations outside of contour
# As this is an example, it expects grid of equidistant angles including both -2 and 2 for performance reasons
# Interpolation done under SPAR representation norm
data_contours_compare <- function(exc_data, contours, phi_grid, gauge, level){
  contours <- contours[[gauge]][[level]]
  exc_data <- exc_data |> mutate(phi_group = (phi+2)/4*(length(phi_grid)-1)+1, period = 0) |>
              mutate(lphi = phi_grid[floor(phi_group)], uphi = phi_grid[ceiling(phi_group)])
  periods <- unique(contours$Period)
  for(i in 1:length(periods)){
    contour <- contours[contours$Period == periods[i], c("phi", "R")]
    inds <- exc_data |> with(R > contour$R[floor(phi_group)] + (contour$R[ceiling(phi_group)] - contour$R[floor(phi_group)])*(phi-lphi))
    exc_data[inds, "period"] <- as.numeric(periods[i])
  }
  
  return(exc_data$period)
}

ret_level_compare_colname <- function(gauge, level){
  paste0("u", level, "_", gauge)
}

for(gauge in gauges){
  for(level in names(levels)){
    col <- ret_level_compare_colname(gauge, levels[level])
    exc_inds <- !is.na(excess[[gauge]][[level]])
    exc_obs <- ekofisk_tscaled[exc_inds, c(paste0("phi_", gauge), paste0("R_", gauge))]
    colnames(exc_obs) <- c("phi", "R") 
    ekofisk_tscaled[, col] <- 0
    ekofisk_tscaled[exc_inds, col] <- data_contours_compare(exc_obs, joint_retlevel_contours, pred_q, gauge, level)
  }
}

print("Plotting joint return-level contours")
for(gauge in gauges){
  for(level in names(levels)){
    plot_joint_return_level(joint_retlevel_contours, level, gauge)
  }
}

bin_probability <- list(l1 = 4/(length(pred_q) - 1)*(angular_density_grid$l1[1:(length(pred_q) -1)]+angular_density_grid$l1[2:length(pred_q)])/2,
                        l2 = 4/(length(pred_q) - 1)*(angular_density_grid$l2[1:(length(pred_q) -1)]+angular_density_grid$l2[2:length(pred_q)])/2)
ekofisk_tscaled <- ekofisk_tscaled |> mutate(bin_l1 = as.numeric(cut(phi_l1, breaks = pred_q[1:length(pred_q)], right=TRUE, labels=1:2000)),
                          bin_l2 = as.numeric(cut(phi_l2, breaks = pred_q[1:length(pred_q)], right=TRUE, labels=1:2000)))
# 100 * sum((ekofisk_tscaled |> group_by(bin_l1) |> filter(u0.8_l1 >= 100) |> summarize(count = n()) |> mutate(count = count * bin_probability$l1[bin_l1]))$count)


# Timeseries plots
timeseries <- data.frame(time = as.POSIXct(ekofisk_wave_surge$time, tz="UTC"), tm2 = ekofisk_wave_surge$tm2, hs = ekofisk_wave_surge$hs, year = format(as.Date(ekofisk_wave_surge$time), "%Y"), doy = ekofisk_wave_surge$doy)
years <- as.character(1980:1981)
timeseries <- data.table::rbindlist(list(timeseries |> filter(year %in% years) |> mutate(response = tm2, param = "tm2"),
                                         timeseries |> filter(year %in% years) |> mutate(response = hs, param = "hs")))

plot <- ggplot(timeseries, mapping=aes(x=doy, y=response)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(param~year) +
  labs(y = "Value", x="Day of Year", title = "Ekofisk Full Year Timeseries")
ggsave("plots/chapter 2/Timeseries/yearly_timeseries.png", plot)

timeseries_month <- data.frame(time = as.POSIXct(ekofisk_wave_surge$time, tz="UTC"), tm2 = ekofisk_wave_surge$tm2, hs = ekofisk_wave_surge$hs, month = format(as.Date(ekofisk_wave_surge$time), "%Y-%m"), doy = ekofisk_wave_surge$doy)
months <- paste(rep(1985,2), paste0("0", 1:2), sep = "-")
timeseries_month <- data.table::rbindlist(list(timeseries_month |> filter(month %in% months) |> mutate(response = tm2, param = "tm2"),
                                               timeseries_month |> filter(month %in% months) |> mutate(response = hs, param = "hs")))

plot <- ggplot(timeseries_month, mapping=aes(x=time, y=response)) +
  geom_line() +
  facet_grid(param~month, scales="free") +
  labs(title = "Ekofisk Month Timeseries", y = "Value", x="Time")
ggsave("plots/chapter 2/Timeseries/monthly_timeseries.png", plot)

peak <- which.max(ekofisk_tscaled$exc_80_l1)
plot <- ggplot() + 
    geom_point(data = ekofisk_tscaled, mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs), size = 0.1) +
    geom_segment(data=data.frame(ekofisk_tscaled[(peak-60):(peak+10),], t = (-60):(10)),
                 mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs, color=t,
                              xend=lead(tm2 * inv_s$tm2 + zero_t$tm2), yend = lead(hs * inv_s$hs + zero_t$hs)),
                 arrow = arrow(length = unit(0.2, "cm")), linewidth=1) +
    scale_color_gradient(low="green", high="red") +
    labs(title = "Ekofisk path around peak exceedance", x = "tm2", y = "hs", color = "Hours offset from peak")
ggsave("plots/chapter 2/Timeseries/peak_path_timeseries.png", plot)

inverse <- list(l1 = list(), l2 = list())
inverse$l1$x <- function(phi, r){
  return(ifelse(phi >= 0, 1 - phi, phi + 1) * r)
}
inverse$l1$y <- function(phi, r){
  return(ifelse(phi >= 0, 1 - abs(1 - phi), -1 + abs(phi + 1)) * r)
}

inverse$l2$x <- function(phi, r){
  return(ifelse(phi>=0, cos(pi * phi / 2), cos(-pi * phi / 2)) * r)
}
inverse$l2$y <- function(phi, r){
  return(ifelse(phi>=0, sqrt(1 - cos(pi * phi / 2)^2), -sqrt(1 - cos(-pi * phi / 2)^2)) * r)
}
transform_scale_inverse <- function(gauge){
  return(list(
    x = function(phi, r){inverse[[gauge]]$x(phi,r) * inv_s$tm2 + zero_t$tm2},
    y = function(phi, r){inverse[[gauge]]$y(phi,r) * inv_s$hs + zero_t$hs}
  ))
}

theme_radial_grid <- function(
    radial = element_line(),
    angular = element_line(),
    ...
) {
  theme(
    panel.grid.radial = radial,
    panel.grid.angular = angular,
    ...
  )
}

radial_grid_layer <- function(
    inverse,
    r_max,
    r_by = 1,
    n_spokes = 12,
    phi_grid_size = 400,
    theme = ggplot2::theme_get()
) {
  radial_el  <- ggplot2::calc_element("panel.grid.radial", theme)
  angular_el <- ggplot2::calc_element("panel.grid.angular", theme)
  
  # Fallback to standard ggplot grid styling
  major_el <- ggplot2::calc_element("panel.grid.major", theme)
  
  radial_colour <- radial_el$colour %||% major_el$colour %||% "grey80"
  radial_lwd    <- radial_el$linewidth %||% major_el$linewidth %||% 0.5
  radial_lty    <- radial_el$linetype %||% major_el$linetype %||% "solid"
  
  angular_colour <- angular_el$colour %||% major_el$colour %||% "grey80"
  angular_lwd    <- angular_el$linewidth %||% major_el$linewidth %||% 0.5
  angular_lty    <- angular_el$linetype %||% major_el$linetype %||% "solid"
  
  r_breaks <- seq(r_by, r_max, by = r_by)
  
  phi_grid <- seq(-2, 2, length.out = phi_grid_size)
  
  circles <- expand_grid(
    r = r_breaks,
    phi = phi_grid
  ) |>
    mutate(
      x = inverse$x(phi, r),
      y = inverse$y(phi, r),
      group = paste0("r_", r),
      type = "radial"
    )
  
  phi_breaks <- seq(-2, 2, length.out = n_spokes + 1)[-(n_spokes + 1)]
  
  spokes <- expand_grid(
    phi = phi_breaks,
    r = seq(0, r_max, length.out = 100)
  ) |>
    mutate(
      x = inverse$x(phi, r),
      y = inverse$y(phi, r),
      group = paste0("phi_", phi),
      type = "angular"
    )
  
  list(
    geom_path(
      data = circles,
      aes(x, y, group = group),
      inherit.aes = FALSE,
      colour = radial_colour,
      linewidth = radial_lwd,
      linetype = radial_lty
    ),
    geom_path(
      data = spokes,
      aes(x, y, group = group),
      inherit.aes = FALSE,
      colour = angular_colour,
      linewidth = angular_lwd,
      linetype = angular_lty
    )
  )
}

plot <- ggplot() + 
  radial_grid_layer(transform_scale_inverse("l1"), r_max = 20, r_by = 1, n_spokes = 100, phi_grid_size = 400) +
  geom_point(data = ekofisk_tscaled, mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs), size = 0.1) +
  geom_segment(data=data.frame(ekofisk_tscaled[(peak-60):(peak+10),], t = (-60):(10)),
               mapping = aes(x = tm2 * inv_s$tm2 + zero_t$tm2, y = hs * inv_s$hs + zero_t$hs, color=t,
                             xend=lead(tm2 * inv_s$tm2 + zero_t$tm2), yend = lead(hs * inv_s$hs + zero_t$hs)),
               arrow = arrow(length = unit(0.2, "cm")), linewidth=1) +
  scale_color_gradient(low="green", high="red") +
  labs(title = "Ekofisk path around peak exceedance", x = "tm2", y = "hs", color = "Hours offset from peak") +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 13)) +
  theme_get() +
  theme_radial_grid(
    radial  = element_line(colour = "grey80", linewidth = 0.3),
    angular = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid = element_blank(),
    validate = FALSE
  )
ggsave("plots/chapter 2/Timeseries/peak_path_timeseries_l1grid.png", plot)
