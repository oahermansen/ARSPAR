# Requires the ARSPAR package
devtools::load_all("ARSPAR")
if(!exists("ekofisk_wave_surge")){load("data/ekofisk_wave_surge.rda")}
library(dplyr)
library(ggplot2)
library(rsample)
library(purrr)
library(evgam)

theme_set(theme_gray(ink="blue"))
# Defining angular gauges
l1_radial <- function(M) {
  M <- as.matrix(M)
  abs(M[, 1]) + abs(M[, 2])
}
l1_angle <- function(M) {
  M <- as.matrix(M)
  s <- abs(M[, 1]) + abs(M[, 2])
  a <- sign(M[, 2])
  a[a == 0] <- 1
  out <- a * (1 - M[, 1] / s)
  out[s == 0] <- NA_real_
  out
}
l2_radial <- function(M){
  M <- as.matrix(M)
  sqrt(M[,1]^2+M[,2]^2)
}
l2_angle <- function(M){
  M <- as.matrix(M)
  s <- abs(M[, 1]) + abs(M[, 2])
  a <- sign(M[,2])
  a[a==0] <- 1
  a*(2/pi)*acos(M[,1]/s)
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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
colSd <- function(X){apply(X, MARGIN = 2, FUN = sd)}

if(!exists("spar")){
  spar <- spar_representation(X = ekofisk_wave_surge[, c("tm2", "hs")], time = ekofisk_wave_surge$time, observation_id = as.numeric(row.names(ekofisk_wave_surge)))
  
  spar <- spar |> spar_build_representation_transform(
    spar_step_mutate(.data = sweep(.data, 2, colMeans(.data), "-"), .name = "center_mean"),
    spar_step_mutate(.data = sweep(.data, 2, colSd(.data), "/"), .name = "scale_sd"),
    name = "center_scale",
    run = TRUE
  )
  # Auto-inverse not yet implemented in ARSPAR, this is a crude implementation not recommended for use.
  spar <- spar |> spar_set_transform_inverse(inverse = function(M){
    M = sweep(M, 2, colSd(spar$data$X_original), "*")
    sweep(M, 2, colMeans(spar$data$X_original), "+")
  }, "center_scale")
  
  angular_domain <- new_spar_angle_domain(type = "cyclical", lower = -2, upper = 2, closed_left = TRUE, closed_right = FALSE, label = "L1 Pseudo-angle")
  
  # Inverse functionality not complete in ARSPAR.
  angular_inverse <- function(M){
    x <- ifelse(M[, "phi"] >= 0, 1 - M[, "phi"], M[, "phi"] + 1) * M[,"R"]
    y <- ifelse(M[, "phi"] >= 0, 1 - abs(1 - M[, "phi"]), -1 + abs(M[, "phi"] + 1)) * M[, "R"]
    return(matrix(c(x,y), ncol = 2, dimnames = list(row.names(M), spar$schema$transformed_names)))
  }
  
  spar <- spar |> spar_build_angular_representation(
    radial_fun = l1_radial,
    angle_fun = l1_angle,
    domain = angular_domain,
    source = "transformed",
    inverse_fun = angular_inverse,
    name = "L1"
  )

}
phi_grid <- spar_angle_grid(n=2001, angular_domain)
tau_grid <- c(0.8, 0.85, 0.9, 0.95)
for(tau in tau_grid){
  name <- paste0("ald_tau_0_",tau*100)
  if(!(name %in% spar$threshold$registry$threshold_id)){
    spar <- spar |> spar_fit_threshold(
      method = "evgam_ald",
      name = name,
      tau = tau,
      phi_grid_n = 2001,
      k = 30,
      trace = 1,
      verbose = TRUE,
      set_active = TRUE,
      apply = TRUE,
      compute_excess = TRUE,
      storage = "compact"
    )
  }
}

spar <- spar |> spar_set_active_transform("center_scale", run = TRUE) |>
  spar_apply_angular_map("transformed", "center_scale") |>
  spar_apply_threshold("ald_tau_0_80") |> spar_update_excess_from_threshold()

# Computing excursion spans and clusters
spar <- spar |> 
  spar_decluster_excursions(gap_rule = 6)

# Obtain exceedance spans
exceedance_spans <- extract_exceedance_spans(spar$excursions$pointwise)

# Extract excursion path groups
path_groups <- spar |> 
  spar_build_excursion_path_group(
    gap_rule = 6,
    space = "transformed",
    keep_prev_next = TRUE,
    store = FALSE
  )
# Build the upper excursion paths
upper_paths <- spar |> spar_build_upper_excursion_paths(
  path_group = path_groups,
  gap_rule = gap_rule,
  space = "transformed",
  store = FALSE
)


gpd_formula <- list(as.formula("excess ~ s(phi, bs='cc', k=25)"), as.formula("excess ~ s(phi, bs='cc', k=35)"))
gpd_fit <- evgam(gpd_formula, tibble(excess = spar$excess$value, phi = spar$angular$phi) |> filter(excess > 0), family = "gpd", trace = 2)
params <- predict(gpd_fit, newdata = data.frame(phi=phi_grid), type = "response")

x <- lapply(upper_paths$paths, FUN = function(x){data.frame(lphi=x[["phi"]][1], uphi=x[["phi"]][length(x[["phi"]])])})
supports <- data.table::rbindlist(x)

inds <- supports$lphi <= supports$uphi
supports <- supports |> mutate(lphi = supports$lphi, uphi = ifelse(inds, supports$uphi, supports$uphi+4)) 
supports$max_exc <- vapply(upper_paths$paths, FUN = function(x){max(x[["excess"]])}, FUN.VALUE = numeric(1))
plot <- ggplot(data = supports, mapping = aes(x=lphi, y=uphi-lphi)) + geom_point(size = 0.1)
plot <- ggplot(data = supports, mapping = aes(x=uphi, y=uphi-lphi)) + geom_point(size = 0.1)


## Obtaining supports
## For each cluster cid obtain first, last
cluster_spans <- list()
for(i in 1:nrow(path_groups$spans$spans)){
  row <- path_groups$spans$spans[i,c("First", "Last")]
  spans <- exceedance_spans$spans |> filter(First >= row$First, Last <= row$Last, Last != First)
  if(nrow(spans) == 0) next
  angle_spans <- list()
  for(j in 1:nrow(spans)){
    range <- spans[j, "First"]:spans[j, "Last"]
    phis <- spar$angular$phi[range]
    angle_spans[[j]] <- spar_angle_range_from_transitions(phis, angular_domain)
  }
  cluster_spans[[i]] <- data.table::rbindlist(lapply(angle_spans, FUN=function(x){data.frame(lphi = x[1], uphi=x[2], cluster_id=i)}))
}
cluster_spans <- data.table::rbindlist(cluster_spans)

nahelp <- function(X){
  ifelse(is.na(X), 0, X)
}

plot <- ggplot(data = cluster_spans |> group_by(cluster_id) |> 
  summarize(lphi_min=min(lphi), uphi_sum=sum(uphi-lphi)), 
  mapping = aes(x=lphi_min, y=uphi_sum)) + geom_point(size = 0.1)
plot <- ggplot(data = cluster_spans |> group_by(cluster_id) |> arrange(lphi) |> mutate(gap = lead(lphi)-uphi) |> 
  summarize(lphi_min=min(lphi), uphi_sum=sum(uphi-lphi)+sum(nahelp(gap))), 
  mapping = aes(x=lphi_min, y=uphi_sum)) + geom_point(size = 0.1)


notunit_ids <- with(spar$excursions$clusters$spans, cluster_id[n_points > 1])
plot <- ggplot(data = cluster_spans |> group_by(cluster_id) |> filter(cluster_id %in% notunit_ids) |> 
         summarize(lphi_min=min(lphi), uphi_sum=sum(uphi-lphi)), 
       mapping = aes(x=lphi_min, y=uphi_sum)) + geom_point(size = 0.1)
plot <- ggplot(data = cluster_spans |> group_by(cluster_id) |> filter(cluster_id %in% notunit_ids) |> arrange(lphi) |> 
  mutate(gap = lead(lphi)-uphi) |> 
  summarize(lphi_min=min(lphi), uphi_sum=sum(uphi-lphi)+sum(nahelp(gap))), 
  mapping = aes(x=lphi_min, y=uphi_sum)) + geom_point(size = 0.1)

# Check amount of paths.
eval <- 0.67
in_range_table_eval <- function(phi, lphi, uphi, domain){
  # input list of clusters' lphi and uphi
  clusters <- list()
  for(i in 1:length(lphi)){
    spans <- list()
    for(j in 1:length(lphi[[i]])){
      ## We store a list indexed by tested spans
      spans[[j]] <- spar_angle_in_range(phi, c(lphi[[i]][j], uphi[[i]][j]), domain)
    }
    clusters[[i]] <- apply(matrix(rbind(unlist(spans)), ncol = length(spans)), MARGIN = 1, FUN = any)
  }
  return(clusters)
}
num_supports <- function(phi, cluster_spans, domain){
  ret <- (cluster_spans |> group_map(~ in_range_table_eval(phi, .x[, "lphi"], .x[, "uphi"], angular_domain)))
  return(apply(matrix(rbind(unlist(ret)), ncol = length(ret)), MARGIN = 1, FUN = sum))
}

years_data <- nrow(ekofisk_wave_surge)/(365.25*24)
sup_counts <- num_supports(phi_grid, cluster_spans |> group_by(cluster_id), angular_domain)
ang_yearly_event_rate <- sup_counts / years_data
lambda_f <- approxfun(x = phi_grid, y = ang_yearly_event_rate)

joint_return_level_rate <- function(phi_grid, ret_period_years, rate_years, params, threshold = "ald_tau_0_80"){
  thresh <- spar$threshold$functions[[threshold]](spar, phi = phi_grid)
  return(lapply(ret_period_years, 
                FUN = function(ret){ 
                  data.frame(phi = phi_grid, 
                             R = thresh + (params$scale/params$shape)*((ret*rate_years)^(params$shape) - 1))
                }))
}

ang_trans_inverse <- function(M){
  M <- spar$transform$inverse(spar$angular$inverse(M))
  colnames(M) <- c("tm2", "hs")
  return(M)
}
params <- list(pointwise = params)
ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), ang_yearly_event_rate, params$pointwise), idcol = "Period")
ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
plot <- ggplot() + geom_point(spar_data(spar, "original", "data.frame"), mapping = aes(x=tm2, y=hs), size = 1) + 
  geom_path(data = ret_level_level_rates, linewidth = 1,
            mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                          y = ang_trans_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) +
  labs(title = "Event-Rate Joint Return-level contours, iid GPD")
ggsave("temp_rate_joint_returnlevel_iid.png", plot)

plot <- ggplot() + 
  geom_point(data = data.frame(phi = phi_grid, rate = ang_yearly_event_rate), mapping = aes(x=phi, y=rate), color = "orange", size = 0.6) + 
  geom_function(fun = lambda_f, color = "red") + labs(title = "Annual event angular support rate")
ggsave("temp_annual_event_support_rate.png", plot)

binned_maxima <- data.table::rbindlist(lapply(path_groups$paths, function(x){
  data.frame(phi = x$phi, excess = x$excess, cluster_id = x$metadata$cluster_id)
})) |> filter(excess > 0 )|> group_by(cluster_id) |> mutate(phi_bin = cut(phi, breaks = phi_grid)) |> group_by(cluster_id, phi_bin) |> summarize(phi = phi[which.max(excess)], excess = max(excess))
binned_maximas <- binned_maxima |> ungroup() |>filter(excess > 0) |> select(phi, excess)
gpd_fits <- list(pointwise = gpd_fit)
gpd_fits$binned_maximas <- evgam(gpd_formula, data = binned_maximas, family = "gpd", trace = 2)
params$binned_maximas <- predict(gpd_fits$binned_maximas, newdata = data.frame(phi = phi_grid), type = "response")

ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), ang_yearly_event_rate, params$binned_maximas), idcol = "Period")
ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
plot <- ggplot() + geom_point(spar_data(spar, "original", "data.frame"), mapping = aes(x=tm2, y=hs), size = 1) + 
  geom_path(data = ret_level_level_rates, linewidth = 1,
            mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                          y = ang_trans_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
  geom_path(data = data.frame(spar$threshold$estimators$ald_tau_0_80$grid), 
            mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=u))[,"tm2"],
                          y = ang_trans_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
  labs(title = "Event-Rate Joint Return-level contours, binned maximas GPD")
ggsave("temp_rate_joint_returnlevel_binned.png", plot)

plot <- ggplot() + geom_point(spar_data(spar, "angular", "data.frame"), mapping = aes(x=phi, y = R)) +
  geom_point(data = data.frame(phi = phi_grid, R = ang_yearly_event_rate), mapping = aes(x = phi, y = R), color = "yellow")+ 
  geom_path(data = data.frame(spar$threshold$estimators$ald_tau_0_80$grid), 
            mapping = aes(x = phi, y = u), color = "red")

longer_inds <- path_groups$spans$spans |> filter(n_points > 5) |> select(cluster_id)
upper_sample_path <- lapply(upper_paths$paths[longer_inds$cluster_id], FUN = function(x){data.frame(phi = x$phi, R = x$R, excess = x$excess)})
upper_sample_path <- data.table::rbindlist(upper_sample_path) |> filter(excess > 0)
gpd_fits$upper_sample <- evgam(gpd_formula, data = upper_sample_path, family = "gpd", trace = 2)
params$upper_sample <- predict(gpd_fits$upper_sample, newdata = data.frame(phi = phi_grid), type = "response")
ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), ang_yearly_event_rate, params$upper_sample), idcol = "Period")
ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
plot <- ggplot() + geom_point(spar_data(spar, "original", "data.frame"), mapping = aes(x=tm2, y=hs), size = 1) + 
  geom_path(data = ret_level_level_rates, linewidth = 1,
            mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                          y = ang_trans_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
  geom_path(data = data.frame(spar$threshold$estimators$ald_tau_0_80$grid), 
            mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=u))[,"tm2"],
                          y = ang_trans_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
  labs(title = "Event-Rate Joint Return-level contours, upper sampled path GPD")
ggsave("temp_rate_joint_returnlevel_uppersampledpath.png", plot)

density_contours <- list(l1 = list(), l2 = list())
total_return_contours <- list(l1 = list(), l2 = list())
density_levels <- 10^-(3:8)
ret_period_years <- 10
obs_year <- 365.25*24



## Log scale
inv_mean <- c(log(mean(spar$data$X_original[,1])), log(mean(spar$data$X_original[,2])))
loginverse <- function(M){
  M[,1] <- exp(M[,1]+inv_mean[1])
  M[,2] <- exp(M[,2]+inv_mean[2])
  return(M)
}
if(!("log_transform" %in% names(spar$transform$chains))){
  spar <- spar |> spar_build_representation_transform(
    spar_step_mutate(tm2 = log(tm2), hs = log(hs), .name = "log_transform"),
    spar_step_mutate(tm2 = tm2 - log(mean(.tm2)), hs = hs - log(mean(.hs)), .name = "center"),
    name = "log_transform",
    run = TRUE
  ) |> spar_set_transform_inverse(inverse = loginverse)
}
spar <- spar |> spar_set_active_transform("log_transform", run = TRUE) |> 
  spar_set_active_angular_map("L1", source = "transformed", transform_id = "log_transform") |>
  spar_apply_angular_map("transformed", "log_transform")

k_grid <- c(40, 40, 45, 45)
for(i in 1:length(tau_grid)){
  tau <- tau_grid[i]
  k <- k_grid[i]
  name <- paste0("log_ald_tau_0_",tau*100)
  if(!(name %in% spar$threshold$registry$threshold_id))
  {
    spar <- spar |> spar_fit_threshold(
      method = "evgam_ald",
      name = name,
      tau = tau,
      k = k,
      trace = 2,
      verbose = TRUE,
      set_active = TRUE,
      apply = TRUE,
      compute_excess = TRUE,
      storage = "compact"
    )
  }
}

inv_mean1p <- c(log1p(mean(spar$data$X_original[,1])), log1p(mean(spar$data$X_original[,2])))
log1pinverse <- function(M){
  M[,1] <- exp(M[,1]+inv_mean1p[1])-1
  M[,2] <- exp(M[,2]+inv_mean1p[2])-1
  return(M)
}

if(!("log1p_transform" %in% names(spar$transform$chains))){
  spar <- spar |> spar_build_representation_transform(
    spar_step_mutate(tm2 = log1p(tm2), hs = log1p(hs), .name = "log_transform"),
    spar_step_mutate(tm2 = tm2 - log1p(mean(.tm2)), hs = hs - log1p(mean(.hs)), .name = "center"),
    name = "log1p_transform",
    run = TRUE
  ) |> spar_set_transform_inverse(inverse = log1pinverse)
}
spar <- spar |> spar_set_active_transform("log1p_transform", run = TRUE) |> 
  spar_set_active_angular_map("L1", source = "transformed", transform_id = "log1p_transform") |>
  spar_apply_angular_map("transformed", "log1p_transform")

k_grid <- c(40, 40, 45, 45)
for(i in 1:length(tau_grid)){
  tau <- tau_grid[i]
  k <- k_grid[i]
  name <- paste0("log1p_ald_tau_0_",tau*100)
  if(!(name %in% spar$threshold$registry$threshold_id))
  {
    spar <- spar |> spar_fit_threshold(
      method = "evgam_ald",
      name = name,
      tau = tau,
      k = k,
      trace = 2,
      verbose = TRUE,
      set_active = TRUE,
      apply = TRUE,
      compute_excess = TRUE,
      storage = "compact"
    )
  }
}

thresh_reg <- spar$threshold$registry
ang_trans_inverse <- function(M){
  M <- spar$transform$inverse((spar$angular$inverse(M)))
  colnames(M) <- c("tm2", "hs")
  return(M)
}

qq_plot_regions <- list(list(name = "All Exceedances", min = -2, max = 2),
                        list(name = "First Quadrant", min = 0, max = 1),
                        list(name = "Second Quadrant", min = 1, max = 2),
                        list(name = "Third Quadrant", min = -2, max = -1),
                        list(name = "Fourth Quadrant", min = -1, max = 0),
                        list(name = "Growth Region", min = 0.7, max = 1.5),
                        list(name = "Growth Asymptotic Region", min = 0.675, max = 0.775))

plot_qq_regions <- function(gpd_fit, data, level, transform, fit_type){
  for(region in qq_plot_regions){
    plot.new()
    png(paste(paste0("plots/chapter 6/",transform,"/qqplots/"), level, region$name, fit_type, "qqplot.png"))
    ang <- data$phi
    indexset <- ang >= region$min & ang <= region$max
    predict(gpd_fit, newdata = data[indexset,], type = "qqplot")
    title(main = paste(region$name, transform, level), line = 0.5, cex.main = 0.75)
    dev.off()
  }
}

gpd_formula <- list(as.formula("excess ~ s(phi, bs='cc', k=35)"), as.formula("excess ~ s(phi, bs='cc', k=20)"))
cluster_counts <- list()
gpd_params <- list()
event_rates <- list()
for(transform in unique(thresh_reg$transform_id)){
  if(length(gpd_params) < 1){
    gpd_params[[transform]] <- list(pointwise = list(), binned_maximas = list(), upper_sample_envelope = list())
  }
  cluster_counts[[transform]] <- list()
  thresh_names <- paste("Tau =", tau_grid)
  names(thresh_names) <- (thresh_reg |> filter(transform_id == transform))$threshold_id
  angular_rates <- list()
  names(tau_grid) <- names(thresh_names)
  params <- gpd_params[[transform]]
  for(threshold in names(thresh_names)){
    if(!is.null(params$upper_sample_envelope[[thresh_names[[threshold]]]])){
      next
    }
    print(paste("Began", transform, thresh_names[threshold]))
    spar <- spar |> spar_set_active_transform(transform, run = TRUE) |> 
      spar_set_active_angular_map("L1", source = "transformed", transform_id = transform) |>
      spar_apply_angular_map("transformed", transform) |>
      spar_apply_threshold(threshold)
    
    # Computing excursion spans and clusters
    print("Declustering")
    spar <- spar |> 
      spar_decluster_excursions(gap_rule = 6)
    
    # Obtain exceedance spans
    print("Extracting exceedance spans")
    exceedance_spans <- extract_exceedance_spans(spar$excursions$pointwise)
    
    # Extract excursion path groups
    print("Building path groups")
    path_groups <- spar |> 
      spar_build_excursion_path_group(
        gap_rule = 6,
        space = "transformed",
        keep_prev_next = TRUE,
        store = FALSE
      )
    # Build the upper excursion paths
    print("Building sample envelope")
    upper_paths <- spar |> spar_build_upper_excursion_paths(
      path_group = path_groups,
      gap_rule = gap_rule,
      space = "transformed",
      store = FALSE
    )
    cluster_counts[[transform]][[threshold]] <- nrow(path_groups$paths)
    
    print("Building cluster spans for rate estimation")
    cluster_spans <- list()
    for(i in 1:nrow(path_groups$spans$spans)){
      row <- path_groups$spans$spans[i,c("First", "Last")]
      spans <- exceedance_spans$spans |> filter(First >= row$First, Last <= row$Last, First != Last)
      if(nrow(spans) == 0) next
      angle_spans <- list()
      for(j in 1:nrow(spans)){
        range <- spans[j, "First"]:spans[j, "Last"]
        phis <- spar$angular$phi[range]
        angle_spans[[j]] <- spar_angle_range_from_transitions(phis, angular_domain)
      }
      cluster_spans[[i]] <- data.table::rbindlist(lapply(angle_spans, FUN=function(x){data.frame(lphi = x[1], uphi=x[2], cluster_id=i)}))
    }
    cluster_spans <- data.table::rbindlist(cluster_spans)
    
    in_range_table_eval <- function(phi, lphi, uphi, domain){
      # input list of clusters' lphi and uphi
      clusters <- list()
      for(i in 1:length(lphi)){
        spans <- list()
        for(j in 1:length(lphi[[i]])){
          ## We store a list indexed by tested spans
          spans[[j]] <- spar_angle_in_range(phi, c(lphi[[i]][j], uphi[[i]][j]), domain)
        }
        clusters[[i]] <- apply(matrix(rbind(unlist(spans)), ncol = length(spans)), MARGIN = 1, FUN = any)
      }
      return(clusters)
    }
    num_supports <- function(phi, cluster_spans, domain){
      ret <- (cluster_spans |> group_map(~ in_range_table_eval(phi, .x[, "lphi"], .x[, "uphi"], angular_domain)))
      return(apply(matrix(rbind(unlist(ret)), ncol = length(ret)), MARGIN = 1, FUN = sum))
    }
    
    years_data <- nrow(ekofisk_wave_surge)/(365.25*24)
    sup_counts <- num_supports(phi_grid, cluster_spans |> group_by(cluster_id), angular_domain)
    angular_rates[[thresh_names[threshold]]] <- tibble(phi = phi_grid, rate = sup_counts / years_data)
    lambda_f <- approxfun(x = phi_grid, y = angular_rates[[thresh_names[threshold]]]$rate)
    
    gpd_fits <- list()
    gpd_fits$pointwise <- evgam(gpd_formula, data.frame(excess = spar$excess$value, phi = spar$angular$phi) |> filter(excess > 0), family = "gpd", trace = 2)
    params$pointwise[[thresh_names[threshold]]] <- predict(gpd_fits$pointwise, newdata = data.frame(phi=phi_grid), type = "response")
    
    print("Computing binned maxima")
    binned_maxima <- data.table::rbindlist(lapply(path_groups$paths, function(x){
      tibble(phi = x$phi, excess = x$excess, cluster_id = x$metadata$cluster_id)
    })) |> filter(excess > 0 )|> group_by(cluster_id) |> mutate(phi_bin = cut(phi, breaks = phi_grid, right = FALSE)) |> group_by(cluster_id, phi_bin) |> summarize(phi = phi[which.max(excess)], excess = max(excess))
    gpd_fits$binned_maximas <- evgam(gpd_formula, data = binned_maxima, family = "gpd", trace = 2)
    params$binned_maximas[[thresh_names[threshold]]] <- predict(gpd_fits$binned_maximas, newdata = data.frame(phi = phi_grid), type = "response")
    
    upper_sample_path <- lapply(upper_paths$paths, FUN = function(x){data.frame(phi = x$phi, R = x$R, excess = x$excess)})
    upper_sample_path <- data.table::rbindlist(upper_sample_path) |> filter(excess > 0)
    gpd_fits$upper_sample <- evgam(gpd_formula, data = upper_sample_path, family = "gpd", trace = 2)
    params$upper_sample_envelope[[thresh_names[[threshold]]]] <- predict(gpd_fits$upper_sample, newdata = data.frame(phi = phi_grid), type = "response")
    
    ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), lambda_f(phi_grid), params$pointwise[[thresh_names[threshold]]], threshold = threshold), idcol = "Period")
    ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
    plot <- ggplot() + geom_point(spar_data(spar, "original", "data.frame"), mapping = aes(x=tm2, y=hs), size = 0.1) + 
      geom_path(data = ret_level_level_rates, linewidth = 1,
                mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                              y = ang_trans_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
      geom_path(data = data.frame(spar$threshold$estimators[[threshold]]$grid), 
                mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=u))[,"tm2"],
                              y = ang_trans_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
      labs(title = "Event-Rate Joint Return-level contours, iid GPD", subtitle = paste0("Using  ", transform," transform, L1 gauge and ", thresh_names[threshold],"."))
    ggsave(paste0("plots/chapter 6/", transform, "/", tau_grid[threshold]*100, "_event-rate_contours_all_exceedances.png"), plot)
    
    ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), lambda_f(phi_grid), params$binned_maximas[[thresh_names[threshold]]], threshold = threshold), idcol = "Period")
    ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
    plot <- ggplot() + geom_point(spar_data(spar, "original", "data.frame"), mapping = aes(x=tm2, y=hs), size = 0.1) + 
      geom_path(data = ret_level_level_rates, linewidth = 1,
                mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                              y = ang_trans_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
      geom_path(data = data.frame(spar$threshold$estimators[[threshold]]$grid), 
                mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=u))[,"tm2"],
                              y = ang_trans_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
      labs(title = "Event-Rate Joint Return-level contours, binned maximas GPD", subtitle = paste0("Using  ", transform," transform, L1 gauge and ", thresh_names[threshold],"."))
    ggsave(paste0("plots/chapter 6/", transform, "/", tau_grid[threshold]*100, "_event-rate_contours_binned_maximas.png"), plot)
    
    ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), lambda_f(phi_grid), params$upper_sample_envelope[[thresh_names[threshold]]], threshold = threshold), idcol = "Period")
    ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
    plot <- ggplot() + geom_point(spar_data(spar, "original", "data.frame"), mapping = aes(x=tm2, y=hs), size = 0.1) +
      geom_path(data = ret_level_level_rates, linewidth = 1,
                mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                              y = ang_trans_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
      geom_path(data = data.frame(spar$threshold$estimators[[threshold]]$grid), 
                mapping = aes(x = ang_trans_inverse(data.frame(phi=phi, R=u))[,"tm2"], 
                              y = ang_trans_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
      labs(title = "Event-Rate Joint Return-level contours, upper sampled path GPD", subtitle = paste0("Using  ", transform," transform, L1 gauge and ", thresh_names[threshold],"."))
    ggsave(paste0("plots/chapter 6/", transform, "/", tau_grid[threshold]*100, "_event-rate_contours_upper_sample_path.png"), plot)
    
    if(transform == "log_transform" || transform == "log1p_transform"){
      ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), lambda_f(phi_grid), params$pointwise[[thresh_names[threshold]]], threshold = threshold), idcol = "Period")
      ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
      plot <- ggplot() + geom_point(spar_data(spar, "transformed", "data.frame"), mapping = aes(x=tm2, y=hs), size = 0.1) + 
        geom_path(data = ret_level_level_rates, linewidth = 1,
                  mapping = aes(x = angular_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                                y = angular_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
        geom_path(data = data.frame(spar$threshold$estimators[[threshold]]$grid), 
                  mapping = aes(x = angular_inverse(data.frame(phi=phi, R=u))[,"tm2"],
                                y = angular_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
        labs(title = "Logscale Event-Rate Joint Return-level contours, all exceedances GPD", subtitle = paste0("Using  ", transform," transform, L1 gauge and ", thresh_names[threshold],"."))
      ggsave(paste0("plots/chapter 6/", transform, "/log/", tau_grid[threshold]*100, "_event-rate_contours_all_exceedances.png"), plot)
      
      ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), lambda_f(phi_grid), params$binned_maximas[[thresh_names[threshold]]], threshold = threshold), idcol = "Period")
      ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
      plot <- ggplot() + geom_point(spar_data(spar, "transformed", "data.frame"), mapping = aes(x=tm2, y=hs), size = 0.1) + 
        geom_path(data = ret_level_level_rates, linewidth = 1,
                  mapping = aes(x = angular_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                                y = angular_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
        geom_path(data = data.frame(spar$threshold$estimators[[threshold]]$grid), 
                  mapping = aes(x = angular_inverse(data.frame(phi=phi, R=u))[,"tm2"],
                                y = angular_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
        labs(title = "Logscale Event-Rate Joint Return-level contours, binned maximas GPD", subtitle = paste0("Using  ", transform," transform, L1 gauge and ", thresh_names[threshold],"."))
      ggsave(paste0("plots/chapter 6/", transform, "/log/", tau_grid[threshold]*100, "_event-rate_contours_binned_maximas.png"), plot)

      ret_level_level_rates <- data.table::rbindlist(joint_return_level_rate(phi_grid, c(10,20,40,60,100), lambda_f(phi_grid), params$upper_sample_envelope[[thresh_names[threshold]]], threshold = threshold), idcol = "Period")
      ret_level_level_rates$Period <- as.factor(c(10,20,40,60,100)[ret_level_level_rates$Period])
      plot <- ggplot() + geom_point(spar_data(spar, "transformed", "data.frame"), mapping = aes(x=tm2, y=hs), size = 0.1) + 
        geom_path(data = ret_level_level_rates, linewidth = 1,
                  mapping = aes(x = angular_inverse(data.frame(phi=phi, R=R))[,"tm2"], 
                                y = angular_inverse(data.frame(phi=phi, R=R))[,"hs"], color=Period)) + 
        geom_path(data = data.frame(spar$threshold$estimators[[threshold]]$grid), 
                  mapping = aes(x = angular_inverse(data.frame(phi=phi, R=u))[,"tm2"], 
                                y = angular_inverse(data.frame(phi=phi, R=u))[,"hs"]), color = "red") +
        labs(title = "Logscale Event-Rate Joint Return-level contours, upper sampled path GPD", subtitle = paste0("Using  ", transform," transform, L1 gauge and ", thresh_names[threshold],"."))
      ggsave(paste0("plots/chapter 6/", transform, "/log/", tau_grid[threshold]*100, "_event-rate_contours_upper_sample_path.png"), plot)
    }
    
    plot.new()
    png(paste(paste0("plots/chapter 6/", transform, "/qqplots/", tau_grid[threshold]*100, "_all_exc_qqplot_GPD_L1.png")))
    predict(gpd_fits$pointwise, newdata = data.frame(phi=spar$angular$phi, excess=spar$excess$value) |> filter(excess > 0), type = "qqplot")
    title(main = paste("All Exceedances QQ-Plot", transform, thresh_names[threshold]), line = 0.5, cex.main = 0.75)
    dev.off()
    plot_qq_regions(gpd_fits$pointwise, data = data.frame(phi=spar$angular$phi, excess=spar$excess$value) |> filter(excess > 0), 
                    level = thresh_names[threshold], transform = transform, fit_type = "All Exceedances")
    plot_qq_regions(gpd_fits$binned_maximas, data = as.data.frame(binned_maxima[, c("phi", "excess")]), 
                    level = thresh_names[threshold], transform = transform, fit_type = "Binned Maxima")
    plot_qq_regions(gpd_fits$upper_sample, data = upper_sample_path, 
                    level = thresh_names[threshold], transform = transform, fit_type = "Sample Envelope")
    plot.new()
    png(paste(paste0("plots/chapter 6/", transform, "/qqplots/", tau_grid[threshold]*100, "_bin_qqplot_GPD_L1.png")))
    predict(gpd_fits$binned_maximas, newdata = as.data.frame(binned_maxima[, c("phi", "excess")]), type = "qqplot")
    title(main = paste("Binned Maxima QQ-Plot", transform, thresh_names[threshold]), line = 0.5, cex.main = 0.75)
    dev.off()
    plot.new()
    png(paste(paste0("plots/chapter 6/", transform, "/qqplots/", tau_grid[threshold]*100, "_sample_env_qqplot_GPD_L1.png")))
    predict(gpd_fits$upper_sample, newdata = upper_sample_path, type = "qqplot")
    title(main = paste("Sample Envelope QQ-Plot", transform, thresh_names[threshold]), line = 0.5, cex.main = 0.75)
    dev.off()
    rm(gpd_fits)
  }
  if(length(angular_rates) < 1){
    next
  }
  angular_rates <- data.table::rbindlist(angular_rates, idcol = "Level")
  angular_rates <- angular_rates
  plot <- ggplot(data = angular_rates) + 
    geom_line(mapping = aes(x=phi, y=rate, color=Level)) + labs(title = "Annual event-rate per angle.", x=expression(phi), y=expression(lambda*" Annual Rate"), subtitle = paste0("For ", transform, " transform using L1 Gauge.")) 
  ggsave(paste0("plots/chapter 6/", transform, "/annual_event-rates.png"), plot)
  param_df <- data.frame()
  for(threshold in thresh_names){
    par <- data.table::rbindlist(list(as.data.frame(params$pointwise[[threshold]]) |> mutate(phi = phi_grid, fit = "Pointwise"),
                                      as.data.frame(params$binned_maximas[[threshold]]) |> mutate(phi = phi_grid, fit = "Binned"),
                                      as.data.frame(params$upper_sample_envelope[[threshold]]) |> mutate(phi = phi_grid, fit = "Sample Envelope")))
    param_df <- data.table::rbindlist(c(list(param_df), list(par |> mutate(Threshold = threshold))))
  }
  plot <- ggplot(data = param_df) + geom_line(mapping = aes(phi, scale, color=Threshold)) +
    facet_grid(cols=vars(fit)) +
    labs(title = expression("Estimated GP Scale Parameters "*sigma*"."), 
         x=expression("Pseudo-Angle "*phi), y=expression("Scale "*sigma), subtitle = paste0("For ", transform, " transform using L1 Gauge."))
  ggsave(paste0("plots/chapter 6/", transform, "/scale_parameters.png"), plot, width = 11, height = 4)
  plot <- ggplot(data = param_df) + geom_line(mapping = aes(phi, shape, color=Threshold)) +
    facet_grid(cols=vars(fit)) +
    labs(title = expression("Estimated GP Shape Parameters "*xi*"."), 
         x=expression("Pseudo-Angle "*phi), y=expression("Shape "*xi), subtitle = paste0("For ", transform, " transform using L1 Gauge."))
  ggsave(paste0("plots/chapter 6/", transform, "/shape_parameters.png"), plot, width = 11, height = 4)
  
  gpd_params[[transform]] <- params
  event_rates[[transform]] <- angular_rates
}

# Bootstrap estimation for log1p, Tau = 0.85
print("Starting bootstrap estimation for log1p transform at threshold level Tau = 0.85.")
spar <- spar |> spar_set_active_transform("log1p_transform", run = TRUE) |> 
  spar_set_active_angular_map("L1", source = "transformed", transform_id = "log1p_transform") |>
  spar_apply_angular_map("transformed", "log1p_transform") |>
  spar_apply_threshold("logp1_ald_tau_0_85")

# Computing excursion spans and clusters
print("Declustering")
spar <- spar |> 
  spar_decluster_excursions(gap_rule = 6)

# Obtain exceedance spans
print("Extracting exceedance spans")
exceedance_spans <- extract_exceedance_spans(spar$excursions$pointwise)

# Extract excursion path groups
print("Building path groups")
path_groups <- spar |> 
  spar_build_excursion_path_group(
    gap_rule = 6,
    space = "transformed",
    keep_prev_next = TRUE,
    store = FALSE
  )
# Build the upper excursion paths
print("Building sample envelope")
upper_paths <- spar |> spar_build_upper_excursion_paths(
  path_group = path_groups,
  gap_rule = gap_rule,
  space = "transformed",
  store = FALSE
)

print("Building cluster spans for rate estimation")
cluster_spans <- list()
for(i in 1:nrow(path_groups$spans$spans)){
  row <- path_groups$spans$spans[i,c("First", "Last")]
  spans <- exceedance_spans$spans |> filter(First >= row$First, Last <= row$Last, First != Last)
  if(nrow(spans) == 0) next
  angle_spans <- list()
  for(j in 1:nrow(spans)){
    range <- spans[j, "First"]:spans[j, "Last"]
    phis <- spar$angular$phi[range]
    angle_spans[[j]] <- spar_angle_range_from_transitions(phis, angular_domain)
  }
  cluster_spans[[i]] <- data.table::rbindlist(lapply(angle_spans, FUN=function(x){data.frame(lphi = x[1], uphi=x[2], cluster_id=i)}))
}
cluster_spans <- data.table::rbindlist(cluster_spans)

in_range_table_eval <- function(phi, lphi, uphi, domain){
  # input list of clusters' lphi and uphi
  clusters <- list()
  for(i in 1:length(lphi)){
    spans <- list()
    for(j in 1:length(lphi[[i]])){
      ## We store a list indexed by tested spans
      spans[[j]] <- spar_angle_in_range(phi, c(lphi[[i]][j], uphi[[i]][j]), domain)
    }
    clusters[[i]] <- apply(matrix(rbind(unlist(spans)), ncol = length(spans)), MARGIN = 1, FUN = any)
  }
  return(clusters)
}
num_supports <- function(phi, cluster_spans, domain){
  ret <- (cluster_spans |> group_map(~ in_range_table_eval(phi, .x[, "lphi"], .x[, "uphi"], angular_domain)))
  return(apply(matrix(rbind(unlist(ret)), ncol = length(ret)), MARGIN = 1, FUN = sum))
}
num_supports_df <- function(phi, sup_data, c_inds, domain){
  return(apply(matrix(rbind(unlist(sup_data[c_inds,])), ncol = length(ret)), MARGIN = 1, FUN = sum))
}
n_resample <- 200
resample_years <- 25
years_data <- nrow(ekofisk_wave_surge)/(365.25*24)
set.seed(123)
event_counts <- rpois(n_resample, resample_years * length(unique(cluster_spans$cluster_id))/years_data)
bstraps <- bootstraps(data = data.frame(ind=unique(cluster_spans$cluster_id)), times = n_resample)
print(paste("Starting bootstrap procedure with n_resample =", n_resample, "T =", resample_years, "years."))

results <- map2(bstraps$splits, event_counts, \(split, count){
  inds <- as.data.frame(split)$ind[1:count]
  cids <- 1:count
  c_spans <- data.table::rbindlist(map2(inds, cids, \(x, y){cluster_spans |> filter(cluster_id == x) |> mutate(cid=y)}))
  
  sup_counts <- num_supports(phi_grid, c_spans |> group_by(cid), angular_domain)
  rate <- tibble(phi = phi_grid, rate = sup_counts / resample_years)
  lambda_f <- approxfun(x = phi_grid, y = rate$rate)

  upper_sample_path <- lapply(upper_paths$paths[inds], FUN = function(x){data.frame(phi = x$phi, R = x$R, excess = x$excess)})
  upper_sample_path <- data.table::rbindlist(upper_sample_path) |> filter(excess > 0)
  gpd_fit <- evgam(gpd_formula, data = upper_sample_path, family = "gpd")
  params <- predict(gpd_fit, newdata = data.frame(phi = phi_grid), type = "response")
  rates <- lambda_f(phi_grid)
  tibble(params=params, rates=rates)
}, .progress = TRUE)

result_df <- data.table::rbindlist(lapply(results, FUN=\(x){tibble(phi = phi_grid, scale = x$params$scale, shape = x$params$shape, rate = x$rates)}), idcol = "id")

saveRDS(results, file = "bootstrap_results_log1p.rds")

ret_levels_all_df <- data.table::rbindlist(result_df |> group_by(id) |> group_map(\(fit,...){
  data.table::rbindlist(
    joint_return_level_rate(phi_grid, c(10,25,50,100,500), 
                            fit$rate, fit[, c("scale", "shape")], 
                            threshold = "logp1_ald_tau_0_85"), idcol = "Period")
}), idcol = "id")


summary_df <- result_df |> group_by(phi) |> 
  summarise(scale_l = quantile(scale, probs=c(0.025)), scale_u = quantile(scale, probs=c(0.975)), 
            shape_l = quantile(shape, probs=c(0.025)), shape_u = quantile(shape, probs=c(0.975)),
            rate_l  = quantile(rate,  probs=c(0.025)), rate_u  = quantile(rate,  probs=c(0.975)),
            scale=mean(scale), shape = mean(shape), rate = mean(rate))

ret_levels_df <- ret_levels_all_df |> 
  summarise(R_l = quantile(R, probs=c(0.1)), R_u = quantile(R, probs=c(0.9)),
            R   = mean(R), .by = c(Period, phi))
ret_levels_df$Period <- as.factor(c(10,25,50,100,500))[ret_levels_df$Period]


ret_levels_df[, c("tm2", "hs")] <- ang_trans_inverse(ret_levels_df |> select(phi, R))
ret_levels_df[, c("tm2_l", "hs_l")] <- ang_trans_inverse(ret_levels_df |> select(phi, R=R_l))
ret_levels_df[, c("tm2_u", "hs_u")] <- ang_trans_inverse(ret_levels_df |> select(phi, R=R_u))

ret_levels_all_df[, c("tm2", "hs")] <- as.data.frame(ang_trans_inverse(as.matrix(ret_levels_all_df |> select(phi, R))))

ret_levels_polygon <- data.table::rbindlist(
  lapply(unique(ret_levels_df$Period), FUN=\(x){
    data.table::rbindlist(
      list(ret_levels_df[ret_levels_df$Period == x, ] |> select(tm2=tm2_l, hs=hs_l, Period),
           (ret_levels_df[ret_levels_df$Period == x, ] |> select(tm2=tm2_u, hs=hs_u, Period))[sum(ret_levels_df$Period==x):1,])
      )}))

plot <- ggplot(data = summary_df) + 
  geom_line(mapping=aes(phi, scale)) + 
  geom_line(mapping=aes(phi, scale_l), linetype = "dotdash") + 
  geom_line(mapping=aes(phi, scale_u), linetype = "dotdash") +
  labs(
    x = expression("Pseudo-angle "*phi),
    y = expression("Scale "*sigma),
    title = expression("Bootstrap estimated GP scale parameter "*sigma*" uncertainty, "*T==25*" years."),
    subtitle = expression("Quantiles "*p %in% paste("[", 0.025, ", ", 0.975, "]")*" for "*log1p*" transformation with "*L[1]*" gauge and threshold level "*tau==0.85*".")
  )
ggsave(paste0("plots/chapter 6/bootstraps/log1p_25year_scale_parameter.png"), plot, width = 7, height = 7)

plot <- ggplot(data = summary_df) + 
  geom_line(mapping=aes(phi, shape)) + 
  geom_line(mapping=aes(phi, shape_l), linetype = "dotdash") + 
  geom_line(mapping=aes(phi, shape_u), linetype = "dotdash") +
  labs(
    x = expression("Pseudo-angle "*phi),
    y = expression("Shape "*xi),
    title = expression("Bootstrap estimated GP shape parameter "*xi*" uncertainty, "*T==25*" years."),
    subtitle = expression("Quantiles "*p %in% paste("[", 0.025, ", ", 0.975, "]")*" for "*log1p*" transformation with "*L[1]*" gauge and threshold level "*tau==0.85*".")
  )
ggsave(paste0("plots/chapter 6/bootstraps/log1p_25year_shape_parameter.png"), plot, width = 7, height = 7)

plot <- ggplot(data=ret_levels_df, mapping=aes(x=tm2, y=hs)) +
  geom_point(data=spar_data(spar, space = "original", format = "tibble"),
             mapping=aes(x=tm2, y=hs), shape = ".") +
  geom_polygon(data=ret_levels_polygon, mapping = aes(x=tm2, y=hs, fill = Period), alpha = 0.2) +
  geom_path(linewidth = 1,
            mapping = aes(x = tm2, y = hs,
                          color = Period)) +
  labs(
    color = "Return period (years)",
    fill = "Return period (years)",
    x = expression(T[m02] ~ "(s)"),
    y = expression(Hs ~ "(m)"),
    title = expression("Bootstrap estimated return-level contour uncertainty, "*T==25*" years."),
    subtitle = expression("Quantiles "*p %in% paste("[", 0.1, ", ", 0.9, "]")*" for "*log1p*" transformation with "*L[1]*" gauge and threshold level "*tau==0.85*".")
  )
ggsave(paste0("plots/chapter 6/bootstraps/log1p_25year_return_levels.png"), plot, width = 8, height = 5.5)

ggplot(data=ret_levels_df, mapping=aes(x=phi, y=R)) +
  geom_point(data=spar_data(spar, space = "angular", format = "tibble"),
             mapping=aes(x=phi, y=R), shape = ".") +
  geom_path(linewidth = 1,
            mapping = aes(x = phi, y = R,
                          color = Period))


n_resample <- 200
resample_years <- 10
set.seed(123)
event_counts <- rpois(n_resample, resample_years * length(unique(cluster_spans$cluster_id))/years_data)
bstraps <- bootstraps(data = data.frame(ind=unique(cluster_spans$cluster_id)), times = n_resample)
print(paste("Starting bootstrap procedure with n_resample =", n_resample, "T =", resample_years, "years."))

results_10y <- map2(bstraps$splits, event_counts, \(split, count){
  inds <- as.data.frame(split)$ind[1:count]
  cids <- 1:count
  c_spans <- data.table::rbindlist(map2(inds, cids, \(x, y){cluster_spans |> filter(cluster_id == x) |> mutate(cid=y)}))
  
  sup_counts <- num_supports(phi_grid, c_spans |> group_by(cid), angular_domain)
  rate <- tibble(phi = phi_grid, rate = sup_counts / resample_years)
  lambda_f <- approxfun(x = phi_grid, y = rate$rate)
  
  upper_sample_path <- lapply(upper_paths$paths[inds], FUN = function(x){data.frame(phi = x$phi, R = x$R, excess = x$excess)})
  upper_sample_path <- data.table::rbindlist(upper_sample_path) |> filter(excess > 0)
  gpd_fit <- evgam(gpd_formula, data = upper_sample_path, family = "gpd")
  params <- predict(gpd_fit, newdata = data.frame(phi = phi_grid), type = "response")
  rates <- lambda_f(phi_grid)
  tibble(params=params, rates=rates)
}, .progress = TRUE)

saveRDS(results_10y, file = "bootstrap_results_log1p_10y.rds")

result_df <- data.table::rbindlist(lapply(results_10y, FUN=\(x){tibble(phi = phi_grid, scale = x$params$scale, shape = x$params$shape, rate = x$rates)}), idcol = "id")

ret_levels_10y <- data.table::rbindlist(result_df |> group_by(id) |> group_map(\(fit,...){
  data.table::rbindlist(
    joint_return_level_rate(phi_grid, c(10,25,50,100,500), 
                            fit$rate, fit[, c("scale", "shape")], 
                            threshold = "logp1_ald_tau_0_85"), idcol = "Period")
}), idcol = "id") 

# Unstable solution for some bootstrap samples
drop_ids <- unique(ret_levels_10y$id[which(ret_levels_10y$R > 4)])

result_df <- result_df |> filter(!(id %in% drop_ids))
ret_levels_10y <- ret_levels_10y |> filter(!(id %in% drop_ids)) |>
  summarise(R_l = quantile(R, probs=c(0.1)), R_u = quantile(R, probs=c(0.9)),
            R   = mean(R), .by = c(Period, phi))
ret_levels_10y$Period <- as.factor(c(10,25,50,100,500))[ret_levels_10y$Period]

ret_levels_10y[, c("tm2", "hs")] <- ang_trans_inverse(ret_levels_10y |> select(phi, R))
ret_levels_10y[, c("tm2_l", "hs_l")] <- ang_trans_inverse(ret_levels_10y |> select(phi, R=R_l))
ret_levels_10y[, c("tm2_u", "hs_u")] <- ang_trans_inverse(ret_levels_10y |> select(phi, R=R_u))

ret_levels_polygon_10y <- data.table::rbindlist(
  lapply(unique(ret_levels_10y$Period), FUN=\(x){
    data.table::rbindlist(
      list(ret_levels_10y[ret_levels_10y$Period == x, ] |> select(tm2=tm2_l, hs=hs_l, Period),
           (ret_levels_10y[ret_levels_10y$Period == x, ] |> select(tm2=tm2_u, hs=hs_u, Period))[sum(ret_levels_10y$Period==x):1,])
    )}))

summary_10y <- result_df |> group_by(phi) |> 
  summarise(scale_l = quantile(scale, probs=c(0.025)), scale_u = quantile(scale, probs=c(0.975)), 
            shape_l = quantile(shape, probs=c(0.025)), shape_u = quantile(shape, probs=c(0.975)),
            rate_l  = quantile(rate,  probs=c(0.025)), rate_u  = quantile(rate,  probs=c(0.975)),
            scale=mean(scale), shape = mean(shape), rate = mean(rate))

plot <- ggplot(data = summary_10y) + 
  geom_line(mapping=aes(phi, scale)) + 
  geom_line(mapping=aes(phi, scale_l), linetype = "dotdash") + 
  geom_line(mapping=aes(phi, scale_u), linetype = "dotdash") +
  labs(
    x = expression("Pseudo-angle "*phi),
    y = expression("Scale "*sigma),
    title = expression("Bootstrap estimated GP scale parameter "*sigma*" uncertainty, "*T==10*" years."),
    subtitle = expression("Quantiles "*p %in% paste("[", 0.025, ", ", 0.975, "]")*" for "*log1p*" transformation with "*L[1]*" gauge and threshold level "*tau==0.85*".")
  )
ggsave(paste0("plots/chapter 6/bootstraps/log1p_10year_scale_parameter.png"), plot, width = 7, height = 7)

plot <- ggplot(data = summary_10y) + 
  geom_line(mapping=aes(phi, shape)) + 
  geom_line(mapping=aes(phi, shape_l), linetype = "dotdash") + 
  geom_line(mapping=aes(phi, shape_u), linetype = "dotdash") +
  labs(
    x = expression("Pseudo-angle "*phi),
    y = expression("Shape "*xi),
    title = expression("Bootstrap estimated GP scale parameter "*xi*" uncertainty, "*T==10*" years."),
    subtitle = expression("Quantiles "*p %in% paste("[", 0.025, ", ", 0.975, "]")*" for "*log1p*" transformation with "*L[1]*" gauge and threshold level "*tau==0.85*".")
  )
ggsave(paste0("plots/chapter 6/bootstraps/log1p_10year_shape_parameter.png"), plot, width = 7, height = 7)

plot <- ggplot(data=ret_levels_10y, mapping=aes(x=tm2, y=hs)) +
  geom_polygon(data=ret_levels_polygon_10y, mapping = aes(x=tm2, y=hs, fill=Period), alpha = 0.2) +
  geom_point(data=spar_data(spar, space = "original", format = "tibble"),
             mapping=aes(x=tm2, y=hs), shape = ".") +
  geom_path(linewidth = 1,
            mapping = aes(x = tm2, y = hs,
                          color = Period)) +
  labs(
    color = "Return period (years)",
    fill = "Return period (years)",
    x = expression(T[m02] ~ "(s)"),
    y = expression(Hs ~ "(m)"),
    title = expression("Bootstrap uncertainty in estimated return-level contours, "*T==10*" years."),
    subtitle = expression("Quantiles "*p %in% paste("[", 0.1, ", ", 0.9, "]")*" for "*log1p*" transformation with "*L[1]*" gauge and threshold level "*tau==0.85*".")
  )
ggsave(paste0("plots/chapter 6/bootstraps/log1p_10year_return_levels.png"), plot, width = 8, height = 5.5)

bootstrap_param_compare <- data.table::rbindlist(list(summary_df |> mutate(time = "25 Years"), 
                                                      summary_10y |> mutate(time = "10 Years")))
bootstrap_param_compare <- data.table::rbindlist(list(bootstrap_param_compare |> mutate(response = scale,
                                                                                        response_l = scale_l,
                                                                                        response_u = scale_u,
                                                                                        param = "Scale"),
                                                      bootstrap_param_compare |> mutate(response = shape,
                                                                                        response_l = shape_l,
                                                                                        response_u = shape_u, param = "Shape")))
plot <- ggplot(data=bootstrap_param_compare) + geom_line(mapping=aes(x=phi, y=response, color = time)) +
  geom_line(mapping=aes(x=phi, y=response_u, color = time), linetype = "longdash") +
  geom_line(mapping=aes(x=phi, y=response_l, color = time), linetype = "longdash") +
  facet_wrap(vars(param), scales = "free") +
  labs(
    color = "Resample period (years)",
    x = expression("Pseudo-Angle "*phi),
    y = expression("Response"),
    title = expression("Bootstrap comparison of estimated GP parameters, 10 vs 25 years."),
    subtitle = expression("Quantiles "*p %in% paste("[", 0.025, ", ", 0.975, "]")*" for "*log1p*" transformation with "*L[1]*" gauge and threshold level "*tau==0.85*".")
  )
ggsave(paste0("plots/chapter 6/bootstraps/log1p_parameter_comparison.png"), plot, width = 11, height = 5.5)
