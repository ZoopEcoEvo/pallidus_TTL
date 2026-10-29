# Function to estimate thermal tolerance landscape from static assays
# General procedure
# Step 1: Calculate CTmax and z from TDT curve
# Step 2: Calculate average log10 time and Ta (mean x and y for interpolation purposes)
# Step 3: Interpolating survival probabilities to make them comparable across treatments 
# Step 4: Overlap all survival curves into a single one by shifting each curve to mean x and y employing z
# Step 5: Build expected survival curve with mean x and y pooling all data
# Step 6: Expand expected curve to multiple Ta (with 0.1ºC difference for predictive purposes)

	tolerance.landscape <- function(ta,time){
	
			data <- data.frame(ta,time)
			data <- data[order(data$ta,data$time),]

			# Step 1: Calculate CTmax and z from TDT curve
			ta <- as.numeric(levels(as.factor(data$ta)))			
			model <- lm(log10(data$time) ~ data$ta); summary(model)
			ctmax <- -coef(model)[1]/coef(model)[2]
			z <- -1/coef(model)[2]

			# Step 2: Calculate average log10 time and Ta (mean x and y for interpolation purposes)
			time.mn <- mean(log10(data$time))
			ta.mn <- mean(data$ta)
	
			# Step 3: Interpolating survival probabilities to make them comparable across treatments 
			time.interpol <- matrix(,1001,length(ta))
			for(i in 1:length(ta)){	
			time <- c(0,sort(data$time[data$ta==ta[i]])); p <- seq(0,100,length.out = length(time))
			time.interpol[,i] <- approx(p,time,n = 1001)$y}			

			# Step 4: Overlap all survival curves into a single one by shifting each curve to mean x and y employing z
			# Step 5: Build expected survival curve with median survival time for each survival probability
			shift <- (10^((ta - ta.mn)/z))
			time.interpol.shift <- t(t(time.interpol)*shift)[-1,]
			surv.pred <- 10^apply(log10(time.interpol.shift),1,median) 	
	
			# Step 6: Expand predicted survival curves to measured Ta (matrix m arranged from lower to higher ta)
			# Step 7: Obtain predicted values comparable to each empirical measurement
			m <- surv.pred*matrix ((10^((ta.mn - rep(ta, each = 1000))/z)), nrow = 1000)
			out <-0
			for(i in 1:length(ta)){
				time <- c(0,data$time[data$ta==ta[i]]); p <- seq(0,100,length.out = length(time))
				out <- c(out,approx(seq(0,100,length.out = 1000),m[,i],xout=p[-1])$y)}
				data$time.pred <- out[-1]
				colnames(m) <- paste("time.at",ta,sep=".")
				m <- cbind(surv.prob=seq(1,0.001,-0.001),m)

			par(mfrow=c(1,2),mar=c(4.5,4,1,1),cex.axis=1.1)
			plot(-10,-10,las=1,xlab="Time (min)",ylab="Survival (%)",col="white",xaxs="i",yaxs="i",xlim=c(0,max(data$time)*1.05),ylim=c(0,105))
				for(i in 1:length(ta)){
				time <- c(0,sort(data$time[data$ta==ta[i]])); p <- seq(100,0,length.out = length(time))
				points(time,p,pch=21,bg="black",cex=0.5)
				time <- c(0,sort(data$time.pred[data$ta==ta[i]]))
				points(m[,i+1],100*m[,1],type="l",lty=2)}
				segments(max(data$time)*0.7,90,max(data$time)*0.8,90,lty=2)
				text(max(data$time)*0.82,90,"fitted",adj=c(0,0.5))
			plot(log10(data$time.pred),log10(data$time),pch=21,bg="black",cex=0.5,lwd=0.7,las=1,xlab="Fitted Log10 time",ylab="Measured Log10 time")
			abline(0,1,lty=2)
			rsq <- round(summary(lm(log10(data$time) ~ log10(data$time.pred)))$r.square,3)
			text(min(log10(data$time.pred)),max(log10(data$time)),substitute("r"^2*" = "*rsq),adj=c(0,1))
			list(ctmax = as.numeric(ctmax), z = as.numeric(z), ta.mn = ta.mn,  S = data.frame(surv=seq(0.999,0,-0.001),time=surv.pred),
			time.obs.pred=cbind(data$time,data$time.pred), rsq = rsq)}					

 

# Function to estimate survival probability from tolerance landscapes and environmental temperature data
#
# Rewritten (2026-10-29) from Rezende et al.'s original loop. The original
# tracked survival minute by minute with repeated approx() calls wrapped in
# try(); once survival ran past the end of the reference curve it returned NA
# instead of 0 for the rest of the series.
#
# Because every survival curve is a time-scaled copy of the reference curve
# (the TDT / AFT assumption), one minute at temperature T is equivalent to
# 1 / shift(T) minutes at the reference temperature ta.mn, where
# shift(T) = 10^((ta.mn - T) / z). Accumulated "reference-equivalent" time is
# therefore cumsum(1 / shift), and survival is the reference curve evaluated
# at that time. This gives the same result as the original algorithm, is
# vectorised, and drops to 0 (rather than NA) once survival is exhausted.
#
# Arguments
#   ta                  temperature for each minute (deg C)
#   tolerance.landscape list with $S (data.frame: surv, time at ta.mn),
#                       $ta.mn and $z (z in log10 units, deg C)
#   plot                draw Rezende's two-panel base plot (default TRUE, as before)
# Value: list(time, ta, alive) with alive = % surviving at the END of each minute

	dynamic.landscape <- function(ta, tolerance.landscape, plot = TRUE){
			S     <- tolerance.landscape$S
			ta.mn <- tolerance.landscape$ta.mn
			z     <- tolerance.landscape$z
			shift <- 10^((ta.mn - ta)/z)              # time multiplier relative to ta.mn
			tau   <- cumsum(1/shift)                  # accumulated reference-equivalent minutes
			# reference curve: survival 1 at time 0, then the tabulated quantiles, 0 beyond the last one
			ref_time <- c(0, S$time)
			ref_surv <- c(1, S$surv)
			alive <- 100 * approx(ref_time, ref_surv, xout = tau, yleft = 1, yright = 0, ties = "ordered")$y
			out <- list(time = seq_along(ta), ta = ta, alive = alive)
			if (plot) {
				par(mar=c(4,4,1,1),mfrow=c(1,2))
				plot(out$time,ta,type="l",xlim=c(0,length(ta)),ylim=range(ta),col="black",lwd=1.5,las=1,
					xlab = "Time (min)", ylab = "Temperature (ºC)")
				plot(out$time,out$alive,type="l",xlim=c(0,length(ta)),ylim=c(0,100),col="black",lwd=1.5,las=1,
					xlab = "Time (min)", ylab = "Survival (%)")}
			out}



# Distribution helpers shared by the AFT-based functions below.
# AFT models give log(T) = location + scale * W, with W a standard logistic
# (loglogistic), normal (lognormal) or minimum extreme value (weibull).

	qstd <- function(dist, p) {
		switch(dist,
			loglogistic = qlogis(p),
			lognormal   = qnorm(p),
			weibull     = log(-log(1 - p)),
			stop("Unsupported distribution: ", dist))}

	pstd <- function(dist, w) {
		switch(dist,
			loglogistic = plogis(w),
			lognormal   = pnorm(w),
			weibull     = 1 - exp(-exp(w)),
			stop("Unsupported distribution: ", dist))}



# Build a tolerance landscape from draws of the AFT parameters.
#
# Model (natural-log time, as fitted by survreg / brms):
#   log(T) = beta0 + slope * (temp - temp_center) + sigma * W
# Here beta0 is the location (log-time) at temp_center.
#
# UNITS: Rezende's z is defined on the log10 scale (the deg C needed for a
# 10-fold change in survival time). With a natural-log slope this is
#   z = -log(10) / slope
# NOT -1 / slope. (The earlier version of model_ttl() used -1/slope, which
# made z 2.3x too small and made dynamic.landscape() exaggerate temperature
# effects away from ta.mn.)
# CTmax is the temperature giving a 1-min location (log T = 0):
#   CTmax = temp_center - beta0 / slope
# For log-logistic and log-normal models the location is the median, so CTmax
# is the 1-min LT50 exactly as in Rezende et al.
#
# Arguments
#   beta0, slope, sigma  numeric vectors of draws (equal length)
#   point                named numeric c(beta0, slope, sigma) used for point
#                        estimates (default: posterior/sampling medians)
#   dist                 "loglogistic", "lognormal" or "weibull"
#   temp_center          temperature at which beta0 applies
#   ta.mn                reference temperature for the S curve (default temp_center)
#   lt_times             LT50 times (min) at which to report the temperature
#   conf.level           width of the reported intervals

	ttl_from_draws <- function(beta0, slope, sigma, dist, temp_center,
	                           point = NULL, ta.mn = temp_center,
	                           lt_times = c(10, 60, 720), conf.level = 0.95) {
		if (is.null(point)) point <- c(beta0 = median(beta0), slope = median(slope), sigma = median(sigma))
		probs <- c((1 - conf.level)/2, 1 - (1 - conf.level)/2)
		qq <- function(x) as.numeric(quantile(x, probs, na.rm = TRUE))

		ctmax_fun <- function(b0, b1) temp_center - b0/b1
		z_fun     <- function(b1) -log(10)/b1
		# temperature at which the median survival time equals t (log-logistic/log-normal)
		temp_at   <- function(b0, b1, t) temp_center + (log(t) - b0)/b1

		sim_ctmax <- ctmax_fun(beta0, slope)
		sim_z     <- z_fun(slope)

		# reference survival curve at ta.mn
		surv_p   <- seq(0.999, 0.001, -0.001)
		loc_pt   <- point[["beta0"]] + point[["slope"]] * (ta.mn - temp_center)
		time_pt  <- exp(loc_pt + point[["sigma"]] * qstd(dist, 1 - surv_p))
		sim_loc  <- beta0 + slope * (ta.mn - temp_center)
		sim_time <- exp(outer(sim_loc, rep(1, length(surv_p))) + outer(sigma, qstd(dist, 1 - surv_p)))
		time.ci  <- apply(sim_time, 2, quantile, probs = probs, na.rm = TRUE)

		lt_temps <- do.call(rbind, lapply(lt_times, function(t) {
			s <- temp_at(beta0, slope, t)
			data.frame(lt50_min = t, temp = temp_at(point[["beta0"]], point[["slope"]], t),
			           lwr = qq(s)[1], upr = qq(s)[2])}))

		list(ctmax = ctmax_fun(point[["beta0"]], point[["slope"]]), ctmax.ci = qq(sim_ctmax),
		     z = z_fun(point[["slope"]]), z.ci = qq(sim_z),
		     sigma = point[["sigma"]], dist = dist,
		     ta.mn = ta.mn, lt_temps = lt_temps,
		     S = data.frame(surv = surv_p, time = time_pt, time.lwr = time.ci[1, ], time.upr = time.ci[2, ]),
		     draws = data.frame(ctmax = sim_ctmax, z = sim_z, beta0 = beta0, slope = slope, sigma = sigma),
		     temp_center = temp_center, conf.level = conf.level)}



# Tolerance landscape from a fitted survreg AFT model (no random effects).
# Uncertainty is propagated by simulating from the asymptotic sampling
# distribution (vcov(model)). `temp_var` must be the (possibly centred)
# temperature term; give `temp_center` if the covariate was centred.

	model_ttl <- function(model, temp_var, temp_center = 0, ta.mn = NULL, nsim = 2000,
	                       conf.level = 0.95, seed = NULL, ...) {
		coefs <- coef(model)
		vc    <- vcov(model)
		if (!("(Intercept)" %in% names(coefs)) || !(temp_var %in% names(coefs)))
			stop("`model` must contain an intercept and a '", temp_var, "' coefficient.")
		has_log_scale <- "Log(scale)" %in% rownames(vc)
		mu <- if (has_log_scale) c(coefs, "Log(scale)" = log(model$scale)) else coefs

		if (is.null(ta.mn)) ta.mn <- temp_center + mean(unique(model.frame(model)[[temp_var]]))

		if (!is.null(seed)) set.seed(seed)
		sims <- MASS::mvrnorm(nsim, mu = mu, Sigma = vc)
		sim_sigma <- if (has_log_scale) exp(sims[, "Log(scale)"]) else rep(model$scale, nsim)

		ttl_from_draws(beta0 = sims[, "(Intercept)"], slope = sims[, temp_var], sigma = sim_sigma,
		               point = c(beta0 = coefs[["(Intercept)"]], slope = coefs[[temp_var]], sigma = model$scale),
		               dist = model$dist, temp_center = temp_center, ta.mn = ta.mn,
		               conf.level = conf.level, ...)}



# Tolerance landscape for one population from a brms AFT model fitted as
#   y | cens(cens, y2) ~ 0 + population + population:temp_c + (1 | run)
# with family lognormal() (or weibull()), where temp_c = temp - temp_center.
# Uses the population-level (fixed) effects, i.e. a typical run; the run
# random effect is treated as experimental (bath/day) noise and is not added
# to the survival curve. Each posterior draw is carried through.

	brms_ttl <- function(fit, population, temp_var = "temp_c", temp_center,
	                     pop_var = "population", ta.mn = temp_center, ndraws = NULL, ...) {
		dr <- as.data.frame(fit)
		if (!is.null(ndraws)) dr <- dr[sample(nrow(dr), min(ndraws, nrow(dr))), ]
		b0_name <- paste0("b_", pop_var, population)
		b1_cand <- c(paste0("b_", pop_var, population, ":", temp_var), paste0("b_", temp_var, ":", pop_var, population))
		b1_name <- b1_cand[b1_cand %in% names(dr)][1]
		if (!(b0_name %in% names(dr)) || is.na(b1_name))
			stop("Could not find ", b0_name, " and a population-specific slope. Fit the model with ",
			     "`0 + population + population:", temp_var, "`.")
		fam <- family(fit)$family
		if (fam == "weibull") {
			# brms weibull: mu is the mean on the identity scale with a log link; convert
			# to the survreg-style location/scale (log T = log(scale_w) + (1/shape) * W)
			shape <- dr$shape
			loc_shift <- -lgamma(1 + 1/shape)   # log(scale_w) = log(mean) - lgamma(1 + 1/shape)
			return(ttl_from_draws(beta0 = dr[[b0_name]] + loc_shift, slope = dr[[b1_name]], sigma = 1/shape,
			                      dist = "weibull", temp_center = temp_center, ta.mn = ta.mn, ...))}
		if (fam != "lognormal") stop("brms_ttl() supports lognormal and weibull families.")
		ttl_from_draws(beta0 = dr[[b0_name]], slope = dr[[b1_name]], sigma = dr$sigma,
		               dist = "lognormal", temp_center = temp_center, ta.mn = ta.mn, ...)}



# Predicted survival through a temperature series with uncertainty: runs
# dynamic.landscape() once per parameter draw of a tolerance landscape built
# by ttl_from_draws() / model_ttl() / brms_ttl().
#
# Value: data.frame(time, ta, alive (point estimate), lwr, upr)

	dynamic_landscape_draws <- function(ta, tl, ndraws = 500, conf.level = 0.95, seed = NULL) {
		if (!is.null(seed)) set.seed(seed)
		d <- tl$draws[sample(nrow(tl$draws), min(ndraws, nrow(tl$draws))), ]
		surv_p <- tl$S$surv
		w <- qstd(tl$dist, 1 - surv_p)
		sims <- vapply(seq_len(nrow(d)), function(i) {
			loc <- d$beta0[i] + d$slope[i] * (tl$ta.mn - tl$temp_center)
			tl_i <- list(S = data.frame(surv = surv_p, time = exp(loc + d$sigma[i] * w)), ta.mn = tl$ta.mn, z = d$z[i])
			dynamic.landscape(ta, tl_i, plot = FALSE)$alive}, numeric(length(ta)))
		probs <- c((1 - conf.level)/2, 1 - (1 - conf.level)/2)
		data.frame(time = seq_along(ta), ta = ta,
		           alive = dynamic.landscape(ta, tl, plot = FALSE)$alive,
		           lwr = apply(sims, 1, quantile, probs[1]),
		           upr = apply(sims, 1, quantile, probs[2]))}



# Survival probability surface S(T, t) implied by an AFT tolerance landscape,
# for plotting (e.g. the 3D landscape figure).
#   log10 LT50(T) = (ctmax - T) / z ; S = 1 - F((log t - log LT50) / sigma)

	ttl_surface <- function(tl, temps, times) {
		outer(times, temps, function(t, T) {
			loc <- log(10) * (tl$ctmax - T) / tl$z
			1 - pstd(tl$dist, (log(t) - loc) / tl$sigma)})}
