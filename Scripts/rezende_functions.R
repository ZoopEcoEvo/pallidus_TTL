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


	dynamic.landscape <- function(ta,tolerance.landscape){
			surv <- tolerance.landscape$S[,2]
			ta.mn <- tolerance.landscape$ta.mn
			z <- tolerance.landscape$z
			shift <- 10^((ta.mn - ta)/z)	
			time.rel <- 0
			alive <- 100
			for(i in 1:length(ta)){			
				if(alive[length(na.omit(alive))] > 0) {							
					alive <- try(c(alive,approx(c(0,shift[i]*surv),seq(100,0,length.out = length(c(0,surv))),xout = time.rel[i] + 1)$y),silent=TRUE)		
					time.rel <- try(c(time.rel,approx(seq(100,0,length.out = length(c(0,surv))),c(0,shift[i + 1]*surv),xout = alive[i + 1])$y),silent=TRUE)}
				else{
					alive <- 0}}				
			out <- data.frame(cbind(ta=ta[1:(length(alive)-1)],time=(1:length(ta))[1:(length(alive)-1)],alive=alive[1:(length(alive)-1)]))
			par(mar=c(4,4,1,1),mfrow=c(1,2))
			plot(1:length(ta),ta,type="l",xlim=c(0,length(ta)),ylim=c(min(ta),max(ta)),col="black",lwd=1.5,las=1,
				xlab = "Time (min)", ylab = "Temperature (ºC)")			
			plot(out$time,out$alive,type="l",xlim=c(0,length(ta)),ylim=c(0,100),col="black",lwd=1.5,las=1,
				xlab = "Time (min)", ylab = "Survival (%)")
			list(time = out$time,ta = out$ta, alive = out$alive)}

 

# Function to estimate a thermal tolerance landscape directly from a fitted
# AFT model (survreg with dist "loglogistic", "lognormal", or "weibull"),
# instead of individual median survival times per temperature/population.
#
# survreg models log(T) = X*beta + scale*W, where W follows a standard
# distribution set by `dist` (logistic, normal, or extreme value). CTmax and
# z come directly from the model's intercept/slope; the reference survival
# curve S is obtained analytically from that distribution's quantile
# function at the mean temperature, rather than by interpolating raw data.
# Coefficient (and, when estimated, log-scale) uncertainty is propagated
# into ctmax, z, and S by simulating from the model's asymptotic sampling
# distribution (vcov(model)) and taking quantiles of the simulated values.

	qstd <- function(dist, p) {
		switch(dist,
			loglogistic = qlogis(p),
			lognormal   = qnorm(p),
			weibull     = log(-log(1 - p)),
			stop("Unsupported distribution: ", dist))}

	model_ttl <- function(model, temp_var, ta.mn = NULL, nsim = 2000,
	                       conf.level = 0.95, seed = NULL) {

		dist  <- model$dist
		coefs <- coef(model)
		vc    <- vcov(model)

		beta0_name <- "(Intercept)"
		if (!(beta0_name %in% names(coefs)) || !(temp_var %in% names(coefs))) {
			stop("`model` must contain an intercept and a '", temp_var, "' coefficient.")}

		has_log_scale <- "Log(scale)" %in% rownames(vc)
		point_params  <- if (has_log_scale) c(coefs, "Log(scale)" = log(model$scale)) else coefs

		if (is.null(ta.mn)) {
			mf <- model.frame(model)
			ta.mn <- mean(unique(mf[[temp_var]]))}

		beta0 <- point_params[[beta0_name]]
		slope <- point_params[[temp_var]]
		scale <- model$scale

		ctmax <- as.numeric(-beta0 / slope)
		z     <- as.numeric(-1 / slope)

		eta_mn    <- as.numeric(beta0 + slope * ta.mn)
		surv_p    <- seq(0.999, 0.001, -0.001)
		log_time  <- eta_mn + scale * qstd(dist, 1 - surv_p)
		time_pred <- exp(log_time)

		# Simulate from the joint asymptotic distribution of the parameters
		# to propagate uncertainty into ctmax, z, and S
		if (!is.null(seed)) set.seed(seed)
		sims <- MASS::mvrnorm(nsim, mu = point_params, Sigma = vc)
		sim_beta0 <- sims[, beta0_name]
		sim_slope <- sims[, temp_var]
		sim_scale <- if (has_log_scale) exp(sims[, "Log(scale)"]) else scale

		sim_ctmax <- -sim_beta0 / sim_slope
		sim_z     <- -1 / sim_slope

		sim_eta_mn   <- sim_beta0 + sim_slope * ta.mn
		sim_log_time <- outer(sim_eta_mn, rep(1, length(surv_p))) +
		                 outer(sim_scale, qstd(dist, 1 - surv_p))
		sim_time <- exp(sim_log_time)

		alpha <- 1 - conf.level
		probs <- c(alpha / 2, 1 - alpha / 2)

		ctmax.ci <- as.numeric(quantile(sim_ctmax, probs = probs, na.rm = TRUE))
		z.ci     <- as.numeric(quantile(sim_z, probs = probs, na.rm = TRUE))
		time.ci  <- apply(sim_time, 2, quantile, probs = probs, na.rm = TRUE)

		S <- data.frame(surv = surv_p, time = time_pred,
		                 time.lwr = time.ci[1, ], time.upr = time.ci[2, ])

		list(ctmax = ctmax, ctmax.ci = ctmax.ci,
		     z = z, z.ci = z.ci,
		     ta.mn = ta.mn, S = S,
		     conf.level = conf.level, nsim = nsim)}

	