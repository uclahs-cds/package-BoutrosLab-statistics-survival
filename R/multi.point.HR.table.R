multi.point.HR.table <- function(all.groups = NULL, all.survtime = NULL, all.survstat = NULL, truncation.thresholds = c(5, 10), covariates = NULL) {

	# sanity checks:
	# remove NA data points
	complete.data.points <- complete.cases(data.frame(all.groups, all.survtime, all.survstat));
	all.groups <- all.groups[complete.data.points];
	all.survtime <- all.survtime[complete.data.points];
	all.survstat <- all.survstat[complete.data.points];

	if (length(all.survtime[all.survtime == 0]) > 0) {
		cat("\nWarning: survtime for some samples is zero, setting it to 0.000001");
		all.survtime[all.survtime == 0] <- 0.000001;
		}

	# make sure max surv time is less than max cut point
	max.surv <- max(all.survtime, na.rm = TRUE);
	if(max.surv < max(truncation.thresholds)){
		new.cut.points <- c(max.surv / 2, max.surv);
		truncation.thresholds <- new.cut.points;
		print(paste("Warning: max survtime was less than max cut point. New cut points are", new.cut.points));
		}

	tmp.data <- as.data.frame( cbind("time" = all.survtime, "status" = all.survstat, "fu.time" = all.survtime) );

	# only extract these stats from each variable in the cox models
	stats.to.keep <- c('HR', 'CI95l', 'CI95h', 'cox.p', 'ph.p');

	# intervalise the dataset
	extended.cox.data <- survSplit(
		formula = Surv(time, status) ~ .,
		data = tmp.data,
		cut = truncation.thresholds,
		end = "time",
		start = "start",
		event = "status"
		);

	# format as desired, stats as columns, and chunks of columns denote different variables from the Cox model
	if(is.null(covariates)){	# no covariates
		out.data <- matrix(
			nrow = length(truncation.thresholds),	# one line per cut point
			ncol = (3 + length(stats.to.keep))	# keep track of all stats for the group variable, and n, wald.p and global PH p
			);
		}
	else{	# need to add extra columns for each covariate
		# determine the number of levels per covariate
		n.levels <- 1;	# for the groups variable
		for(x in 1:ncol(covariates)){
			if(length(levels(covariates[,x])) > 0){
				n.levels <- n.levels + length(levels(covariates[,x])) - 1;
				covariates[,x] <- as.factor(covariates[,x]);
				}
			else{
				n.levels <- n.levels + 1;
				}
			}

		# keep track of all stats for the group variable and for each covariate level, and n, wald.p and global PH p
		out.data <- matrix(
			nrow = length(truncation.thresholds),	# one line per cut point
			ncol = (3 + length(stats.to.keep) * (n.levels))	
			);
		}

	# create groups vector by adjusting for patients that have an event SO FAR
	revised.groups <- all.groups;
	revised.covariates <- covariates;
	to.remove <- NULL;
	step.size <- 0;

	for (index in 1:length(truncation.thresholds)) {
		truncation.threshold <- truncation.thresholds[index];
		# start time is previous step
		start.time <- step.size;

		interval.data <- extended.cox.data[
			which(extended.cox.data$start == start.time & extended.cox.data$time <= truncation.threshold),
			];
		surv.time <- interval.data$time - interval.data$start;
		surv.stat <- interval.data$status;

		if (length(to.remove) > 0) {
			# cat("REMOVING: ", to.remove); cat(" BEFORE: ", length(revised.groups));
			revised.groups <- revised.groups[-to.remove];
			if(!is.null(revised.covariates)){
				revised.covariates <- revised.covariates[-to.remove , ];
				}
			}

		to.remove <- which(interval.data$status == 1 & interval.data$time <= truncation.threshold);
		to.remove <- c(to.remove, which(interval.data$status == 0 & interval.data$fu.time <= truncation.threshold));

		# fit cox model and keep the actual model so that we can extract info on covariates if necessary
		if(is.null(covariates)){
			# no covariates, very easy just fit cox model add results to table
			coxfit <- fit.coxmodel(
				groups = revised.groups,
				survobj = Surv(surv.time, surv.stat),
				rounding = 4,
				return.cox.model = TRUE
				);
			coxmodel <- summary(coxfit);

			# Perform cox.zph test on cox.model
			cox.zph.test <- cox.zph(coxfit, transform = "identity");
		
			# extract summary characteristics from the Cox fits and PH test
			n 			<-  length(revised.groups[!is.na(revised.groups)]);
			wald.p 		<- coxmodel$waldtest[3];
			global.ph.p <- cox.zph.test$table[nrow(cox.zph.test$table), "p"];	# save the GLOBAL p-value.

			out.data[index , ] <- c(n, wald.p, global.ph.p, coxmodel$conf.int[c(1,3,4)], coxmodel$coef[5], cox.zph.test$table[1, 'p']);
			}
		# we do have covariates, need to do some reformatting of the Cox output to put everything on one line
		else{
			if(length(unique(surv.stat)) == 1){		# all patients are in the same survival state, set  everything to NA
				out.data[index,] <- rep(NA, ncol(out.data));
				}
			else{
				coxfit <- fit.coxmodel(
					groups = revised.groups,
					survobj = Surv(surv.time, surv.stat),
					other.data = as.data.frame(revised.covariates),
					rounding = 4,
					return.cox.model = TRUE
					);
				coxmodel <- summary(coxfit);

				# Perform cox.zph test on cox.model
				cox.zph.test <- cox.zph(coxfit, transform = "identity");

				# extract summary characteristics from the Cox fits and PH test
				n 			<-  length(revised.groups[!is.na(revised.groups)]);
				wald.p 		<- coxmodel$waldtest[3];
				global.ph.p <- cox.zph.test$table[nrow(cox.zph.test$table), "p"];	# save the GLOBAL p-value.

				cox.data <- data.frame(
					hr		= coxmodel$conf.int[,1],
					ci95l	= coxmodel$conf.int[,3],
					ci95u	= coxmodel$conf.int[,4],
					pval	= coxmodel$coef[,5],
					ph.p 	= cox.zph.test$table[-nrow(cox.zph.test$table), "p"]	# take all but last element (This is the GLOBAL p-value)
					);

				# format as desired, stats as columns, and chunks of columns denote different variables from the Cox model
				new.row <- c(n, wald.p, global.ph.p);
				for(i in 1:nrow(cox.data)){
					new.row <- as.vector(c(new.row, as.vector(cox.data[i,])));
					}
				# add current stats to table
				out.data[index , ] <- unlist(new.row);
				}
			}

		# update step size for next iteration
		step.size <- truncation.threshold;
		}

	# add row names to table and return it
	rownames(out.data) <- c(paste0('Time', truncation.thresholds));
	if(is.null(covariates)){
		# make column names specific to the variable and stat combination
		colnames(out.data) <- 	c(
			'n',
			'wald.p',
			'global.ph.p',
			paste('group', stats.to.keep, sep = "_")
			);
		}
	else{
		colnames(out.data) <- 	c(
			'n',
			'wald.p',
			'global.ph.p',
			paste(
				unlist(lapply(rownames(coxmodel$conf.int), function(f) rep(f, length(stats.to.keep)))),
				stats.to.keep,
				sep = "_"
				)
			);
		}

	return(out.data);
	}
