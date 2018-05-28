# The BoutrosLab.statistics.survival package is copyright (c) 2010 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

get.survival.rate <- function(surv.obj, groups, cut.points){
	survival.rates <- matrix(nrow = length(levels(groups)), ncol = length(cut.points));
	survfit.obj <- summary(survfit(surv.obj ~ groups));
	print(survfit.obj);
	# make sure cut points are less than the max survival time
	for(x in 1:length(cut.points)){
		if(cut.points[x] > max(as.vector(surv.obj))){
			print(cut.points[x]);
			print(max(as.vector(surv.obj)));
			stop("Cut points are larger than the maximum survival time");
			}
		}

	for(i in 1:nrow(survival.rates)){
		indices <- survfit.obj$strata == levels(survfit.obj$strata)[i];

		# case where no patients have an event in the current group i
		if(length(which(indices)) == 0){
			survival.rates[i, ] <- 1;
			}
		else{
			for(j in 1:ncol(survival.rates)){
				# check if there is a time-point exactly at the cut point 
				perfect.match <- which(survfit.obj$time[indices] == cut.points[j]);
				if(length(perfect.match) > 0){
					survival.rates[i,j] <- survfit.obj$surv[indices][perfect.match];
					}
				# possible that cut point is smaller than the first time point. If so, return 1.
				else if(survfit.obj$time[indices][1] > cut.points[j]){
					survival.rates[i, j] <- 1;
					}
				else{
					# choose the first time point greater than the cut point, and substract one ( we want the point just before the cut point)
					cur.time.pt <- which(survfit.obj$time[indices] > cut.points[j])[1];
					if(! is.na(cur.time.pt)){	# cut point is within range of survival object
						survival.rates[i,j] <- survfit.obj$surv[indices][cur.time.pt - 1];
						}
					# use the last index then
					else{
						survival.rates[i, j] <- survfit.obj$surv[indices][length(which(indices))];
						}
					}
				}
			}
		}
	survival.rates <- as.data.frame(survival.rates);
	rownames(survival.rates) <- levels(survfit.obj$strata);
	colnames(survival.rates) <- paste0('time', cut.points);

	return(survival.rates);
	}
