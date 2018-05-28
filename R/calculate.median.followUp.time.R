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

calculate.median.followUp.time <- function(fu.time, vital.status, method){
	median.fu.time <- NULL;
	if(method == 1){	# use median survival of alive patients
		alive <- which(vital.status == 0);
		median.fu.time <- median(fu.time[alive]);
		}
	else if(method == 2){	# use Kaplan-Meier potential-follow-up method
		survobj <- Surv(fu.time, (vital.status==0));
		km.reverse <- survfit(survobj ~ 1);
		median.fu.time <- summary(km.reverse)$table['median'];
		}
	else{
		print("Incorrect parameter value for 'method'. Must specify 1 or 2.");
		}
	return(median.fu.time);
	}
