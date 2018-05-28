# The BoutrosLab.statistics.survival package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### ph.fails.R ######################################################################
# Description:
#       Tests cox.zph assumption for a cox model, and returns warnings and/or residual plots.
# Input variables:
#       
# Output variables: 
#	
#     

ph.fails <- function(cox.model, cox.zph.threshold = 0.1, pvalues = FALSE) {
	
	# Perform cox.zph test on cox.model
	cox.zph.test <- cox.zph(cox.model, transform = "identity");

	# Get number of betas estimated in the cox model, in order to be able to pull out
	# appropriate pvalues later in the code
	if ("GLOBAL" %in% row.names(cox.zph.test$table)) {
		number.of.betas <- nrow(cox.zph.test$table) - 1;
		row.number.for.GLOBAL <- nrow(cox.zph.test$table);
		}
	else {
		number.of.betas <- nrow(cox.zph.test$table);
		row.number.for.GLOBAL <- NA;
		}

	# Get vector of all pvalues computed for cox.zph test, one for each beta and one for the overall model (global)
	pvalues.zph <- cox.zph.test$table[,'p'];

	# Get vector of names of covariates that fail the PH assumption;
	covariates.failing.PH <- rownames(cox.zph.test$table)[pvalues.zph < cox.zph.threshold]	

	if (TRUE == pvalues){
		return(pvalues.zph);
		}
	else {
		# If one or more covariates fails the PH assumption, we return TRUE
		return(length(covariates.failing.PH) > 0);
		}

	}

		


