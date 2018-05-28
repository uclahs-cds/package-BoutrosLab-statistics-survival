# The BoutrosLab.statistics.survival package is copyright (c) 2013 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

fit.coxmodel <- function(groups, survobj, other.data = NULL, stratification.factor = NULL, stratification.value = NULL, rounding = 3, return.cox.model = FALSE) {

	# verify that the variables to be passed to the coxph function have the same length
	if (nrow(survobj) != length(groups)) {
		stop('The objects passed to survobj and groups in the fit.coxmodel function do not have consistent dimensions.')
		}
	if (!is.null(other.data)){
		if (nrow(other.data) != nrow(survobj)){
			stop('The objects passed to survobj and other.data in the fit.coxmodel function do not have consistent dimensions.')
			}
		}
	if (!is.null(stratification.factor)){
		if (length(stratification.factor) != nrow(survobj)) {
			stop('The objects passed to survobj and stratification.factor in the fit.coxmodel function do not have consistent dimensions.')
			}
		}
	if (!is.null(stratification.factor) & !is.null(stratification.value)){
		if (stratification.value > max(stratification.factor) | stratification.value < min(stratification.factor)){
			stop('The stratification value is outside the range of the variable');
			}
		}
	
	# define data frame containing all variables for which we want to estimate HRs
	if (!is.null(other.data)){
		variables <- data.frame(groups = groups, other.data);	
		}
	else {
		variables <- data.frame(groups = groups);
		}

	# fit cox model
	if (!is.null(stratification.factor) & !is.null(stratification.value)){
		cox.model <- coxph(formula = survobj ~ . + strata(stratification.factor > stratification.value), data = variables);
		}
	else if (!is.null(stratification.factor)){
		cox.model <- coxph(formula = survobj ~ . + strata(stratification.factor), data = variables);
		}
	else {
		cox.model <- coxph(formula = survobj ~ ., data = variables);
		}

	summary.cox.model <- summary(cox.model);
	this.hr <- summary.cox.model$conf.int[1,1];
	this.95l <- summary.cox.model$conf.int[1,3];
	this.95u <- summary.cox.model$conf.int[1,4];
	this.pval <- 2*pnorm(-abs(summary.cox.model$coef[1,4]));
	this.n <- length(groups[!is.na(groups)]);

	# round the major values to a few decimal places
	if (rounding) {
		this.hr  <- round(this.hr,  digits = rounding);
		this.95l <- round(this.95l, digits = rounding);
		this.95u <- round(this.95u, digits = rounding);
		}

	if (TRUE == return.cox.model) {
		return(cox.model);
		}
	else {
		return(c(this.hr, this.95l, this.95u, this.pval, this.n));
		}

	}
