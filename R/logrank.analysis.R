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

logrank.analysis <- function(survival.object, groups){

	# perform logrank to assess survival differences
	stats <- survdiff(formula = survival.object ~ groups);

	# define statistical parameters
	this.pval <- pchisq(stats$chisq, df = (length(levels(as.factor(groups))) - 1), lower.tail = FALSE);
	expected.groups <- stats$exp[1:length(levels(as.factor(groups)))];

	# prepare summary statistical values to return, if requested
	summary.stats <- data.frame(
		statistical.method = "logrank",
		pvalue = this.pval,
		groups = levels(as.factor(groups)),
		observed.groups = stats$obs,
		expected.groups = stats$exp
		);

	# Remove ret.stats rownames to make returned value look nicer
#	rownames(summary.stats) <- NULL;

	# Return summary.stats
	return(summary.stats)
	}
