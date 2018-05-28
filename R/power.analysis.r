power.analysis <- function(event.number, HR, alpha = 0.05) {

	warning("This function has been moved to BoutrosLab.statistics.power and renamed power.cox.univariate");
	zalpha <- qnorm(1 - alpha/2);
	zbeta <- sqrt(event.number) * log(HR) / 2 - zalpha;
	power <- pnorm(zbeta);

	return(power)

	}
