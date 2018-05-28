parse.cox.table <- function(cox.data, table.rownames = NULL, zph.p = NULL){
	summary.cox <- summary(cox.data);
	if(is.null(zph.p)){
		cox.zph.table = cox.zph(cox.data, transform = 'identity')$table;
		if(1 < nrow(cox.zph.table)){
			zph.p <- cox.zph.table[-nrow(cox.zph.table), 'p'];	# take all but last element (This is the GLOBAL p-value)
			}
		else{
			zph.p <- cox.zph.table[1, 'p'];	# take all but last element (This is the GLOBAL p-value)
			}
		}
	table.out <- data.frame(
		HR = summary.cox$conf.int[, 1],
		CI95low = summary.cox$conf.int[, 3],
		CI95high = summary.cox$conf.int[,4],
		p = summary.cox$coef[, 5],
		ph.p = zph.p
		);

	# add the Wald and Log rank p-values to the last row
	table.out <- rbind(table.out,
					c('wald.p', summary.cox$waldtest[3], 'logrank', summary(cox.data)$logtest[3], '')
					);

	if(! is.null(table.rownames)){
		rownames(table.out) <- table.rownames;
		}

	return(table.out);
	}
