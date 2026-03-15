functions {
	real bivarwbl_log(matrix xy, real th1, real th2, real b1, real b2, real dlt) {
        	vector[rows(xy)] prob;
		real out;
		for (i in 1:rows(xy)) {
			prob[i] <- ((b1/th1) * (xy[i, 1]/th1)^((b1/dlt)-1)) * ((b2/th2) * (xy[i, 2]/th2)^((b2/dlt)-1)) * ((xy[i, 1]/th1)^(b1/dlt) + (xy[i, 2]/th2)^(b2/dlt))^(dlt-2) * ((((xy[i, 1]/th1)^(b1/dlt) + (xy[i, 2]/th2)^(b2/dlt))^(dlt)) + (1/dlt)-1) * exp(-(((xy[i, 1]/th1)^(b1/dlt) + (xy[i, 2]/th2)^(b2/dlt))^(dlt)));
		}	
		out <- sum(log(prob));
		return out;
	}
}

data {
	int LENGTH;
	matrix[LENGTH, 2] xy;
}

parameters {
	real<lower=0> theta1;
    	real<lower=0> theta2;
    	real<lower=0> theta3;
    	real<lower=0> theta4;
    	real<lower=0> theta5;
}
model {
	real lambdaI;
	real lambdaII;
	real lambdaIII;
	real lambdaIV;
	real lambdaV;
	lambdaI <- 40;
	lambdaII <- 3;
	lambdaIII <- 9;
	lambdaIV <- 1;
	lambdaV <- 0.2;
    	theta1 ~ exponential(lambdaI);
    	theta2 ~ exponential(lambdaII);
    	theta3 ~ exponential(lambdaIII);
    	theta4 ~ exponential(lambdaIV);
    	theta5 ~ exponential(lambdaIV);
    	XY ~ bivarwbl(theta1, theta2, theta3, theta4, theta4);
}