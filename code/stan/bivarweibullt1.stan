functions {
	real bivarwbl_log(matrix xy, real th1, real th2, real b1, real b2, real dlt) {
        	vector[rows(xy)] prob;
		real out;
		for (i in 1:rows(xy)) {
			prob[i] <- ((b1/th1) * ((xy[i, 1]/th1)^((b1/dlt)-1))) * ((b2/th2) * ((xy[i, 2]/th2)^((b2/dlt)-1))) * ((((xy[i, 1]/th1)^(b1/dlt)) + ((xy[i, 2]/th2)^(b2/dlt)))^(dlt-2)) * (((((xy[i, 1]/th1)^(b1/dlt)) + ((xy[i, 2]/th2)^(b2/dlt)))^(dlt)) + (1/dlt)-1) * exp(-((((xy[i, 1]/th1)^(b1/dlt)) + ((xy[i, 2]/th2)^(b2/dlt)))^(dlt)));
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
	real<lower=20> theta1;
    	real<lower=0> theta2;
    	real<lower=0> theta3;
    	real<lower=0> theta4;
    	real<lower=0> theta5;
}

model {
	real alpha1;
	real beta1;
	real alpha2;
	real beta2;
	real alpha3;
	real beta3;
	real alpha4;
	real beta4;
	real alpha5;
	real beta5;

	alpha1 <- 20;
	beta1 <- 60;
	alpha2 <- 0;
	beta2 <- 20;
	alpha3 <- 0;
	beta3 <- 20;
	alpha4 <- 0;
	beta4 <- 8;
	alpha5 <- 0;
	beta5 <- 3;

    	theta1 ~ uniform(alpha1, beta1);
	theta2 ~ uniform(alpha2, beta2);
	theta3 ~ uniform(alpha3, beta3);
	theta4 ~ uniform(alpha4, beta4);
	theta5 ~ uniform(alpha5, beta5);

    	xy ~ bivarwbl(theta1, theta2, theta3, theta4, theta5);
}