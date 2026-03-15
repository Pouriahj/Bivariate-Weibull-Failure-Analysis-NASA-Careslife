clc;
clear all;
close all;
%%
tic;

Data=xlsread('newdata.xlsx');

d1 = Data(:, 1);
d2 = Data(:, 2);
dnom = length(d1);

T = 5000;                                          % "T" exhibits the Iteration porcess
threshold = 15;                                 % defining threshold value

t=1;                                                   % define a temp value of "t" for using further in the loop process
theta1 = 40;                                      % defining the initial value of "thet1" parameter
theta2 = 3;                                        % defining the initial value of "thet2" parameter
theta3 = 9;                                        % defining the initial value of "Beta1" parameter
theta4 = 1;                                        % defining the initial value of "Beta2" parameter
theta5 = 0.2;                                     % defining the initial value of "delta" parameter

lambda1 = 20;                                   % Lambda1 exhibits the parameter of the exponential dist for "thet1" ("thet1" proposal dist)
lambda2 = 1.6;                                  % Lambda2 exhibits the parameter of the exponential dist for "thet2" ("thet2" proposal dist) 
lambda3 = 4.5;                                  % Lambda3 exhibits the parameter of the exponential dist for "Beta1" ("Beta1" proposal dist)                     
lambda4 = 0.5;                                  % Lambda4 exhibits the parameter of the exponential dist for "Beta2" ("Beta2" proposal dist) 
lambda5 = 0.07;                                % Lambda5 exhibits the parameter of the exponential dist for "delta" ("delta" proposal dist) 

tau1 = 5;                                            % tau1 exhibits the 2nd parameter of the gamma dist for "thet1" ('thet1" proposal dist) 
tau2 = 5;                                            % tau2 exhibits the 2nd parameter of the gamma dist for "thet2" ('thet2" proposal dist) 
tau3 = 5;                                            % tau3 exhibits the 2nd parameter of the gamma dist for "thet3" ('thet3" proposal dist) 
tau4 = 5;
tau5 = 5;

rng('default');                                    % It is suggested to define the "initial seed" in randomness processes, here I used the default seed value of the matlab 

% We defining a THETA MATRIX for our parameter with "row" = Iteration process and "column" = parameters size(3)

Theta = zeros(5, T);                          

% For bring up that weather a parameter accpeted or not we define a DECISIONS MATRIX which the "num 1" show that it is accpeted and "num 0" shows that it is rejected
% NOTE: every row might has its own particular "0 and 1 arrays" arrangement 'cause it is used "component wise accpet - reject process"

decisions = zeros(5, T);                      

Theta(1, t) = theta1;                           % leave the inital value of "thet1" in the matrix THETA
Theta(2, t) = theta2;                           % leave the inital value of "thet2" in the matrix THETA                  
Theta(3, t) = theta3;                           % leave the inital value "Beta1" in the matrix THETA
Theta(4, t) = theta3;                           % leave the inital value "Beta2" in the matrix THETA
Theta(5, t) = theta3;                           % leave the inital value "delta" in the matrix THETA

decisions(1, t) = 1;                              % it is assumed that the initial vlaue is accpeted so we leave num 1 for 1st column of the DECISION MATRIX
decisions(2, t) = 1;
decisions(3, t) = 1;
decisions(4, t) = 1;
decisions(5, t) = 1;

int1 = [10, 70];
int2 = [3, 20];

Freq = [];

% a WHILE loop till the temporary value of "t" reaches iteration process "T"

while t<T
    
    t=t+1;
    
    % "thet1" generated value via its proposal gamma dist
    % "theta1_star" stands for generated value of "thet1" - NEW ONE
    
    theta1_star = gamrnd(theta1*tau1, 1/tau1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUDGEMENT PROCESS for NEW generated "sigma0" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1st ratio is the probability (PDF) of the OLD "sigma0" to NEW "sigma0"
    
    asymmetric_ratio = gampdf(theta1, theta1_star*tau1, 1/tau1)/gampdf(theta1_star, theta1*tau1, 1/tau1);   
    
    % 2nd ratio is the probability (PDF) of the prior dist of the NEW "sigma0" to OLD "sigma0"
    
    ratio1 = exppdf(theta1_star, lambda1)/exppdf(theta1, lambda1);
    
    % 3rd ratio is the probability (PDF) of the likelihood dist (Modified Weibull) of the NEW "sigma0" to OLD "sigma0"
    % NOTE1: ratio2 is vector
    % NOTE2: "mod_wblpdf" is user defined function
    
    ratio2 = mod_wblpdf2(d1, d2, theta1_star, theta2, theta3, theta4, theta5)./mod_wblpdf2(d1, d2, theta1, theta2, theta3, theta4, theta5);
    
    % Multiplying ratios
    
    posterior_ratio = prod(ratio2).*ratio1;
    
    % calculating the alpha value
    
    alpha = min(1, posterior_ratio*asymmetric_ratio);
    
    % choosing a random num between 0 and 1
    
    u =rand;
    
    % accpet - reject of "sigma0"
    
    if alpha>u
        theta1 = theta1_star;
        decisions(1, t) = 1;
    end
    
    % "Beta" generated value via its proposal gamma dist
    % "theta2_star" stands for generated value of "Beta" - NEW ONE
       
    theta2_star = gamrnd(theta2*tau2, 1/tau2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUDGEMENT PROCESS for NEW generated"BETA" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1st ratio is the probability (PDF) of the OLD "Beta" to NEW "Beta"
    
    asymmetric_ratio = gampdf(theta2, theta2_star*tau2, 1/tau2)/gampdf(theta2_star, theta2*tau2, 1/tau2);
    
    % 2nd ratio is the probability (PDF) of the prior dist of the NEW "Beta" to OLD "Beta"
    
    ratio1 = exppdf(theta2_star, lambda2)/exppdf(theta2, lambda2);
    
    % 3rd ratio is the probability (PDF) of the likelihood dist (Modified Weibull) of the NEW "Beta" to OLD "Beta"
        
    ratio2 = mod_wblpdf2(d1, d2, theta1, theta2_star, theta3, theta4, theta5)./mod_wblpdf2(d1, d2, theta1, theta2, theta3, theta4, theta5);
    
    % Multiplying ratios
    
    posterior_ratio = prod(ratio2).*ratio1;
    
    % calculating the alpha value
     
    alpha = min(1, posterior_ratio*asymmetric_ratio);
    
    % choosing a random num between 0 and 1
    
    u =rand;
    
    % calculating the alpha value
    
    if alpha>u
        theta2 = theta2_star;
        decisions(2, t) = 1;
    end
    
    % "m" generated value via its proposal gamma dist
    % "theta3_star" stands for generated value of "m" - NEW ONE
    
    theta3_star = gamrnd(theta3*tau3, 1/tau3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUDGEMENT PROCESS for NEW generated "m" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1st ratio is the probability (PDF) of the OLD "m" to NEW "m"
    
    asymmetric_ratio = gampdf(theta3, theta3_star*tau3, 1/tau3)/gampdf(theta3_star, theta3*tau3, 1/tau3);
    
    % 2nd ratio is the probability (PDF) of the prior dist of the NEW "m" to OLD "m"
    
    ratio1 = exppdf(theta3_star, lambda3)/exppdf(theta3, lambda3);
    
    % 3rd ratio is the probability (PDF) of the likelihood dist (Modified Weibull) of the NEW "m" to OLD "m"
    
    ratio2 = mod_wblpdf2(d1, d2, theta1, theta2, theta3_star, theta4, theta5)./mod_wblpdf2(d1, d2, theta1, theta2, theta3, theta4, theta5);
    
    % Multiplying ratios
    
    posterior_ratio = prod(ratio2).*ratio1;
    
    % calculating the alpha value
    
    alpha = min(1, posterior_ratio*asymmetric_ratio);
    
    % choosing a random num between 0 and 1
    
    u =rand;
    
    % accpet - reject of "m"
    
    if alpha>u
        theta3 = theta3_star;
        decisions(3, t) = 1;
    end
    
    
    % "Beta" generated value via its proposal gamma dist
    % "theta2_star" stands for generated value of "Beta" - NEW ONE
    
    theta4_star = gamrnd(theta4*tau4, 1/tau4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUDGEMENT PROCESS for NEW generated"BETA" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1st ratio is the probability (PDF) of the OLD "Beta" to NEW "Beta"
    
    asymmetric_ratio = gampdf(theta4, theta4_star*tau4, 1/tau4)/gampdf(theta4_star, theta4*tau4, 1/tau4);
    
    % 2nd ratio is the probability (PDF) of the prior dist of the NEW "Beta" to OLD "Beta"
    
    ratio1 = exppdf(theta4_star, lambda4)/exppdf(theta4, lambda4);
    
    % 3rd ratio is the probability (PDF) of the likelihood dist (Modified Weibull) of the NEW "Beta" to OLD "Beta"
        
    ratio2 = mod_wblpdf2(d1, d2, theta1, theta2, theta3, theta4_star, theta5)./mod_wblpdf2(d1, d2, theta1, theta2, theta3, theta4, theta5);
    
    % Multiplying ratios
    
    posterior_ratio = prod(ratio2).*ratio1;
    
    % calculating the alpha value
     
    alpha = min(1, posterior_ratio*asymmetric_ratio);
    
    % choosing a random num between 0 and 1
    
    u =rand;
    
    % calculating the alpha value
    
    if alpha>u
        theta4 = theta4_star;
        decisions(4, t) = 1;
    end
    
    % "m" generated value via its proposal gamma dist
    % "theta3_star" stands for generated value of "m" - NEW ONE
    
    theta5_star = gamrnd(theta5*tau5, 1/tau5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUDGEMENT PROCESS for NEW generated "m" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1st ratio is the probability (PDF) of the OLD "m" to NEW "m"
    
    asymmetric_ratio = gampdf(theta5, theta5_star*tau5, 1/tau5)/gampdf(theta5_star, theta5*tau5, 1/tau5);
    
    % 2nd ratio is the probability (PDF) of the prior dist of the NEW "m" to OLD "m"
    
    ratio1 = exppdf(theta5_star, lambda5)/exppdf(theta5, lambda5);
    
    % 3rd ratio is the probability (PDF) of the likelihood dist (Modified Weibull) of the NEW "m" to OLD "m"
    
    ratio2 = mod_wblpdf2(d1, d2, theta1, theta2, theta3, theta4, theta5_star)./mod_wblpdf2(d1, d2, theta1, theta2, theta3, theta4, theta5);
    
    % Multiplying ratios
    
    posterior_ratio = prod(ratio2).*ratio1;
    
    % calculating the alpha value
    
    alpha = min(1, posterior_ratio*asymmetric_ratio);
    
    % choosing a random num between 0 and 1
    
    u =rand;
    
    % accpet - reject of "m"
    
    if alpha>u
        theta5 = theta5_star;
        decisions(5, t) = 1;
    end
        
    % Leave the NEW value of the "sigma0", "Beta", "m" (weather it is accepted or rejected) in the THETA matrix
    
    Theta(1, t) = theta1;
    Theta(2, t) = theta2;
    Theta(3, t) = theta3;
    Theta(4, t) = theta4;
    Theta(5, t) = theta5;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ABC Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generating new strength value with parameters (Y) 
    % "my_webiull_dist" is a user-defined function  
    
    func = @(x, y) (theta3./theta1) .* ((x./theta1).^((theta3./theta5)-1)) .* (theta4./theta2) .* ((y./theta2).^((theta4./theta5)-1)) .* ((((x./theta1).^(theta3./theta5)) + ((y./theta2).^(theta4./theta5))).^(theta5-2)) .* (((((x./theta1).^(theta3./theta5)) + ((y./theta2).^(theta4./theta5))).^(theta5)) + (1./theta5)-1) .* exp(-(((x./theta1).^(theta3./theta5))+((y./theta2).^(theta4./theta5))).^theta5);
    [Xgen, Ygen] = sample_2D(func, [int1, int2], dnom);  
    % sorting strength value for every specfic length and volume from small to large one 
    
    % calculating the norm of the differnece vector of the real strength (Strength) and generated strength (Y)
    % compare the norm of the differnece vector to threshold and accpe - reject it
    
    norm_vec = norm([sort(Xgen), sort(Ygen)] - [sort(d1), sort(d2)]);
    if norm_vec<threshold
        Freq = [Freq, [theta1; theta2; theta3; theta4; theta5]];
        sprintf("***************************** we passed step %s****************************", num2str(t))
    end
end

toc;
%% Save the FREQ Matrix
Freqq = Freq';
save('ooutFreq.mat', 'Freqq')
% labels = {'theta1'; 'theta2'; 'theta3'; 'theta4'; 'theta5'}; % labels for each row
% labeled_Freq = [labels' num2cell(Freq)]; % combine labels with the original matrix

%% Data - Visualization
% for 5k run
% Load the .csv file into Matlab as a table
table_from_fit = readtable("matrix_from_fit.csv");

% Convert the table to a matrix
matrix_from_fit = table2array(table_from_fit);
%%
% Generate a kernel density estimate of the data
[fone1, xone1] = ksdensity(Freq(1, :)');
[fone2, xone2] = ksdensity(matrix_from_fit(:, 2));

[ftwo1, xtwo1] = ksdensity(Freq(3, :)');
[ftwo2, xtwo2] = ksdensity(matrix_from_fit(:, 4));

[fthree1, xthree1] = ksdensity(Freq(2, :)');
[fthree2, xthree2] = ksdensity(matrix_from_fit(:, 3));

[ffour1, xfour1] = ksdensity(Freq(4, :)');
[ffour2, xfour2] = ksdensity(matrix_from_fit(:, 5));

[ffive1, xfive1] = ksdensity(Freq(5, :)');
[ffive2, xfive2] = ksdensity(matrix_from_fit(:, 6));

% Plot the data and the two distributions
figure
subplot('Position', [LEFT(1, 1), BOTTOM(1, 1), subplot_width, subplot_height]);
plot(xone1, fone1, 'color', '#0d32f2', 'LineWidth', 1.6); hold on;
plot(xone2, fone2, 'color', '#f60997', 'linewidth', 1.6);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 13);
ylabel("Probability Density", 'fontname', 'times', 'fontsize', 13);
set(gca, 'fontname', 'times', 'fontsize', 14)

subplot('Position', [LEFT(1, 2), BOTTOM(1, 2), subplot_width, subplot_height]);
plot(xtwo1, ftwo1, 'color', '#0d32f2', 'linewidth', 1.6); hold on;
plot(xtwo2, ftwo2, 'color', '#f60997', 'linewidth', 1.6);
xlabel("{\beta}_{1}", 'fontname', 'times', 'fontsize', 13);
ylabel("Probability Density", 'fontname', 'times', 'fontsize', 13);
set(gca, 'fontname', 'times', 'fontsize', 14)

subplot('Position', [LEFT(1, 3), BOTTOM(1, 3), subplot_width, subplot_height]);
plot(xthree1, fthree1, 'color', '#0d32f2', 'linewidth', 1.6); hold on;
plot(xthree2, fthree2, 'color', '#f60997', 'linewidth', 1.6);
xlabel("{\theta}_{2}", 'fontname', 'times', 'fontsize', 13);
ylabel("Probability Density", 'fontname', 'times', 'fontsize', 13);
set(gca, 'fontname', 'times', 'fontsize', 14)

subplot('Position', [LEFT(2, 1), BOTTOM(2, 1), subplot_width, subplot_height]);
plot(xfour1, ffour1, 'color', '#0d32f2', 'linewidth', 1.6); hold on;
plot(xfour2, ffour2, 'color', '#f60997', 'linewidth', 1.6);
xlabel("{\beta}_{2}", 'fontname', 'times', 'fontsize', 13);
ylabel("Probability Density", 'fontname', 'times', 'fontsize', 13);
set(gca, 'fontname', 'times', 'fontsize', 14)

subplot('Position', [LEFT(2, 2), BOTTOM(2, 2), subplot_width, subplot_height]);
plot(xfive1, ffive1, 'color', '#0d32f2', 'linewidth', 1.6); hold on;
plot(xfive2, ffive2, 'color', '#f60997', 'linewidth', 1.6);
xlabel("{\delta}", 'fontname', 'times', 'fontsize', 13);
ylabel("Probability Density", 'fontname', 'times', 'fontsize', 13);
set(gca, 'fontname', 'times', 'fontsize', 14)

%bubblechart(Freq(2, :), Freq(3, :), rand(1,length(Freq)))

%%
Freq = Freq';

figure
% First one
subplot(4, 5, 1)
scatter(Freq(:, 4), Freq(:, 5), 25, 'markeredgecolor', 'k', 'markerfacecolor', '#f2cd0d', 'markerfacealpha', 1);
xlabel("{\beta}_{2}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\delta}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 20)
scatter(matrix_from_fit(:, 5), matrix_from_fit(:, 6), 25, 'markeredgecolor', 'k', 'markerfacecolor', '#f2cd0d', 'markerfacealpha', 1);
xlabel("{\beta}_{2}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\delta}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)

% Second one
subplot(4, 5, 2)
scatter(Freq(:, 1), Freq(:, 2), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 1);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\theta}_{2}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 19)
scatter(matrix_from_fit(:, 2), matrix_from_fit(:, 3), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 1);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\theta}_{2}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)

subplot(4, 5, 3)
scatter(Freq(:, 1), Freq(:, 3), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 0.5);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\beta}_{1}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 18)
scatter(matrix_from_fit(:, 2), matrix_from_fit(:, 4), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 0.5);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\beta}_{1}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)

subplot(4, 5, 4)
scatter(Freq(:, 1), Freq(:, 4), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 0.15);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\beta}_{2}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 17)
scatter(matrix_from_fit(:, 2), matrix_from_fit(:, 5), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 0.15);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\beta}_{2}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)

subplot(4, 5, 5)
scatter(Freq(:, 1), Freq(:, 5), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 0.005);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\delta}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 16)
scatter(matrix_from_fit(:, 2), matrix_from_fit(:, 6), 'markeredgecolor', 'k', 'markerfacecolor', '#0d32f2', 'markerfacealpha', 0.005);
xlabel("{\theta}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\delta}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)

% Third one
subplot(4, 5, 8)
scatter(Freq(:, 2), Freq(:, 3), 'markeredgecolor', 'none', 'markerfacecolor', '#f60997', 'markerfacealpha', 1);
xlabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 13)
scatter(matrix_from_fit(:, 3), matrix_from_fit(:, 4), 'markeredgecolor', 'none', 'markerfacecolor', '#f60997', 'markerfacealpha', 1);

subplot(4, 5, 9)
scatter(Freq(:, 2), Freq(:, 4), 'markeredgecolor', 'none', 'markerfacecolor', '#f60997', 'markerfacealpha', 0.5);
xlabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 12)
scatter(matrix_from_fit(:, 3), matrix_from_fit(:, 5), 'markeredgecolor', 'none', 'markerfacecolor', '#f60997', 'markerfacealpha', 0.5);

subplot(4, 5, 10)
scatter(Freq(:, 2), Freq(:, 5), 'markeredgecolor', 'none', 'markerfacecolor', '#f60997', 'markerfacealpha', 0.1);
xlabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 11)
scatter(matrix_from_fit(:, 3), matrix_from_fit(:, 6), 'markeredgecolor', 'none', 'markerfacecolor', '#f60997', 'markerfacealpha', 0.1);

% Fourth one
subplot(4, 5, 14)
scatter(Freq(:, 3), Freq(:, 4), 'markeredgecolor', 'none', 'markerfacecolor', '#09F668', 'markerfacealpha', 1);
xlabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 7)
scatter(matrix_from_fit(:, 4), matrix_from_fit(:, 5), 'markeredgecolor', 'k', 'markerfacecolor', '#09F668', 'markerfacealpha', 1);

subplot(4, 5, 15)
scatter(Freq(:, 3), Freq(:, 5), 'markeredgecolor', 'k', 'markerfacecolor', '#09F668', 'markerfacealpha', 0.7);
xlabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
ylabel("{\theta1}_{1}", 'fontname', 'times', 'fontsize', 10);
set(gca, 'fontname', 'times', 'fontsize', 10)
subplot(4, 5, 6)
scatter(matrix_from_fit(:, 4), matrix_from_fit(:, 6), 'markeredgecolor', 'k', 'markerfacecolor', '#09F668', 'markerfacealpha', 0.7);


