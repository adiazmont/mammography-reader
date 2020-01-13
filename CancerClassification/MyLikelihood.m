function L = MyLikelihood(X,k,W,M,V)
% Compute L based on K. V. Mardia, "Multivariate Analysis", Academic Press, 1979, PP. 96-97
% to enchance computational speed
L = Likelihood(X,k,W,M,V);
return;

log_expectation = 0;
for i=1:length(X)
    gmm_one_gaussian = 0;
    gmm_one_gaussian_tmp = 0;
    for j=1:k
        gmm_one_gaussian_tmp = normpdf(X(i),M(j),V(j));
        gmm_one_gaussian = gmm_one_gaussian + W(j)*gmm_one_gaussian_tmp;
    end;
    log_expectation = log_expectation+log(gmm_one_gaussian);
end;

% divide by Total number of samples according to mit paper http://www.ll.mit.edu/IST/pubs/000101_Reynolds.pdf
log_expectation = log_expectation/length(X);
L=log_expectation;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Likelihood %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%