function verify_likelihood
% Check ordinary and far-tail log densities without stochastic simulations.
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
old_path = path;
cleanup = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(root);
assert(abs(loglikelihood_TRN(0,0,2,0)+log(2*sqrt(2*pi))) < 1e-12);
assert(abs(loglikelihood_TRN(0,0,sqrt(2),1)+log(2)) < 1e-12);
assert(abs(loglikelihood_TRN(0,100,1,0)+5000+0.5*log(2*pi)) < 1e-10);
assert(abs(loglikelihood_TRN(0,1000,sqrt(2),1)+1000+log(2)) < 1e-10);
assert(loglikelihood_TRN(0,3,2,0) == loglikelihood_TRN(0,-3,2,0));
fprintf('Log-likelihood checks passed.\n');
end
