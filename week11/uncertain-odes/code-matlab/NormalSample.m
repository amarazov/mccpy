function X = NormalSample(mu, Sigma, N)
% draw a sample of N random vectors
% with normal distribution mean mu and covariance
% matrix Sigma.

mu = mu(:);

dim = length(mu);

if size(Sigma) ~= [dim,dim]
    fprintf('Dimentions of Sigma does not correspond to the dimention of mu');
    return
end

X = randn(dim,N);

[U, T] = schur(Sigma,'real');

X = U*(T.^0.5)*X+repmat(mu,[1,N]);

return

end