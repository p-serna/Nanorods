function xavg = running_mean(x, N)
cumsumt = cumsum(cat(1,0,x));
xavg = (cumsumt(N+1:end) - cumsumt(1:end-N)) / N;
end