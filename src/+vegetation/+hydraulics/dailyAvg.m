function avg_data = dailyAvg(data, windowsize)
    n = length(data);
    n_windows = floor(n / windowsize);
    data = data(1:n_windows * windowsize);
    reshaped_data = reshape(data, windowsize, n_windows);
    avg_data = mean(reshaped_data, 1, 'omitnan')'; % Column vector
end
