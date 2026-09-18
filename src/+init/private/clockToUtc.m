function utcTime = clockToUtc(clockTime, clockBasis, utcOffset)
%CLOCKTOUTC  Convert a forcing-clock datetime to the model's internal UTC time.
%   clockBasis 'local_standard': UTC = clock - utc_offset hours (e.g. Beijing +8)
%   clockBasis 'utc'           : UTC = clock (global reanalysis, ERA5)
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    if strcmpi(clockBasis, 'local_standard')
        utcTime = clockTime - hours(utcOffset);
    else
        utcTime = clockTime;
    end
end
