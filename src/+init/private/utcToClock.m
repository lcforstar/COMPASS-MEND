function clockTime = utcToClock(utcTime, clockBasis, utcOffset)
%UTCTOCLOCK  Convert internal UTC time back to the forcing-clock convention.
%   local_standard: clock = UTC + utc_offset hours; utc: clock = UTC.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    if strcmpi(clockBasis, 'local_standard')
        clockTime = utcTime + hours(utcOffset);
    else
        clockTime = utcTime;
    end
end
