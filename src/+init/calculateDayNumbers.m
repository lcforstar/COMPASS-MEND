function dateAndDayNumbers = calculateDayNumbers(start_time, end_time)
%CALCULATEDAYNUMBERS  Build the legacy day-of-year lookup table.
%
%   dateAndDayNumbers = init.calculateDayNumbers(start_time, end_time)
%   returns one row per calendar day from start_time through end_time, inclusive,
%   preserving the time-of-day carried by start_time. The table is used by
%   legacy vegetation/forcing code as JN.
%
%   Inputs:
%     start_time  [datetime] first forcing-clock day to include
%     end_time    [datetime] final forcing-clock day to include
%
%   Returns:
%     dateAndDayNumbers  table with Date [datetime] and DayNumber [1..366].
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    allDates = start_time:end_time;
    dayNumbers = zeros(length(allDates), 1);

    for i = 1:length(allDates)
        currentDate = allDates(i);
        yearVal = year(currentDate);
        monthVal = month(currentDate);
        dayVal = day(currentDate);
        dayNumbers(i) = dayOfYear(monthVal, dayVal, yearVal);
    end

    dateAndDayNumbers = table(allDates', dayNumbers, 'VariableNames', {'Date', 'DayNumber'});
end

function J = dayOfYear(M, D, year)
%DAYOFYEAR  Day number within year, matching the legacy formula.
    J = floor(275 * M / 9 - 30 + D) - 2;

    if M < 3
        J = J + 2;
    end

    if isLeapYear(year) && M > 2
        J = J + 1;
    end
end

function is_leap = isLeapYear(year)
%ISLEAPYEAR  True for Gregorian leap years.
    if mod(year, 4) == 0
        if mod(year, 100) == 0
            if mod(year, 400) == 0
                is_leap = true;
            else
                is_leap = false;
            end
        else
            is_leap = true;
        end
    else
        is_leap = false;
    end
end
