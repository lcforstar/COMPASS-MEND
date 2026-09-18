function dateAndDayNumbers = calculate_day_numbers(start_time, end_time)
%CALCULATE_DAY_NUMBERS  Legacy alias for init.calculateDayNumbers.
%
%   New initialization code should call init.calculateDayNumbers. This wrapper
%   remains for older code paths that still refer to helpers.calculate_day_numbers.
%
%   Inputs:
%     start_time  [datetime] first forcing-clock day to include
%     end_time    [datetime] final forcing-clock day to include
%
%   Returns:
%     dateAndDayNumbers  table with Date [datetime] and DayNumber [1..366].
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    dateAndDayNumbers = init.calculateDayNumbers(start_time, end_time);
end
