function R_O2_sink = calcO2Sink(Snit, Smin_aerobic, NL)
%calcO2Sink  Convert nitrogen reaction rates to O2 consumption rate.
%
%   Assembles the R_O2_sink term [NL x 1] needed by solveO2Transport from
%   the nitrification and aerobic mineralisation rates already computed by
%   the STEMMUS nitrogen module.
%
%   Stoichiometry used:
%     Nitrification:          NH4+ + 2 O2 -> NO3- + H2O + 2 H+
%                             2 mol O2 per mol N  =  2/14 mol O2 per g N
%     Aerobic mineralisation: C-H-O + O2 -> CO2 + H2O
%                             1 mol O2 per mol C  =  1/12 mol O2 per g C
%                             (pass zeros if not yet coupled)
%
%   INPUTS
%   ------
%   Snit         [NL x 2]  nitrification rate from CalculateNitrificationorginal,
%                          ug N / (cm3_bulk * s)   (negative = consumption of NH4)
%   Smin_aerobic [NL x 2]  aerobic mineralisation rate, ug C / (cm3_bulk * s)
%                          Pass zeros(NL,2) if not yet available.
%   NL           scalar    number of soil layers
%
%   OUTPUT
%   ------
%   R_O2_sink  [NL x 1]  O2 consumption rate, mol O2 / (cm3_bulk * s)
%                        (always >= 0)
%
%   UNIT CONVERSIONS applied internally:
%     ug/s  ->  g/s  :  * 1e-6 
%     g N   ->  mol O2:  / 14 * 2  (nitrification)
%     g C   ->  mol O2:  / 12 * 1  (mineralisation)

% conversion factors
ug_per_d_to_g_per_s = 1e-6 ;

% element-average rates (mean over two nodes, same as Theta_LLBAR convention)
Snit_avg  = (abs(Snit(:,1))  + abs(Snit(:,2)))  * 0.5;   % ug N/(cm3*s)
Smin_avg  = (abs(Smin_aerobic(:,1)) + abs(Smin_aerobic(:,2))) * 0.5;

% convert to mol O2/(cm3*s)
R_O2_nit  = Snit_avg * ug_per_d_to_g_per_s / 14 * 2;   % nitrification O2 demand
R_O2_min  = Smin_avg * ug_per_d_to_g_per_s / 12 * 1;   % mineralisation O2 demand

R_O2_sink = max(R_O2_nit + R_O2_min, 0);

end
