%%   Arbitrage Comparison Script - Day (24) Ahead - Meeting Building loads & Exporting Excess Energy

 % REAL refers to the total energy bought from the grid to charge the
 % battery and the REAL energy exported or the net energy that is used to
 % meet building loads.
 
 % power_purchased refers to the power transfered from the grid, even when
 % the electricity prices are negative (if any).
 
 %% All data must include 8760 values for a year. Leap years not supported.
 % Data that contain more than 8760 values are NOT acceptable

%% INPUT FOR ELECTRICITY SPOT PRICES (ideally tab-delimited .txt file)
 % Both numbers that use commar or period as a decimal seperator are
 % supported. This is necessary as NordPool european data use the comma as
 % the decimal seperator.
 
 % Important: The timezone on the NordPool's N2EX hourly data must be
 % double checked as in 2016 and 2017, the price series start at 31/12
 % 23:00 of the previous year i.e. the first hour must be replaced by the
 % last price of the data to be consistent with the UK time (CET - 1).
 
 %% Load Building Energy consumption data and RTP input if necessary. Disable if variables are already loaded.
% Load building data from data_analysis_oneyear.m and RTP data from RTP_NEW_ONEYEAR.mat
% load data_analysis_oneyear.mat
% load RTP_NEW_ONEYEAR.mat

%% Input Efficiency Data of the inverter.
% An exponential fit (exp2) is created to calculate the function(P/Prated)
power_ratio_input = [0.1 0.2 0.4 0.6 0.8 1]';
pbidi_eff_input = [0.85 0.90 0.92 0.93 0.94 0.96]';
pbidi_fit = fit(power_ratio_input,pbidi_eff_input,'exp2');
pbidi_fit_coeff = coeffvalues(pbidi_fit);
figure966856 = figure('visible','off');
plot(pbidi_fit,power_ratio_input,pbidi_eff_input);
xlabel('Inverter Output Power (P/Prated)')
ylabel('Inverter Efficiency')
savefig('pbidi_fit_eff_var.fig')
saveas(gcf,'pbidi_fit_eff_var.png')

% The following line to be used for the pbidi_eff_variable of the storage
% model post-loop.
% pbidi_fit_coeff(1)*exp(pbidi_fit_coeff(2)*X) + pbidi_fit_coeff(3)*exp(pbidi_fit_coeff(4)*X)

%% Basic input Data for the Battery and the converter
bat_cap = 200; %Battery capacity in kWh
pbidi_cap_disch = -90; % inverter rated AC power (kW)
pbidi_cap_ch = 65 ; % rectifier DC power (kW)
bat_volt = 51.8; % battery voltage (V)

%% Qualitative constraints to eliminate charging/discharging in peak/inactive times
% minHourLimit & maxHourLimit not implemented for the first day of the year.
% minHourLimit calculated as the average load of the first day (bank holiday).

minHourPowerLimit_enforce = true; % necessary to avoid higher peak building loads
maxHourPowerLimit_enforce = true; % the battery will be significantly more inactive and exported energy will be minimised.
exclusions = true; % This will enable the exclusion of hours from charging or discharging based on the respective building loads in order to avoid higher peak loads and discharging during off-peak periods.
% If this option is disabled (false), there is no need for input for the Building's Activity Basic Profile below.
% Exclusions are NOT applicable to the first day of the year.
pbidi_cap_disch_PowerLimit_enforce = true; % 
minHourPowerLimit_margin = 5; % allow x kW as a margin error
maxHourPowerLimit_margin = 5; % allow x kW as a margin error
minHourPowerLimit = mean(totalenergyLIFE(1:24)) + minHourPowerLimit_margin; % set to the average value of the 1st day loads.
% It is possible to use instead the hour range during which the building is
% not operational (e.g. midnight - 7pm and 7pm-midnight).
maxHourPowerLimit = minHourPowerLimit; % equal for this occasion.

%% Building Activity Basic Profile (Alternative Strategy, not needed if exclusions == true)

% This is set primarily in DesignBuilder but it must be introduced here as
% well in order to improve the qualitative constraints. The Activity mostly
% consists of the opening hours and the preheating period.
starting_hour = 9; % meaning 08:00 in DB i.e. the time slot 8-9am. It is a convention as the first hour of the day is midnight-1am
preheating_period = 1; % how many hours before the beggining of opening hours, setpoint-based heating starts.
closing_hour = 18; % the closing time of the building e.g. for 18, the building will close at 17:00-18:00 (6pm). 
% That means that there will be no loads at the 19th hour of the day
% (6-7pm) and the last active period will be the 18th hour (5-6pm aka 6pm sharp)
% 
% In each case, the number of the starting/closing hour corresponds to the
% number of the period during the day. For example:

% Period 1: midnight - 1 am
% Period 2: 1 am - 2 am
% Period 3: 3 am - 3 am
% Period 19: 19 - 20 (pm)
% Period 24: 23 - midnight

% IMPORTANT. Occupation/heating/cooling/lighting/DHW/RoomEquipment profiles 
% to be double-checked in DesignBuilder before running the building simulation!

% The battery is not allowed to charge during the preheating hour. For
% example, if starting_hour = 9 and preheating_period = 1, the battery will
% not charge for the 8th period (7-8am) to avoid achieving a higher peak.

% If the closing hour is 18 and the building closes at 5-6pm, then the
% battery must be prohibited from discharging the following hour (i.e.
% closing hour + 1, when no loads will take place).

%% Battery and Converter Specs
delta = 0; % hourly battery self-discharge coefficient (%)
SOC_min = 0.1; % Depth of dishcarge - minimum SOC %
SOC_max = 1; % Maximum battery SOC
SOC_start = 0.1; % Initial SOC for t=1
bat_net_cap = bat_cap*(1-SOC_min); % useful battery capacity
imaxd = 750; % maximum battery discharge current (A). Optional
imaxc = 1500; % maximum battery charge current (A). Optional
nbattch = 0.97; % Battery charging efficiency
nbattd = 0.97; % Battery discharging efficiency 
pbidi_eff_RECTIFIER = 0.96; % assuming the rectifier efficiency is constant
pbidi_eff = 0.96; % If the efficiency is calculated as the mean value of the input at the beggining of the code, disable this line and enable the following one
% pbidi_eff = mean(pbidi_eff_input); % Only enable if the efficiency is calculated from the input at the beggining of the code. Disable line 109 if so.
pbidi_eff_CONSTANT = true; % switch. When true, it has a constant value equal to pbidi_eff. When false, it follows the exponential fit.
% Post-loop, the real efficincies of the inverter will be applied if the efficiency is not considered constant. This is
% possible as only the RECTIFIER efficiency is included in the bottlenecks.
roundtrip_eff = nbattch * nbattd * pbidi_eff_RECTIFIER * pbidi_eff;
pbidi_cap_ch = pbidi_cap_ch * nbattch * pbidi_eff_RECTIFIER;  % The rectifier size is reduced to take into account the losses.
Ncycles90 = 5000; % Lifetime of the battery for a DOD of 90%
Ncycles = (1-SOC_min)*Ncycles90; % number of equivalent full cycles until failure
%Ncycles = 5000; % Use ONLY if the lifecycle of the battery for the equivalent full cycles is given.

%% Miscellaneous
pmax = 99999; % max power that can be bought from the grid (kW)
% life_system = 10; % lifetime of the project in years. This is now SET in data_analysis.m instead.

%% Economic parameters & Costs
interest = 0.05; % interest rate
inflation = 0.02; % inflation rate of electricity
inflation_el = 0.02; % inflation rate of electricity 
inflation_OM = 0.02; % inflation rate of O&M
VAT = 0; % Value Added Tax

% Costs
cost_fixed_elec_annual = 0; % Fixed cost per year, paid to the Network Operator
cost_battery_exVAT = 390 * bat_cap;
cost_battery_final = cost_battery_exVAT * (1+VAT);
cost_converter_excVAT = 170 * abs(pbidi_cap_disch);  % Sunny Island 8.0H
cost_converter_final = cost_converter_excVAT * (1+ VAT); % 6 refers to the discharge capacity of the Sunny Island 8.0H
cost_OM_annual = 0; % annual O&M cost (£)
cost_capital_initial = cost_battery_final + cost_converter_final; %Initial capital cost
% END OF INPUT %

%% Energy content (Assuming the battery is initially at its minimun SOC).
initial_energy_stored = bat_cap * SOC_min;
energy_stored = initial_energy_stored * ones(length(RTP_retail),1);
energy_exchange = zeros(length(RTP_retail),1);
totalprofit = zeros(length(RTP_retail),1);
profitperkWh = totalprofit;
number_of_loops = 0;
HourIndexExamined = ones(length(RTP_retail),1); 
% 1 means that the hour has yet to be examined and removed from the series.
% 0 means that the price has been already removed (reason may vary).


%% FIRST DAY LOOP

while any(HourIndexExamined(1:24)) == 1 % while there are values in the time series not removed, perform the while loop.
number_of_loops = number_of_loops + 1; % counter of the loops

% Identify maxHourIndex in the distribution of the RTP prices and the range of minHour
NON_REMOVED_HOURS = find(HourIndexExamined(1:24) == 1); % first line of the matrix includes the non-removed hours
PRICES_OF_NON_REMOVED_HOURS = RTP_retail(HourIndexExamined(1:24) == 1); % second line of the matric includes the respective prices
COUNTER = find(PRICES_OF_NON_REMOVED_HOURS == max(PRICES_OF_NON_REMOVED_HOURS),1,'first'); % find the highest price of the day.
maxHourIndex = NON_REMOVED_HOURS(COUNTER);
maxHourPrice = PRICES_OF_NON_REMOVED_HOURS(COUNTER);
clear NON_REMOVED_HOURS
clear PRICES_OF_NON_REMOVED_HOURS
clear COUNTER;

   if exclusions == true % then specific hours will be excluded from charging or discharging based on the respective building loads       
   % this will ensure that hours will be removed from the DISCHARGING phase if during a weekday, the building load is lower than the set limit. Weekends are NOT excluded from discharging.
        if (isweekend(datestampLIFE(maxHourIndex)) == 0 && ismember(datestamp2LIFE(maxHourIndex),BankHolidays) == false  && totalenergyLIFE(maxHourIndex) <= maxHourPowerLimit && maxHourPowerLimit_enforce == 1)
        HourIndexExamined(maxHourIndex) = 0;
        continue
        end
        
        if ((isweekend(datestampLIFE(maxHourIndex)) == 1)  || (ismember(datestamp2LIFE(maxHourIndex),BankHolidays) == true))
        HourIndexExamined(maxHourIndex) = 0;
        continue
        end
   end    

% Identify the last time before maxHourIndex that the battery was fully
% charged plus one hour.
minRangeIndex = find(energy_stored(1:maxHourIndex) == bat_cap,1,'last') + 1;
if(isempty(minRangeIndex) == 1)  % If the battery is never full
    minRangeIndex = 1;           % the earliest hour it can charge is the first period in the series 
end

if minRangeIndex > maxHourIndex   % If the battery is already full during maxHourIndex
   minRangeIndex = maxHourIndex;  % Ensure minRangeIndex is not greater than maxHourIndex
end
    
% The latest the battery can charge is the hour before the battery reached its minimum SOC.
% '-1' is used to get the previous period and 'maxHourIndex -1' to convert the index found through the find function.
maxRangeIndex = find(energy_stored(maxHourIndex:24)== SOC_min * bat_cap,1,'first')+ (maxHourIndex - 1) - 1;
% modified to include SOC_min * bat_cap as the minimum battery capacity allowed.
% It is noted again that the (maxHourIndex -2) is used to convert the 'find' function to the index system.
if(isempty(maxRangeIndex) == 1) % if the battery is not found to be at SOC_min
    maxRangeIndex = 24; % the range will include all the values of the first day.
end

if maxRangeIndex < maxHourIndex    % True if the battery is already at the minimum SOC during maxHourIndex
   maxRangeIndex = maxHourIndex;   % Ensure maxRangeIndex is not less than maxHourIndex
end
% The if statement above ensures that a number of periods are excluded from the beggining/end of the time series
% where the discharging/charging capacities reached their maximum potential.
% Identify minHourIndex when charging takes place
RTP_MINRANGE_TO_MAXRANGE = RTP_retail(minRangeIndex:maxRangeIndex); % only prices included
HourIndexExamined_MINRANGE_TO_MAXRANGE = HourIndexExamined(minRangeIndex:maxRangeIndex); % hours yet to be removed in the the min-max hour price range.
% In case there are no hours left that have not been removed 
% (meaning there are no hours that meet the condition HourIndexExamined = 1), maxHourIndex is removed
if(isempty(RTP_MINRANGE_TO_MAXRANGE(HourIndexExamined_MINRANGE_TO_MAXRANGE == 1)) == 1)
    HourIndexExamined(maxHourIndex) = 0;
else
    NON_REMOVED_HOURS = find(HourIndexExamined_MINRANGE_TO_MAXRANGE == 1); % first line of the matrix contains indeces of non-removed time series
    PRICES_OF_NON_REMOVED_HOURS = RTP_MINRANGE_TO_MAXRANGE(HourIndexExamined_MINRANGE_TO_MAXRANGE == 1); % second line contains the respective prices
    COUNTER = find(PRICES_OF_NON_REMOVED_HOURS == min(PRICES_OF_NON_REMOVED_HOURS),1,'first');
    minHourIndex = NON_REMOVED_HOURS(COUNTER) + (minRangeIndex - 1); % (minRangeIndex - 1) is used to convert the index found through the 'find' function
    minHourValue = PRICES_OF_NON_REMOVED_HOURS(COUNTER);
    clear NON_REMOVED_HOURS;
    clear PRICES_OF_NON_REMOVED_HOURS;
    clear COUNTER;
        
    if exclusions == true % then specific hours will need to be excluded from charging.
        % Hours will be removed from charging if it is a weekeday and the building loads exceed the set limit at that time. Charging takes place in weekends.
            if (isweekend(datestampLIFE(minHourIndex)) == 0 && ismember(datestamp2LIFE(minHourIndex),BankHolidays) == false && totalenergyLIFE(minHourIndex) >= minHourPowerLimit && minHourPowerLimit_enforce == 1)
            HourIndexExamined(minHourIndex) = 0;
            continue
            end
        % Hours will be excluded from charging if it is during a weekend or a bank holiday 
            if ((isweekend(datestampLIFE(minHourIndex)) == 1)  || (ismember(datestamp2LIFE(minHourIndex),BankHolidays) == true))
            HourIndexExamined(minHourIndex) = 0;
            continue
            end
    end
    
    % Calculation of the marginal operating cost of the battery (Capital Costs not included)
    % Initial assumption to consider MCdischarge and MCcharge to be zero.
    MCdischarge = 0; % pounds per kWh, based on the kWh discharged during the lifetime of the battery in cycles
    % This is basically the marginal operating cost of generation (OPEX).
    MCcharge = 0; % This is the marginal operating cost (OPEX) for charging the battery storage system
    MCmargin = 0; % Unit: pounds. This works as an input in order to constrain the battery utilisation if needed.
    
    % Provision for negative electricity prices.
     if RTP_retail(minHourIndex) >= 0
           MCproduction = MCdischarge + (RTP_retail(minHourIndex)+ MCcharge)/(pbidi_eff_RECTIFIER * nbattch * pbidi_eff * nbattd);  % £/kWh of the real discharge, after all losses and four efficiencies.
     else
           MCproduction = MCdischarge + ((RTP_retail(minHourIndex) + MCcharge) * nbattd * pbidi_eff)/(pbidi_eff_RECTIFIER * nbattch);
     end
        
     if (MCproduction + MCmargin) < (RTP_retail(maxHourIndex)) && (minHourIndex ~= maxHourIndex) 
        % if it is cost-effective to use the battery and minHourIndex is not the same with maxHourIndex
        % The bottleneck will be the minimum of the three bottlebeck conditions
        bottleneck = zeros(1,3);
        bottleneck(1) = energy_exchange(maxHourIndex)- pbidi_cap_disch; % available 'turbine' (discharging) capacity at maxHourIndex
        bottleneck(2) = pbidi_cap_ch - energy_exchange(minHourIndex); % available 'pumping' (charging) capacity at minHourIndex
          if maxHourIndex > minHourIndex % charging takes place first
             storageFree = bat_cap - max(energy_stored(minHourIndex:maxHourIndex)); % minimum free storage space
             bottleneck(3) = storageFree; % This does not allow the battery to charge more than a SOC of 100%
          else % discharging takes place first
             storageLeft = min(energy_stored(maxHourIndex:minHourIndex)) - SOC_min * bat_cap; % minimum storage content
             bottleneck(3) = storageLeft; % so that the battery will not get discharged to a SOC lower than the SOC_min
          end
        final_bottleneck = min(bottleneck);

        % Updating the operation of the battery system and energy exchanges
        % that take place. Energy_stored does not have to updated for
        % minHour and maxHour, only the energy_exchange changes
        energy_exchange(minHourIndex) = energy_exchange(minHourIndex) + final_bottleneck  ; % charging, including losses
        energy_exchange(maxHourIndex) = energy_exchange(maxHourIndex) - final_bottleneck ; % discharging, including losses
        if maxHourIndex > minHourIndex
            energy_stored(minHourIndex:maxHourIndex-1) = energy_stored(minHourIndex:maxHourIndex-1) + final_bottleneck; 
            % as discharging takes place later and energy increases after charging at minHourIndex
        else % if maxHourIndex takes place prior to minHourIndex
            energy_stored(maxHourIndex:minHourIndex-1) = energy_stored(maxHourIndex:minHourIndex-1) - final_bottleneck;
        end

        % Remove time periods when the charging and discharging capacities
        % reached their full potential.
        if(energy_exchange(maxHourIndex) <= pbidi_cap_disch)
            HourIndexExamined(maxHourIndex) = 0;
        end
        if(energy_exchange(minHourIndex) >= pbidi_cap_ch)
            HourIndexExamined(minHourIndex) = 0;
        end
    else
        % If it is not cost-effective to operate the battery, remove the
        % tarriffs from the time series.
        HourIndexExamined(maxHourIndex) = 0;
        HourIndexExamined(minHourIndex) = 0;
    end
end

end

%% WHILE LOOP FOR THE REST 364 DAYS of the calendar year.

for pp = 2:length(RTP_retail)/24
    energy_stored(24*pp-24+1:24*pp) = energy_stored(24*pp-24+1-1); % This will update the storage content based on what happened the previous day.
    
    while any(HourIndexExamined(24*pp-24+1:24*pp)) == 1 % while there are values in the time series not removed, perform the while loop.

    number_of_loops = number_of_loops + 1;
    %% Identify maxHourIndex in the distribution of the RTP prices and the range of minHour
    NON_REMOVED_HOURS = find(HourIndexExamined(24*pp-24+1:24*pp)==1) + 24 * pp - 24 + 1 - 1; % first line of the matrix includes the non-removed hours
    PRICES_OF_NON_REMOVED_HOURS = RTP_retail(NON_REMOVED_HOURS); % second line of the matric includes the respective prices
    COUNTER = find(PRICES_OF_NON_REMOVED_HOURS  == max(PRICES_OF_NON_REMOVED_HOURS),1,'first'); % find the highest price of the second line
    maxHourIndex = NON_REMOVED_HOURS(COUNTER);
    maxHourPrice = PRICES_OF_NON_REMOVED_HOURS(COUNTER);
    clear NON_REMOVED_HOURS;
    clear PRICES_OF_NON_REMOVED_HOURS;
    clear COUNTER;
    
    % Constraint to remove discharging hours according to the building's
    % activity in order to avoid higher peak loads during the day and
    % discharging of the battery during the off-peak hours.
    
   if exclusions == true % then specific hours will be excluded from charging or discharging based on the respective building loads       
   % this will ensure that hours will be removed from the DISCHARGING phase if during a weekday, the building load is lower than the set limit. Weekends are NOT excluded from discharging.
        if (isweekend(datestampLIFE(maxHourIndex)) == 0 && ismember(datestamp2LIFE(maxHourIndex),BankHolidays) == false  && totalenergyLIFE(maxHourIndex) <= maxHourPowerLimit && maxHourPowerLimit_enforce == 1)
        HourIndexExamined(maxHourIndex) = 0;
        continue
        end
        
        if ((isweekend(datestampLIFE(maxHourIndex)) == 1)  || (ismember(datestamp2LIFE(maxHourIndex),BankHolidays) == true))
        HourIndexExamined(maxHourIndex) = 0;
        continue
        end
   end       

   % FROM NOW ON, the code is mostly a repetition of Day 1, with the exception of exclusions for the charging phase.
   
    % Identify the last time before maxHourIndex that the battery was fully
    % charged plus one hour.
    minRangeIndex = find(energy_stored(24*pp-24+1:maxHourIndex) == bat_cap,1,'last') + (24*pp - 24 + 1) - 1 + 1;
    if(isempty(minRangeIndex) == 1)  % If the battery is never full
        minRangeIndex = 24*pp-24+1;           % the earliest hour it can charge is the first period on the day 
    end

    if minRangeIndex > maxHourIndex   % If battery is already full during maxHour
       minRangeIndex = maxHourIndex;  % Ensure minRangeIndex is not greater than maxHourIndex
    end

    % The latest the battery can charge is the previous period from when the
    % battery is at minimum SOC. '-1' is used to get the previous period and 'maxHourIndex -1' to
    % convert the index found through the find function.
    maxRangeIndex = find(energy_stored(maxHourIndex:24*pp)== SOC_min * bat_cap,1,'first')+ (maxHourIndex - 1) - 1; % modified to include SOC_min * bat_cap instead of 0 kWh
    if(isempty(maxRangeIndex) == 1) % if the battery is not found to be at SOC_min
        maxRangeIndex = 24*pp; % the range will include all the values of the time series
    end

    if maxRangeIndex < maxHourIndex      % If the battery is already at SOC_min during maxHour
       maxRangeIndex = maxHourIndex;     % Ensure maxRangeIndex is not less than maxHourIndex
    end

    %% Identify minHourIndex when charging takes place
    RTP_MINRANGE_TO_MAXRANGE = RTP_retail(minRangeIndex:maxRangeIndex); % only prices included
    HourIndexExamined_MINRANGE_TO_MAXRANGE = HourIndexExamined(minRangeIndex:maxRangeIndex); % hours yet to be removed in the the min-max hour price range.

    % In case there are no hours left that have not been removed 
    % (meaning there are no hours that meet the condition HourIndexExamined = 1), maxHourIndex is removed
    if(isempty(RTP_MINRANGE_TO_MAXRANGE(HourIndexExamined_MINRANGE_TO_MAXRANGE == 1)) == 1)
        HourIndexExamined(maxHourIndex) = 0;
    else
        NON_REMOVED_HOURS = find(HourIndexExamined_MINRANGE_TO_MAXRANGE == 1) + (minRangeIndex - 1); % first line of the matrix contains indeces of non-removed time series
        PRICES_OF_NON_REMOVED_HOURS = RTP_MINRANGE_TO_MAXRANGE(HourIndexExamined_MINRANGE_TO_MAXRANGE == 1); % second line contains the respective prices
        COUNTER = find(PRICES_OF_NON_REMOVED_HOURS == min(PRICES_OF_NON_REMOVED_HOURS),1,'first');
        minHourIndex =  NON_REMOVED_HOURS(COUNTER);  % is used to convert the index found through the 'find' function
        minHourValue = PRICES_OF_NON_REMOVED_HOURS(COUNTER);
        clear NON_REMOVED_HOURS;
        clear PRICES_OF_NON_REMOVED_HOURS;
        clear COUNTER;

        % Constraint to remove minHourIndeces that will induce higher peak
        % loads than the original ones.
   
        if exclusions == true % then specific hours will need to be excluded from charging.
        % Hours will be removed from charging if it is a weekeday and the building loads exceed the set limit at that time. Charging takes place in weekends.
            if (isweekend(datestampLIFE(minHourIndex)) == 0 && ismember(datestamp2LIFE(minHourIndex),BankHolidays) == false && totalenergyLIFE(minHourIndex) >= minHourPowerLimit && minHourPowerLimit_enforce == 1)
            HourIndexExamined(minHourIndex) = 0;
            continue
            end
        % Hours will be excluded from charging if it is during a weekend or a bank holiday 
            if ((isweekend(datestampLIFE(minHourIndex)) == 1)  || (ismember(datestamp2LIFE(minHourIndex),BankHolidays) == true))
            HourIndexExamined(minHourIndex) = 0;
            continue
            end
        end
       
        % Calculation of the marginal operating cost of the battery (Capital Costs not included)
        % Initial assumption to consider MCdischarge and MCcharge to be
        % zero but they can be changed if deemed necessary.
        MCdischarge = 0; % pounds per kWh, based on the kWh discharged during the lifetime of the battery in cycles.
        % This is the OPEX (marginal cost of generation for the battery
        % storage system in pounds per kWh.
        MCcharge = 0; % the marginal cost (OPEX) of 'pumping' (charging) the battery storage system.
        MCmargin = 0; % Unit: pounds. This works as an input in order to constrain the battery utilisation.
        
        % Provision for negative electricity prices.
        if RTP_retail(minHourIndex) >= 0
           MCproduction = MCdischarge + (RTP_retail(minHourIndex)+ MCcharge)/(pbidi_eff_RECTIFIER * nbattch * pbidi_eff * nbattd);  % £/kWh of the real discharge, after all losses and four efficiencies.
        else
           MCproduction = MCdischarge + ((RTP_retail(minHourIndex) + MCcharge) * nbattd * pbidi_eff)/(pbidi_eff_RECTIFIER * nbattch);
        end
        
        if (MCproduction + MCmargin) < (RTP_retail(maxHourIndex)) && (minHourIndex ~= maxHourIndex) 
            % if it is cost-effective to use the battery
            % The bottleneck will be the minimum of the three bottlebeck conditions
            bottleneck = zeros(1,3);
            bottleneck(1) = energy_exchange(maxHourIndex)- pbidi_cap_disch ;
            bottleneck(2) = pbidi_cap_ch - energy_exchange(minHourIndex);
              if maxHourIndex > minHourIndex % charging takes place first 
                 storageFree = (bat_cap - max(energy_stored(minHourIndex:maxHourIndex)));
                 bottleneck(3) = storageFree; % This does not allow the battery to charge more than a SOC of 100%
              else % discharging takes place first
                 storageLeft = min(energy_stored(maxHourIndex:minHourIndex)) - SOC_min * bat_cap;
                 bottleneck(3) = storageLeft; % so that the battery will not get discharged to a SOC lower than the SOC_min
              end
            final_bottleneck = min(bottleneck);

            % Updating the operation of the battery system and energy exchanges
            % that take place. Energy_stored does not have to updated for
            % minHour and maxHour, only the energy_exchange changes
            energy_exchange(minHourIndex) = energy_exchange(minHourIndex) + final_bottleneck ;
            energy_exchange(maxHourIndex) = energy_exchange(maxHourIndex) - final_bottleneck;
            if maxHourIndex > minHourIndex
                energy_stored(minHourIndex:maxHourIndex-1) = energy_stored(minHourIndex:maxHourIndex-1) + final_bottleneck; % efficiencies are included as when charging, fewer energy is stored into the battery
                 % as discharging takes place later and energy increases after charging at minHourIndex
            else % if maxHourIndex takes place prior to minHourIndex
                energy_stored(maxHourIndex:minHourIndex-1) = energy_stored(maxHourIndex:minHourIndex-1) - final_bottleneck; % efficiencies not used after disharging as if pbidi_cap_disch = -10, then 10 kWh will be lost.
            end

            % Remove time periods when the charging and discharging capacities
            % reached their full potential.
            if (energy_exchange(maxHourIndex) <= pbidi_cap_disch)
                HourIndexExamined(maxHourIndex) = 0;
            end

            if(energy_exchange(minHourIndex) >= pbidi_cap_ch) 
                HourIndexExamined(minHourIndex) = 0;
            end
      else
            % If it is not cost-effective to operate the battery, remove the
            % tarriffs from the time series.
            HourIndexExamined(maxHourIndex) = 0;
            HourIndexExamined(minHourIndex) = 0;
        end
    end

    end
end

%% POST-LOOP CALCULATIONS

% The REAL refers to the real (net) energy that can be UTILISED By the
% building when discharging and the gross energy used to charge the battery
% (grid purchases)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

%% Scenario 1: Using the Battery

% Calculate the real inverter efficiencies based on the exponential fit of
% the inverter manufacturer sheet.
% pbidi_fit_coeff(1)*exp(pbidi_fit_coeff(2)*X) + pbidi_fit_coeff(3)*exp(pbidi_fit_coeff(4)*X)
pbidi_eff_variable = zeros(length(RTP_retail),1);
power_ratio = zeros(length(RTP_retail),1);
power_ratio(energy_exchange >=0) = 0;
power_ratio(energy_exchange <0) = energy_exchange(energy_exchange <0)/pbidi_cap_disch;
pbidi_eff_variable(energy_exchange >=0) = 0;

if  pbidi_eff_CONSTANT == false % then efficiency will change every time the battery discharges.
    pbidi_eff_variable(energy_exchange <0) =  pbidi_fit_coeff(1)*exp(pbidi_fit_coeff(2)*power_ratio(energy_exchange <0)) + pbidi_fit_coeff(3)*exp(pbidi_fit_coeff(4)*power_ratio(energy_exchange <0));
else % pbidi_eff_variable will have the value of pbidi_eff
    pbidi_eff_variable(energy_exchange <0) = pbidi_eff;
end

energy_exchange_REAL = zeros(length(energy_exchange),1);
energy_exchange_REAL(energy_exchange < 0) = energy_exchange(energy_exchange < 0) .* pbidi_eff_variable(energy_exchange < 0) * nbattd; 
% incorporating the real variable inverter efficiency
energy_exchange_REAL(energy_exchange > 0) = energy_exchange(energy_exchange > 0) / (pbidi_eff_RECTIFIER * nbattch);


% energy_stored already refers to the net energy stored that will be used
% when charging and discharging. Therefore, there is no need to calculate
% the real one.

SOC = energy_stored/bat_cap;
battery_charge = zeros(length(RTP_retail),1); battery_discharge = zeros(length(RTP_retail),1);
battery_charge(energy_exchange > 0) =  energy_exchange(energy_exchange >  0);
battery_discharge(energy_exchange < 0) = - energy_exchange(energy_exchange < 0);
battery_charge_REAL = zeros(length(RTP_retail),1);
battery_charge_REAL(energy_exchange_REAL > 0) = energy_exchange_REAL(energy_exchange_REAL > 0); % Including losses
battery_discharge_REAL = zeros(length(RTP_retail),1);
battery_discharge_REAL(energy_exchange_REAL < 0) = - energy_exchange_REAL(energy_exchange_REAL < 0);

% Current calculation, assuming battery voltage is constant.
current = energy_exchange*1000/bat_volt; % optional
current_REAL = energy_exchange_REAL*1000/bat_volt; % optional 

% Seperating battery discharge to both cover building loads and export the
% remaining energy content.

remaining_loads = zeros(length(RTP_retail),1);
exported_energy = zeros(length(RTP_retail),1);

% Remaining loads refer ONLY to the building's remaining loads when
% discharging the battery. For example, if the battery charges at 18 kW,
% this is not included in the remaining loads. Therefore, the total
% power_purchased and energy_shifted are calculated below.

for ee = 1:length(RTP_retail)
if battery_discharge_REAL(ee) >= totalenergyLIFE(ee)
    remaining_loads(ee) = 0;
    exported_energy(ee) = battery_discharge_REAL(ee) - totalenergyLIFE(ee);
else
    remaining_loads(ee) = totalenergyLIFE(ee) - battery_discharge_REAL(ee);
    exported_energy(ee) = 0;
end
end

    battery_charge_annual = battery_charge(1:8760); 
    battery_discharge_annual = battery_discharge(1:8760); 
    battery_charge_REAL_annual = battery_charge_REAL(1:8760); 
    battery_discharge_REAL_annual = battery_discharge_REAL(1:8760);
    battery_output_annual_scalar = sum(battery_discharge(1:8760));
    exported_energy_annual = exported_energy(1:8760);
    remaining_loads_annual = remaining_loads(1:8760);
    exported_energy_annual_scalar = sum(exported_energy(1:8760));
    remaining_loads_annual_scalar = sum(remaining_loads(1:8760));


battery_lifetime = (bat_cap * Ncycles)/battery_output_annual_scalar;
% Ncycles refers to the number of full equivalent cycles.
% This will constitute an estimation of the battery lifetime, based on its
% 1st year of of operation.

% Calculation of electricity costs purchased, O&M & revenues
cost_payg_WITH_matrix = zeros(length(RTP_retail),1);
energy_from_the_grid = zeros(length(RTP_retail),1);

for tt = 1:length(RTP_retail)
    
if energy_exchange_REAL(tt) >= 0 && RTP_retail(tt) >=0
    cost_payg_WITH_matrix(tt) = RTP_retail(tt) *(totalenergyLIFE(tt) + energy_exchange_REAL(tt)); % energy_exchange_REAL is needed as it refers to the gross grid purchases.
    energy_from_the_grid(tt) = totalenergyLIFE(tt) + energy_exchange_REAL(tt);
elseif energy_exchange_REAL(tt) >= 0 && RTP_retail(tt) < 0
    cost_payg_WITH_matrix(tt) = 0; % no cost as there is a negative spot electricity price.
    energy_from_the_grid(tt) = totalenergyLIFE(tt) + energy_exchange_REAL(tt); % power bought for free considered.
elseif energy_exchange_REAL(tt) < 0 
% It is assumed that the battery will only discharge during a positive RTP price
    cost_payg_WITH_matrix(tt) = RTP_retail(tt) * remaining_loads(tt); % removing loads met by the battery
    energy_from_the_grid(tt) = remaining_loads(tt);
end

end
cost_payg_WITH_scalar = sum(cost_payg_WITH_matrix);

revenues_WITH_matrix = zeros(length(RTP_retail),1);
% to include revenues for charging the battery (and meeting building loads) on negative electricity prices
% and for any energy exported back to the grid.
% It is assumed that the building is awarded the final retail electricity
% price in order to export to the grid during the times of expensive
% wholesale prices.
for tt = 1:length(RTP_retail) % Revenues are based on the wholesale energy prices, NOT the retail electricity prices.
    if RTP_retail(tt) < 0
           revenues_WITH_matrix(tt) = - RTP_retail(tt) .* (energy_from_the_grid(tt)); % assuming that discharge is not possible during negative RTP prices
    elseif RTP_retail(tt) >=0
           revenues_WITH_matrix(tt) = RTP_retail(tt) .* exported_energy(tt); 
    end        
end
revenues_WITH_scalar = sum(revenues_WITH_matrix);

%% Alternative Case - Exports are awarded with the WHOLESALE electricity price only
% This does not change the technical part of the algorithm but it is used
% purely to demonstrate the difference between retail and wholesale revenues.

revenues_WITH_alternative_matrix = zeros(length(RTP_retail),1);
for tt = 1:length(RTP_retail)
    if RTP_retail(tt) < 0
           revenues_WITH_alternative_matrix(tt) = - RTP_wholesale_kWh(tt) .* (energy_from_the_grid(tt)); % assuming that discharge is not possible during negative RTP prices
    elseif RTP_retail(tt) >=0
           revenues_WITH_alternative_matrix(tt) = RTP_wholesale_kWh(tt) .* exported_energy(tt); 
    end        
end
revenues_WITH_alternative_scalar = sum(revenues_WITH_alternative_matrix);


%% Basic economic parameters for the 1st year
energy_from_the_grid_annual_matrix = energy_from_the_grid(1:8760); % energy bought
energy_from_the_grid_annual_scalar = sum(energy_from_the_grid_annual_matrix);
cost_payg_WITH_annual_matrix = cost_payg_WITH_matrix(1:8760);
cost_payg_WITH_annual_scalar = sum(cost_payg_WITH_annual_matrix); % NOT ADJUSTED FOR INFLATION/INTEREST RATE
revenues_WITH_annual_matrix = revenues_WITH_matrix(1:8760);
revenues_WITH_annual_scalar = sum(revenues_WITH_annual_matrix); % NOT ADJUSTED FOR INFLATION/INTEREST RATE

% Calculating the energy shifted  (meeting building loads)
power_shifted = zeros(length(RTP_retail),1);

for ee = 1:length(RTP_retail)
    if energy_from_the_grid(ee) < totalenergyLIFE(ee) % When the power purchased is less than than original demand
       power_shifted(ee) = totalenergyLIFE(ee) - energy_from_the_grid(ee);
    else 
        power_shifted(ee) = 0; 
        % power_shifted does not take into account the energy shift when charging the battery
    end
end

power_shifted_annual_matrix = power_shifted(1:8760);
power_shifted_annual_scalar = sum(power_shifted_annual_matrix);

%% Scenario 2 BAU: Not Using the battery
cost_payg_WITHOUT_matrix = zeros(length(RTP_retail),1);
cost_payg_WITHOUT_matrix(RTP_retail >=0 ) = RTP_retail(RTP_retail >=0) .* totalenergyLIFE(RTP_retail >=0);
cost_payg_WITHOUT_matrix(RTP_retail < 0) = 0;
cost_payg_WITHOUT_matrix_annual = cost_payg_WITHOUT_matrix(1:8760); % for the 1st year
cost_payg_WITHOUT_scalar_annual = sum(cost_payg_WITHOUT_matrix_annual);


% Revenues from purchasing electricity with negative prices (if any) are
% included in the WITHOUT BATTERY scenario. Retail price considered
revenues_WITHOUT_matrix = zeros(length(RTP_retail),1);
revenues_WITHOUT_matrix(RTP_retail >= 0) = 0;
revenues_WITHOUT_matrix(RTP_retail < 0) = RTP_retail(RTP_retail<0) .* totalenergyLIFE(RTP_retail<0);
revenues_WITHOUT_scalar = sum(revenues_WITHOUT_matrix);


%% Fix the datestampLIFE. This can happen only post-loop.
% This is purely symbolic.
datestampLIFE = datestamp;

 
%% NEW FIGURES WITHOUT LIMITS FOR VALIDATION PURPOSES

figure905531 = figure('visible','off');
yyaxis right; xlabel('Date'); plot(datestampLIFE,energy_stored); ylabel('Energy Stored (kWh)'); yyaxis left; plot(datestampLIFE,RTP_retail); ylabel('RTP Price (£/kWh)');
savefig('PHES_RTP_Arbitrage_8760_stored_energy.fig')
saveas(gcf,'PHES_RTP_Arbitrage_8760_stored_energy.png')

figure9455555931 = figure('visible','off');
yyaxis right; xlabel('Date'); plot(datestampLIFE,SOC*100); ylabel('SOC (%)'); yyaxis left; plot(datestampLIFE,RTP_retail); ylabel('RTP Price (£/kWh)');
 savefig('PHES_RTP_Arbitrage_8760_SOC.fig')
saveas(gcf,'PHES_RTP_Arbitrage_8760_SOC.png')

figure5956254596 = figure('visible','off');
yyaxis right; xlabel('Date'); plot(datestampLIFE,energy_exchange_REAL,'-.o','MarkerSize', 7); ylabel('Battery Charge (+) and Battery Discharge (-) (kW)'); yyaxis left; plot(datestampLIFE,RTP_retail); ylabel('RTP Price (£/kWh)')
xlim auto; ylim auto
savefig('PHES_RTP_Arbitrage_8760_battery_power.fig')
saveas(gcf,'PHESRTPArbitrage8760batterypower2.png')
%%
figure14156244496 = figure('visible','on');
yyaxis right; xlabel('Date'); plot(datestampLIFE,totalenergyLIFE,'k-.'); hold on; plot(datestampLIFE, - exported_energy, 'color', [0.0000 0.4900 0.0000]); hold on; plot(datestampLIFE, -battery_discharge_REAL, 'r', datestampLIFE, battery_charge_REAL,'r', datestampLIFE, energy_from_the_grid,'k'); 
ylabel('Loads (kW)');  yyaxis left; plot(datestampLIFE,RTP_retail,'-.'); ylabel('RTP Price (£/kWh)');% ylim([0.05 inf]) 
lgd = legend('Real-time Price','Original Loads','Exported Energy','Battery Charge/Discharge', '.', 'Loads using Storage','Location', 'north','Orientation', 'horizontal'); lgd.FontSize = 11;
savefig('PHES_RTP_Arbitrage_8760_overall.fig') 
saveas(gcf,'PHES_RTP_Arbitrage_8760_overall.png')
fprintf('Total energy dicharged by the battery (excluding losses): %g kWh\nTotal energy exported to the grid: %g kWh\nTotal energy discharged by the battery to meet building loads: %g kWh\n',sum(battery_discharge),sum(exported_energy),sum(totalenergyLIFE) - sum(remaining_loads));
fprintf('The lifetime of the battery is %g years based on the number of the full equivalent cycles.\n',battery_lifetime)
clear lgd
%%
figuret541419654558 = figure('visible','off');
yyaxis right; xlabel('Date'); plot(datestampLIFE,exported_energy,'LineStyle','none','Marker','o','MarkerSize',10,'LineWidth',1); ylabel('Exported Energy (kW)'); yyaxis left; plot(datestampLIFE,RTP_retail,'-.'); ylim([0.03 0.20]); ylabel('RTP Price (£/kWh)')
savefig('PHES_RTP_Arbitrage_8760_battery_power.fig')
saveas(gcf,'PHES_RTP_Arbitrage_8760_battery_power.png')

%% Clear unwanted variables
clear rr pp tt ee raw

%% Further post-processing for the Annual Comparison Script
totalenergy_WITHOUT_permonth = zeros(12,1);
totalcost_WITHOUT_permonth = zeros(12,1);
revenues__WITHOUT_permonth = zeros(12,1);
totalenergy_WITH_permonth = zeros(12,1);
totalcost_WITH_permonth = zeros(12,1);
revenues_WITH_permonth = zeros(12,1);
for ii = 1:12
% Monthly Variables without battery storage
totalenergy_WITHOUT_permonth(ii) = sum(totalenergy(month(datestamp)==ii)); % building loads without battery storage
totalcost_WITHOUT_permonth(ii) = sum(cost_payg_WITHOUT_matrix(month(datestamp)==ii)); % energy cost without battery storage
revenues__WITHOUT_permonth(ii) = sum(revenues_WITHOUT_matrix(month(datestamp) == ii)); % revenuews without battery storage (if negative electricity prices are present)
% Monthly Variables with battery storage
totalenergy_WITH_permonth(ii) = sum(energy_from_the_grid(month(datestamp)==ii)); % building loads without battery storage
totalcost_WITH_permonth(ii) = sum(cost_payg_WITH_matrix(month(datestamp)==ii)); % energy cost without battery storage
revenues_WITH_permonth(ii) = sum(revenues_WITH_matrix(month(datestamp) == ii)); % revenues without battery storage (if negative electricity prices are present)
end
clear ii

%% Total annual energy loads and final costs
totalenergy_WITHOUT_annual = sum(totalenergy_WITHOUT_permonth);
totalcost_WITHOUT_annual = sum(totalcost_WITHOUT_permonth);
revenues_WITHOUT_annual = revenues_WITHOUT_scalar;
net_cost_WITHOUT_annual = totalcost_WITHOUT_annual - revenues_WITHOUT_annual; % assuming retail prices
totalenergy_WITH_annual = sum(totalenergy_WITH_permonth);
totalcost_WITH_annual = sum(totalcost_WITH_permonth);
exports_annual = sum(exported_energy);
revenues_WITH_annual = sum(revenues_WITH_permonth); % assuming retail prices
net_cost_WITH_annual = totalcost_WITH_annual - revenues_WITH_annual; % assuming retail prices
net_cost_WITH_alternative_annual = totalcost_WITH_annual - revenues_WITH_alternative_scalar; % assuming wholesale prices.

%% Calculate the NPC and LCOE with and without Battery Storage.

% Asigning the life duration of the project based on the battery lifetime

life_system = 10; % Set to 10 years

% Calculation of NPC/LCOE without Battery Storage
cost_payg_WITHOUT_annual_IG = zeros(life_system,1);
revenues_WITHOUT_annual_IG = zeros(life_system,1);
    
 for rr = 1:life_system
 cost_payg_WITHOUT_annual_IG(rr) = cost_payg_WITHOUT_scalar_annual * ((1 + inflation_el)/(1 + interest))^rr;
 revenues_WITHOUT_annual_IG(rr) = revenues_WITHOUT_scalar * ((1 + inflation_el)/(1 + interest))^rr;
 end;
  
cost_payg_WITHOUT_lifetime = sum(cost_payg_WITHOUT_annual_IG);
revenues_WITHOUT_lifetime = sum(revenues_WITHOUT_annual_IG);
    
% Calculation with Battery Storage

cost_payg_WITH_annual_IG = zeros(life_system,1); % cost of the electricity 'fuel' per year
revenues_WITH_annual_IG = zeros(life_system,1); % revenues per year
revenues_WITH_alternative_IG = zeros(life_system,1); % revenues per year, considering the wholesale price.
cost_OM_annual_IG = zeros(life_system,1);  % O&M costs per year
cost_fixed_elec_annual_IG = zeros(life_system,1); % fixed electricity costs per year, if any.
       
for rr = 1:life_system
    cost_payg_WITH_annual_IG(rr) = cost_payg_WITH_annual_scalar * ((1 + inflation_el)/(1 + interest))^rr;
    revenues_WITH_annual_IG(rr) = revenues_WITH_annual_scalar * ((1 + inflation_el)/(1 + interest))^rr;
    revenues_WITH_alternative_IG(rr) = revenues_WITH_alternative_scalar * ((1 + inflation_el)/(1 + interest))^rr;
    cost_fixed_elec_annual_IG(rr) = cost_fixed_elec_annual * ((1 + inflation_el)/(1 + interest))^rr;
    cost_OM_annual_IG(rr) = cost_OM_annual * ((1 + inflation_el)/(1 + interest))^rr;
end;
      
cost_payg_WITH_lifetime = sum(cost_payg_WITH_annual_IG);
revenues_WITH_lifetime = sum(revenues_WITH_annual_IG);
revenues_WITH_alternative_lifetime = sum(revenues_WITH_alternative_IG);
cost_OM_lifetime = sum(cost_OM_annual_IG);
cost_fixed_elec_lifetime = sum(cost_fixed_elec_annual_IG);
clear rr

%% NPC and LCOE

NPC_with_storage = cost_battery_final + cost_converter_final + cost_OM_lifetime + cost_fixed_elec_lifetime...
    + cost_payg_WITH_lifetime - revenues_WITH_lifetime;

NPC_with_storage_alternative = cost_battery_final + cost_converter_final + cost_OM_lifetime + cost_fixed_elec_lifetime...
    + cost_payg_WITH_lifetime - revenues_WITH_alternative_lifetime;

NPC_without_storage  = cost_fixed_elec_lifetime + cost_payg_WITHOUT_lifetime - revenues_WITHOUT_lifetime;


% Calculating LCOE. The total NPC will be divided by the number of loads
% i.e. sum(totalenergy(1:8760)) for one year.

LCOE_with_storage = NPC_with_storage/((sum(totalenergy)+sum(exported_energy))*life_system); % considering retail prices as revenues.
LCOE_with_storage_alternative = NPC_with_storage_alternative/((sum(totalenergy)+sum(exported_energy))*life_system); % considering wholesale prices
LCOE_without_storage = NPC_without_storage/(sum(totalenergyLIFE)*life_system);

financial_motive_shifted = (NPC_with_storage - NPC_without_storage)/(sum(power_shifted)*life_system);
financial_motive_shifted_exports = (NPC_with_storage - NPC_without_storage)/((sum(power_shifted) + sum(exported_energy)) * life_system) ;
financial_motive_shifted_alt = (NPC_with_storage_alternative - NPC_without_storage)/(sum(power_shifted)*life_system);
financial_motive_shifted_exports_alt = (NPC_with_storage_alternative - NPC_without_storage)/((sum(power_shifted) + sum(exported_energy)) * life_system) ;
financial_motive_exports = (NPC_with_storage - NPC_without_storage)/(sum(exported_energy) * life_system) ;
financial_motive_exports_alt = (NPC_with_storage_alternative - NPC_without_storage)/(sum(exported_energy) * life_system) ;

%% Print the final comparison values.
fprintf('Energy Consumption (kWh) Energy Cost (£) Exports (kWh) Revenues (£) Net Costs (£) Battery Usable Discharges (kWh) Energy Shifted \n')
fprintf('With Battery Storage: %.2f kWh £%.2f %.2f kWh £%.2f £%.2f %.2f kWh %.2f kWh \n',totalenergy_WITH_annual,totalcost_WITH_annual,exports_annual,revenues_WITH_annual,net_cost_WITH_annual,sum(battery_discharge_REAL),power_shifted_annual_scalar);
fprintf('Without Battery Storage: %.2f kWh £%.2f 0 kWh £%.2f £%.2f 0 kWh 0 kWh \n',totalenergy_WITHOUT_annual,totalcost_WITHOUT_annual,revenues_WITHOUT_annual,net_cost_WITHOUT_annual);
fprintf('When considering the exporting scenario and wholesale prices ONLY, the revenues are: £%.2f and the net cost is £%.2f \n',revenues_WITH_alternative_scalar,net_cost_WITH_alternative_annual);
fprintf('When considering the exporting scenario and final retail prices, the revenues are: £%.2f and the net cost is £%.2f \n',revenues_WITH_scalar,net_cost_WITH_annual);

aaa = [totalenergy_WITH_annual,totalcost_WITH_annual,exports_annual,revenues_WITH_annual,net_cost_WITH_annual,sum(battery_discharge_REAL),power_shifted_annual_scalar,battery_lifetime,(power_shifted_annual_scalar/sum(totalenergy))*100,(power_shifted_annual_scalar/peak_loads_scalar)*100, 0.02, net_cost_WITH_alternative_annual/2500, NPC_without_storage, NPC_with_storage, NPC_with_storage_alternative, LCOE_without_storage, LCOE_with_storage, LCOE_with_storage_alternative, financial_motive_shifted, financial_motive_shifted_exports, financial_motive_shifted_alt, financial_motive_shifted_exports_alt]';
%%  Save all the results (Data_Analysis/RTP/Arbitrage)
fname3 = sprintf('ARBITRAGE_NPC_COMPARISON_1b_%.0f_%.0f_%.0f_%s.mat',bat_cap,pbidi_cap_ch/(nbattch*pbidi_eff_RECTIFIER),-pbidi_cap_disch,ScenarioID);
save(fname3)
save ARBITRAGE_NPC_COMPARISON.mat
clear fname3