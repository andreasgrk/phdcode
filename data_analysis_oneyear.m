%% Analysis of Building data for One year, imported from DesignBuilder to Microsoft Excel.

%% Building & Location Characteristics for the results filename (.mat) (Part 1)
Shape = 'Rectangular';
ScenarioID = '2rN_BIR_2017_55BP';
Year = 2017; % Input for the calendar year
% It should be pointed out that the analysis below refers to the building
% loads WITHOUT taking into account any additional loads that may take
% place when the battery storage system is operational and extra
% electricity will be needed to charge the battery.

% Mean hourly profiles per month are calculated by exluding weekends (by
% using the isweekend function). This assumes that the building has
% typical commercial building opening hours, which are set in DesignBuilder
% through the activity tab.

% The building is active only on weekdays (8am - 6pm) with one hour of
% preheating; therefore, energy activity can be seen as early as 7am,
% always depending on the weather conditions. This can be easily changed
% through the activity schedule in the current model.

% The building does not operate during the UK bank holidays. Due to the difficulties
% recognising their respective dates, as Matlab only has embedded 
% functions for US public holidays, they are included to the calculations of
% mean monthly profiles but do not affect the overall results (only 8 days
% per calendar year).

% Due to the different number of days per month, it was not possible to use
% a single matrix with 12 columns; therefore, seperate matrices for each
% month were used.

% The first day of the year can and MUST be set in DesignBuilder simulation
% options, under the Location Tab and the Hourly Weather File Options.
% While the results and the dates are correct, the days are reported incorrectly inside DB 
% and the north american syntax is used (MM/dd/yyyy).
% Therefore, an option to override the DB days was added and a manual datestamp 
% can be created by using the datetime function and the specified year.
% 2017 starts on a Sunday and 2015 on a Thursday. 2016 is avoided as it's
% a leap year.

% Bank Holidays must also be configured inside DB, from the Activity Tab as
% there is no way to import them into Matlab or export them from Matlab to
% DB. If this does not happen, there will be incosistencies as the building
% loads will be matched with the wrong dates (e.g. full loads on a bank
% holiday). BankHolidays must be set in MATLAB as well.

% Daylight Saving Time (DST) MUST BE turned off, at the Location Options,
% inside DB, otherwise there are going to be inconsistencies regarding the
% energy profile of the building.

% This script is used for comparison cases between different buildings and
% therefore includes only data from ONE calendar year. However, as the
% storage model still uses lifetime matrices, the same names will be used
% for the final variables but their values will be set equal to the annual
% ones (e.g. totalenergyLIFE = totalenergy(1:8760).

%% Import the data, extracting spreadsheet dates in Excel serial date format
xls_building_data_filename = input('Enter the name for the xlsx file: ', 's');
[~, ~, raw] = xlsread(xls_building_data_filename,'','B2:G8761',''); % change Excel cells location if needed
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
raw = raw(:,[1,2,3,4,5,6]);

%% INPUT FOR THE YEAR

% As DesignBuilder does not support leap years, only common years can be
% examined. The first day of the year (Monday etc.) can be set from the weather file
% settings in DesignBuilder while the year number from the variable below.

% YEAR has been set at the beggining of this script

% life_system = 10; % not used as this script is for annual comparisons
% only. Calculated at the end of the storage model based on the battery
% lifetime
duplicate = false; % not used as this script is for annual comparisons only

%% Create datetime table for bank holidays depending on the year.
% This covers only the years 2015 and 2017 (2016 is a leap year and therefore avoided)
if Year == 2015
    BankHolidays = transpose([datetime(Year,1,1) datetime(Year,4,3) datetime(Year,4,6) datetime(Year,5,4) datetime(Year,5,25) datetime(Year,8,31) datetime(Year,12,25) datetime(Year,12,28)]);
elseif Year == 2017 
    BankHolidays =  transpose([datetime(Year,1,2) datetime(Year,4,10) datetime(Year,12,25) datetime(Year,12,26)]);
elseif Year == 2018 
    BankHolidays =  transpose([datetime(Year,1,1) datetime(Year,4,10) datetime(Year,12,25) datetime(Year,12,26)]);
else 
    error('The number of the Year is neither 2015/2017 nor 2018. Check again.')
end
BankHolidays.Format = 'dd-MM-yyyy 00:00';  

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create temporary output variable
data = reshape([raw{:}],size(raw));

%% IMPORTANT: OVERRIDE DESIGNBUILDER DATESTAMP AND CREATE A NEW ONE
% To make sure that the datestamp used is correct.
override_designbuilder_datestamp = true;

if override_designbuilder_datestamp == true
    tt1 = datetime(Year,1,1,0,0,0);
    tt2 = datetime(Year,12,31,23,0,0);
    datestamp = transpose(tt1:hours(1):tt2);
    datestamp.Format = 'eee dd/MM/yyyy HH:00';
    clear tt1 tt2
end

% Create a datestamp2 vector without hours or minutes (using the midnight
% hour of the day) in order to be able to directly compare datetimes later
% in the storage model.

datestamp2 = dateshift(datestamp,'start','day','current');
datestamp2.Format = 'dd/MM/yyyy HH:00';

%% Identify the datetimes when DST is observed and the respectice indeces.
% Identify the date and time when DST starts and ends (last Sunday of March
% and October. 1 am (1-2am) for both occasions.)
dst1 = calendar(Year,3); dst2 = calendar(Year,10);
last_sunday_dst1 = dst1(dst1(:,1)~= 0,1); last_sunday_dst2 = dst2(dst2(:,1)~= 0,1); % because the function 'calendar', the first day of the week is Sunday.
last_sunday_dst1 = last_sunday_dst1(end);
last_sunday_dst2 = last_sunday_dst2(end);
datetime_dst1 = datetime(Year,03,last_sunday_dst1,01,00,00);
datetime_dst2 = datetime(Year,10,last_sunday_dst2,01,00,00);

% Identify the datetime_pp1 index in the datestamp
index_dst1 = find(datestamp == datetime_dst1); % index 2018 for the year 2017
index_dst2 = find(datestamp == datetime_dst2);% index 7226 for year 2017 
clear dst1 dst2

%% Continue with assigning energy sectors to variables.
roomElectricity = data(:,1); % refers to equipment loads
lighting = data(:,2);
auxiliary = data(:,3); % parasitic energy (fans/pumps/controls)
heating = data(:,4);

% Set cooling equal to zero if the building is naturally ventilated
if ~isnan(data(:,6))
    cooling = data(:,5) ;
    DHW = data(:,6);
    fprintf('\nThe Building is  mechanically ventilated and extra cooling is provided through the Heat Pump.\n')
elseif isnan(data(:,6))
cooling = zeros(8760,1);
DHW = data(:,5);
fprintf('\nThe Building is naturally ventilated and therefore cooling loads are zero.\n')
end

%% Clear temporary variables
clearvars R;

%% Building & Location Characteristics for the results filename (.mat) (Part 2)
if cooling == zeros(8760,1)
    Ventilation = 'natural';
else
    Ventilation = 'cooling';
end
%% Set the days included in every month (leap years not supported but the code is written anyway)

% Not necessary as only common years will be examined. Keeping for potential future
% use.

days_per_month = [31; 0; 31; 30; 31; 30; 31; 31; 30; 31; 30; 31];
if numel(datestamp(month(datestamp) == 2)) == 672 % for common years
    days_per_month(2) = 28;
    Year_type = 'common';
elseif numel(datestamp(month(datestamp == 2))) == 696 % for leap years
    days_per_month(2) = 29;
    Year_type = 'leap';
end

%% Create the sum of all the energy consumption types and save it to file (ptional)
totalenergy = auxiliary + cooling + DHW + heating + lighting + roomElectricity;
% baseFileName = sprintf('totalenergy_%s.txt', xls_building_data_filename); % Enable these three lines if needed
% dlmwrite(baseFileName,totalenergy,'delimiter',';','precision',4);
% clear baseFileName 

%% Reshape the total energy table to 24 lines and 365 columns
total_energyR = reshape(totalenergy, [24 365]);
datestampR = reshape(datestamp, [24 365]);

%% Create a column vector with the respective months and its respective datetime vector
monthstamp = month(datestamp);
months = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September' 'October', 'November', 'December'}';
monthsno = datetime(months,'InputFormat','MM');
monthsno.Format = 'MMMM';

%% Create DAYSTAMP for one entire year
t1 = datetime(Year,1,1,0,0,0); % Year is set at the input of this file
t2 = datetime(Year,12,31,23,0,0); % Year is set at the input of this file
daystamp = t1:caldays(1):t2;
daystamp = daystamp';
daystamp.Format = 'eee dd-MMM-yyyy';
clear t1 t2

%% Plot total energy against datestamp
figure1 = figure('visible','off');
plot(datestamp,totalenergy, 'LineWidth',1);
grid on
title('Total Annual Electricity Consumption')
xlabel('time'); ylabel('kWh')
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
savefig('total_annual_energy.fig')
saveas(gcf,'total_annual_energy.png')

%% Calculate the total energy consumption per sector for the entire first  year (kWh)
total_AUX = sum(auxiliary); total_COOL = sum(cooling); total_DHW = sum(DHW); total_HEAT = sum(heating); 
total_LIGHT = sum(lighting); total_ROOM = sum(roomElectricity); TOTAL = sum(totalenergy);

%% Calculate the energy consumption of each sector for each month (kWh)
% Pre-allocating matrices size for speed

month_AUX = zeros(12,1); month_COOL = zeros(12,1); month_HEAT = zeros(12,1); month_LIGHT = zeros(12,1);
month_TOTAL = zeros(12,1); month_ROOM = zeros(12,1); month_DHW = zeros(12,1);

% Continuing with the calculation
for ii=1:12
    month_AUX(ii,1) = sum(auxiliary(monthstamp == ii));
    month_COOL(ii,1) = sum(cooling(monthstamp == ii));
    month_LIGHT(ii,1) = sum(lighting(monthstamp == ii));
    month_ROOM(ii,1) = sum(roomElectricity(monthstamp == ii));
    month_TOTAL(ii,1) = sum(totalenergy(monthstamp == ii));
    month_DHW(ii,1) = sum(DHW(monthstamp == ii));
    month_HEAT(ii,1) = sum(heating(monthstamp == ii));
end
clear ii

%% Plot the sector energy consumptions per month and total energy consumption per month
figure2 = figure('visible','off');
plot(monthsno, month_AUX, monthsno, month_COOL, monthsno, month_DHW, monthsno, month_HEAT, monthsno, month_LIGHT, monthsno, month_ROOM,'LineWidth',3)
title('Monthly Electricity consumption per sector')
grid on
xlabel('month'); ylabel('kWh') ;
legend('Auxiliary', 'Cooling','DHW', 'Heating', 'Lighting', 'Room Electricity');
set(legend,'Location','best','Orientation','horizontal');
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
savefig('monthy_energy_per_sector.fig')
saveas(gcf,'monthly_energy_per_sector.png')

figure3 = figure('visible','off');
plot(monthsno, month_TOTAL,'LineWidth',3)
grid on
title('Total Monthly Electricity consumption')
xlabel('month'); ylabel('kWh') ;
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
savefig('monthly_energy_total.fig')
saveas(gcf,'monthly_energy_total.png')

%% Create a bar chart with monthly energy consumption per sector and percentages

% Calculate the respective percentages
AUX_percentage = (total_AUX/TOTAL)*100; COOL_percentage = (total_COOL/TOTAL)*100;
DHW_percentage = (total_DHW/TOTAL)*100; HEAT_percentage = (total_HEAT/TOTAL)*100;
LIGHT_percentage = (total_LIGHT/TOTAL)*100; ROOM_percentage = (total_ROOM/TOTAL)*100;
SECTOR_percentages = [AUX_percentage COOL_percentage DHW_percentage HEAT_percentage LIGHT_percentage ROOM_percentage 100];
sectors = { 'Auxiliary', 'Cooling', 'DHW', 'Heating', 'Lighting', 'Equipment', 'Total'};

figure4 = figure('visible','off');
bar(1:length(sectors), [total_AUX total_COOL total_DHW total_HEAT total_LIGHT total_ROOM TOTAL]);
barvalues;
title('Annual Electricity Consumption per sector')
ylabel('kWh')
xlabel('Sector')
xticklabels(sectors)
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
savefig('annual_energy_bar.fig')
saveas(gcf,'annual_energy_bar.png')
clear tt

figure85hj = figure('visible','off');
bar(1:length(sectors), SECTOR_percentages);
barvalues;
title('Annual Electr. Consumption per sector (%)')
xlabel('Sector')
xticklabels(sectors)
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
savefig('annual_energy_bar_percent.fig')
saveas(gcf,'annual_energy_bar_percent.png')

%% Calculate and plot mean daily electricity profile (from all 12 months).

temp = zeros(1,24);
for ii=1:24
 temp(ii) = (mean(totalenergy(ii:24:end)));
 daily_mean_total = transpose(temp);
end
clear temp
clear ii

figure5 = figure('visible','off');
plot(0:23, daily_mean_total,'LineWidth',3)
grid on
title('Mean Daily Electricity Profile')
xlabel('Hour of the day'); ylabel('kWh') ;
xlim([0 23]);
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xticks(0:2:23)
savefig('daily_mean_total.fig')
saveas(gcf,'daily_mean_total.png')

%% Identify the most intensive month in terms of total energy spent, sort its loads and the respective datestamp
intensive_month = monthsno(month_TOTAL == max(month_TOTAL));
fprintf('The most energy intensive month is %s with a total energy consumption of %.2f kWh.\n', intensive_month, max(month_TOTAL))
intensive_month_load = (totalenergy(month(datestamp) == month(intensive_month)));
days_of_peak_month = numel(intensive_month_load)/24;
t1 = datetime(Year,month(intensive_month),1,0,0,0);
t2 = datetime(Year,month(intensive_month),days_of_peak_month,23,0,0);
intensive_monthstamp = transpose(t1:hours(1):t2);
intensive_monthstamp.Format = 'eee dd/MM/yyyy HH:mm';
clear t1 t2

%% Plot the most intensive month in terms of total energy spent
figure6 = figure('visible','off');
plot(intensive_monthstamp, intensive_month_load,'LineWidth',3)
grid on
title(sprintf('Total loads of the most intensive month (%s)',intensive_month))
xlabel('Date'); ylabel('kWh') ;
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
savefig('intensive_month_loads.fig')
saveas(gcf,'intensive_month_loads.png')


%% Identify the datestamp and energy consumption for each month, per hour.
january_datestamp = datestamp(month(datestamp)==1);
february_datestamp = datestamp(month(datestamp)==2);
march_datestamp = datestamp(month(datestamp)==3);
april_datestamp = datestamp(month(datestamp)==4);
may_datestamp = datestamp(month(datestamp)==5);
june_datestamp = datestamp(month(datestamp)==6);
july_datestamp = datestamp(month(datestamp)==7);
august_datestamp = datestamp(month(datestamp)==8);
september_datestamp = datestamp(month(datestamp)==9);
october_datestamp = datestamp(month(datestamp)==10);
november_datestamp = datestamp(month(datestamp)==11);
december_datestamp = datestamp(month(datestamp)==12);

january_month_load = totalenergy(month(datestamp)==1);
february_month_load = totalenergy(month(datestamp)==2);
march_month_load = totalenergy(month(datestamp)==3);
april_month_load = totalenergy(month(datestamp)==4);
may_month_load = totalenergy(month(datestamp)==5);
june_month_load = totalenergy(month(datestamp)==6);
july_month_load = totalenergy(month(datestamp)==7);
august_month_load = totalenergy(month(datestamp)==8);
september_month_load = totalenergy(month(datestamp)==9);
october_month_load = totalenergy(month(datestamp)==10);
november_month_load = totalenergy(month(datestamp)==11);
december_month_load = totalenergy(month(datestamp)==12);

%% Calculate the mean WORKING hourly profile of each month (excluding weekends)
datestamp_weekdays = datestamp(~isweekend(datestamp)); % identify the dates of working days (weekdays)
datestamp_weekdaysR = reshape(datestamp_weekdays,[24 numel(datestamp_weekdays)/24]);
weekdays_number = day(datestamp_weekdaysR,'dayofyear');
weekdays_number = transpose(weekdays_number(1,:)); % identify the weekends day number in a year (e.g. 6,7,13,14 etc.)
weekdays_hourly_indicesR = zeros(24, numel(datestamp_weekdays)/24); % hourly index in a year (1 to 8760 possible)

for trt = 1:length(weekdays_number)
    weekdays_hourly_indicesR(:,trt) = (weekdays_number(trt)-1)*24 +1 :24*weekdays_number(trt);
end
clear trt

weekdays_hourly_indices = reshape(weekdays_hourly_indicesR, [numel(weekdays_hourly_indicesR) 1]);

% Setting 12 temporary matrices to identify logically if the month when the
% weekdays take place is January, February etc.

koukou1 =  month(datestamp(weekdays_hourly_indices)) == 1;
koukou2 =  month(datestamp(weekdays_hourly_indices)) == 2;
koukou3 =  month(datestamp(weekdays_hourly_indices)) == 3;
koukou4 =  month(datestamp(weekdays_hourly_indices)) == 4;
koukou5 =  month(datestamp(weekdays_hourly_indices)) == 5;
koukou6 =  month(datestamp(weekdays_hourly_indices)) == 6;
koukou7 =  month(datestamp(weekdays_hourly_indices)) == 7;
koukou8 =  month(datestamp(weekdays_hourly_indices)) == 8;
koukou9 =  month(datestamp(weekdays_hourly_indices)) == 9;
koukou10 = month(datestamp(weekdays_hourly_indices)) == 10;
koukou11 = month(datestamp(weekdays_hourly_indices)) == 11;
koukou12 = month(datestamp(weekdays_hourly_indices)) == 12;

% Total energy for the weekdays per month
totalenergy_weekdays_january = totalenergy(weekdays_hourly_indices(koukou1 == 1));
totalenergy_weekdays_february = totalenergy(weekdays_hourly_indices(koukou2 == 1));
totalenergy_weekdays_march = totalenergy(weekdays_hourly_indices(koukou3 == 1));
totalenergy_weekdays_april = totalenergy(weekdays_hourly_indices(koukou4 == 1));
totalenergy_weekdays_may = totalenergy(weekdays_hourly_indices(koukou5 == 1));
totalenergy_weekdays_june = totalenergy(weekdays_hourly_indices(koukou6 == 1));
totalenergy_weekdays_july = totalenergy(weekdays_hourly_indices(koukou7 == 1));
totalenergy_weekdays_august = totalenergy(weekdays_hourly_indices(koukou8 == 1));
totalenergy_weekdays_september = totalenergy(weekdays_hourly_indices(koukou9 == 1));
totalenergy_weekdays_october = totalenergy(weekdays_hourly_indices(koukou10 == 1));
totalenergy_weekdays_november = totalenergy(weekdays_hourly_indices(koukou11 == 1));
totalenergy_weekdays_december = totalenergy(weekdays_hourly_indices(koukou12 == 1));

clear koukou1 koukou2 koukou3 koukou4 koukou5 koukou6 koukou7 koukou8 koukou9 koukou10 koukou11 koukou12

% Total energy for the weekdays per month RESHAPED
totalenergy_weekdays_januaryR = reshape(totalenergy_weekdays_january,[24 numel(totalenergy_weekdays_january)/24]);
totalenergy_weekdays_februaryR = reshape(totalenergy_weekdays_february,[24 numel(totalenergy_weekdays_february)/24]);
totalenergy_weekdays_marchR = reshape(totalenergy_weekdays_march,[24 numel(totalenergy_weekdays_march)/24]);
totalenergy_weekdays_aprilR = reshape(totalenergy_weekdays_april,[24 numel(totalenergy_weekdays_april)/24]);
totalenergy_weekdays_mayR = reshape(totalenergy_weekdays_may,[24 numel(totalenergy_weekdays_may)/24]);
totalenergy_weekdays_juneR = reshape(totalenergy_weekdays_june,[24 numel(totalenergy_weekdays_june)/24]);
totalenergy_weekdays_julyR = reshape(totalenergy_weekdays_july,[24 numel(totalenergy_weekdays_july)/24]);
totalenergy_weekdays_augustR = reshape(totalenergy_weekdays_august,[24 numel(totalenergy_weekdays_august)/24]);
totalenergy_weekdays_septemberR = reshape(totalenergy_weekdays_september,[24 numel(totalenergy_weekdays_september)/24]);
totalenergy_weekdays_octoberR = reshape(totalenergy_weekdays_october,[24 numel(totalenergy_weekdays_october)/24]);
totalenergy_weekdays_novemberR = reshape(totalenergy_weekdays_november,[24 numel(totalenergy_weekdays_november)/24]);
totalenergy_weekdays_decemberR = reshape(totalenergy_weekdays_december,[24 numel(totalenergy_weekdays_december)/24]);

mean_hourly_weekdays_january_profile = mean(totalenergy_weekdays_januaryR,2);
mean_hourly_weekdays_february_profile = mean(totalenergy_weekdays_februaryR,2);
mean_hourly_weekdays_march_profile = mean(totalenergy_weekdays_marchR,2);
mean_hourly_weekdays_april_profile = mean(totalenergy_weekdays_aprilR,2);
mean_hourly_weekdays_may_profile = mean(totalenergy_weekdays_mayR,2);
mean_hourly_weekdays_june_profile = mean(totalenergy_weekdays_juneR,2);
mean_hourly_weekdays_july_profile = mean(totalenergy_weekdays_julyR,2);
mean_hourly_weekdays_august_profile = mean(totalenergy_weekdays_augustR,2);
mean_hourly_weekdays_september_profile = mean(totalenergy_weekdays_septemberR,2);
mean_hourly_weekdays_october_profile = mean(totalenergy_weekdays_octoberR,2);
mean_hourly_weekdays_november_profile = mean(totalenergy_weekdays_novemberR,2);
mean_hourly_weekdays_december_profile = mean(totalenergy_weekdays_decemberR,2);

mean_hourly_weekdays_ALLMONTHS_profile = [mean_hourly_weekdays_january_profile mean_hourly_weekdays_february_profile mean_hourly_weekdays_march_profile mean_hourly_weekdays_april_profile mean_hourly_weekdays_may_profile mean_hourly_weekdays_june_profile mean_hourly_weekdays_july_profile mean_hourly_weekdays_august_profile mean_hourly_weekdays_september_profile mean_hourly_weekdays_october_profile mean_hourly_weekdays_november_profile mean_hourly_weekdays_december_profile];
mean_hourly_weekdays_overall_profile = mean(mean_hourly_weekdays_ALLMONTHS_profile,2);

% Plot the mean hourly profile of a working day for the entire year
figuree5e588poe7 = figure('visible','off');
plot(0:23, mean_hourly_weekdays_overall_profile,'LineWidth',3)
grid on
title('Mean Hourly Profile of a working day for the entire year');
xlabel('Hour of the day'); ylabel('kWh') ;
xlim([0 23]);
xticks(0:2:23)
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
savefig('overall_hourly_mean_working_profile.fig')
saveas(gcf,'overall_hourly_mean_working_profile.png')

%% Calculate the average profile of the weekends for each month and the average weekend load per month.

datestamp_weekend = datestamp(isweekend(datestamp) == 1); % identify the dates of Saturdays and Sundays
datestamp_weekendR = reshape(datestamp_weekend,[24 numel(datestamp_weekend)/24]);
weekend_days_number = day(datestamp_weekendR,'dayofyear');
weekend_days_number = transpose(weekend_days_number(1,:)); % identify the weekends day number in a year (e.g. 6,7,13,14 etc.)
weekend_hourly_indicesR = zeros(24, numel(datestamp_weekend)/24); % hourly index in a year (1 to 8760 possible)

for trt = 1:length(weekend_days_number)
    weekend_hourly_indicesR(:,trt) = (weekend_days_number(trt)-1)*24 +1 :24*weekend_days_number(trt);
end
clear trt

weekend_hourly_indices = reshape(weekend_hourly_indicesR,[numel(datestamp_weekend) 1]);
weekend_totalenergy = totalenergy(weekend_hourly_indices); % total energy in weekends.

% Setting 12 temporary matrices to identify logically if the month when the
% weekend takes place is January, February etc.

koukou1 =  month(datestamp(weekend_hourly_indices)) == 1;
koukou2 =  month(datestamp(weekend_hourly_indices)) == 2;
koukou3 =  month(datestamp(weekend_hourly_indices)) == 3;
koukou4 =  month(datestamp(weekend_hourly_indices)) == 4;
koukou5 =  month(datestamp(weekend_hourly_indices)) == 5;
koukou6 =  month(datestamp(weekend_hourly_indices)) == 6;
koukou7 =  month(datestamp(weekend_hourly_indices)) == 7;
koukou8 =  month(datestamp(weekend_hourly_indices)) == 8;
koukou9 =  month(datestamp(weekend_hourly_indices)) == 9;
koukou10 = month(datestamp(weekend_hourly_indices)) == 10;
koukou11 = month(datestamp(weekend_hourly_indices)) == 11;
koukou12 = month(datestamp(weekend_hourly_indices)) == 12;

totalenergy_weekend_january = totalenergy(weekend_hourly_indices(koukou1));
totalenergy_weekend_february = totalenergy(weekend_hourly_indices(koukou2));
totalenergy_weekend_march = totalenergy(weekend_hourly_indices(koukou3));
totalenergy_weekend_april = totalenergy(weekend_hourly_indices(koukou4));
totalenergy_weekend_may = totalenergy(weekend_hourly_indices(koukou5));
totalenergy_weekend_june = totalenergy(weekend_hourly_indices(koukou6));
totalenergy_weekend_july = totalenergy(weekend_hourly_indices(koukou7));
totalenergy_weekend_august = totalenergy(weekend_hourly_indices(koukou8));
totalenergy_weekend_september = totalenergy(weekend_hourly_indices(koukou9));
totalenergy_weekend_october = totalenergy(weekend_hourly_indices(koukou10));
totalenergy_weekend_november = totalenergy(weekend_hourly_indices(koukou11));
totalenergy_weekend_december = totalenergy(weekend_hourly_indices(koukou12));

totalenergy_weekend_januaryR = reshape(totalenergy_weekend_january,[24 numel(totalenergy_weekend_january)/24]);
totalenergy_weekend_februaryR = reshape(totalenergy_weekend_february,[24 numel(totalenergy_weekend_february)/24]);
totalenergy_weekend_marchR = reshape(totalenergy_weekend_march,[24 numel(totalenergy_weekend_march)/24]);
totalenergy_weekend_aprilR = reshape(totalenergy_weekend_april,[24 numel(totalenergy_weekend_april)/24]);
totalenergy_weekend_mayR = reshape(totalenergy_weekend_may,[24 numel(totalenergy_weekend_may)/24]);
totalenergy_weekend_juneR = reshape(totalenergy_weekend_june,[24 numel(totalenergy_weekend_june)/24]);
totalenergy_weekend_julyR = reshape(totalenergy_weekend_july,[24 numel(totalenergy_weekend_july)/24]);
totalenergy_weekend_augustR = reshape(totalenergy_weekend_august,[24 numel(totalenergy_weekend_august)/24]);
totalenergy_weekend_septemberR = reshape(totalenergy_weekend_september,[24 numel(totalenergy_weekend_september)/24]);
totalenergy_weekend_octoberR = reshape(totalenergy_weekend_october,[24 numel(totalenergy_weekend_october)/24]);
totalenergy_weekend_novemberR = reshape(totalenergy_weekend_november,[24 numel(totalenergy_weekend_november)/24]);
totalenergy_weekend_decemberR = reshape(totalenergy_weekend_december,[24 numel(totalenergy_weekend_december)/24]);

mean_hourly_weekend_january_profile = mean(totalenergy_weekend_januaryR,2);
mean_hourly_weekend_february_profile = mean(totalenergy_weekend_februaryR,2);
mean_hourly_weekend_march_profile = mean(totalenergy_weekend_marchR,2);
mean_hourly_weekend_april_profile = mean(totalenergy_weekend_aprilR,2);
mean_hourly_weekend_may_profile = mean(totalenergy_weekend_mayR,2);
mean_hourly_weekend_june_profile = mean(totalenergy_weekend_juneR,2);
mean_hourly_weekend_july_profile = mean(totalenergy_weekend_julyR,2);
mean_hourly_weekend_august_profile = mean(totalenergy_weekend_augustR,2);
mean_hourly_weekend_september_profile = mean(totalenergy_weekend_septemberR,2);
mean_hourly_weekend_october_profile = mean(totalenergy_weekend_octoberR,2);
mean_hourly_weekend_november_profile = mean(totalenergy_weekend_novemberR,2);
mean_hourly_weekend_december_profile = mean(totalenergy_weekend_decemberR,2);

mean_hourly_weekend_january = mean(mean_hourly_weekend_january_profile);
mean_hourly_weekend_february = mean(mean_hourly_weekend_february_profile);
mean_hourly_weekend_march = mean(mean_hourly_weekend_march_profile);
mean_hourly_weekend_april = mean(mean_hourly_weekend_april_profile);
mean_hourly_weekend_may = mean(mean_hourly_weekend_may_profile);
mean_hourly_weekend_junee = mean(mean_hourly_weekend_june_profile);
mean_hourly_weekend_july = mean(mean_hourly_weekend_july_profile);
mean_hourly_weekend_august = mean(mean_hourly_weekend_august_profile);
mean_hourly_weekend_september = mean(mean_hourly_weekend_september_profile);
mean_hourly_weekend_october = mean(mean_hourly_weekend_october_profile);
mean_hourly_weekend_november = mean(mean_hourly_weekend_november_profile);
mean_hourly_weekend_december = mean(mean_hourly_weekend_december_profile);

clear koukou1 koukou2 koukou3 koukou4 koukou5 koukou6 koukou7 koukou8 koukou9 koukou10 koukou11 koukou12

%% Plot the mean hourly energy profile of each month (excluding weekends)
% peak_mean_hourly_value is the highest hourly value per month.
% max_mean_hourly_value is the maximum load per month, noticed in the mean
% hourly electricity profile for weekdays.

peak_hourly_value_per_month = [max(january_month_load(:)); max(february_month_load(:)); max(march_month_load(:));...
    max(april_month_load(:)); max(may_month_load(:)); max(june_month_load(:)); max(july_month_load(:));...
    max(august_month_load(:)); max(september_month_load(:)); max(october_month_load(:));...
    max(november_month_load(:)); max(december_month_load(:))];

peak_hourly_value = max([january_month_load(:); february_month_load(:); march_month_load(:);...
    april_month_load(:); may_month_load(:); june_month_load(:); july_month_load(:);...
    august_month_load(:); september_month_load(:); october_month_load(:);...
    november_month_load(:); december_month_load(:)]);

max_mean_hourly_value = [max(mean_hourly_weekdays_january_profile) max(mean_hourly_weekdays_february_profile)...
     max(mean_hourly_weekdays_march_profile) max(mean_hourly_weekdays_april_profile) max(mean_hourly_weekdays_may_profile)...
     max(mean_hourly_weekdays_june_profile) max(mean_hourly_weekdays_july_profile) max(mean_hourly_weekdays_august_profile)...
     max(mean_hourly_weekdays_september_profile) max(mean_hourly_weekdays_october_profile) max(mean_hourly_weekdays_november_profile)...
     max(mean_hourly_weekdays_december_profile)]';
 
% Print the peak mean hourly value
fprintf('The peak mean hourly electricity load for the entire year is %.2f kW.\n',peak_hourly_value)
% Print the highest hourly value observed in totalenergy
datestamp_max = datestamp(totalenergy == max(totalenergy));
fprintf('The highest power value observed is %.2f kW. \n',max(totalenergy))

for ii = 1:numel(datestamp_max)
   fprintf('The highest power consumption takes place at %s.\n',datestamp_max(ii))
end
clear ii

figure8 = figure('visible','off');
x = 0:23;
subplot(3,4,1)
plot(x, mean_hourly_weekdays_january_profile,'LineWidth',3)
title('January')
grid on
ylabel('kWh')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23]);
ylim([0 peak_hourly_value]);
set(gca,'xticklabel',{[]}) 

subplot(3,4,2)
plot(x, mean_hourly_weekdays_february_profile,'LineWidth',3)
grid on
title('February')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23])
set(gca,'xticklabel',{[]}) 

subplot(3,4,3)
plot(x, mean_hourly_weekdays_march_profile,'LineWidth',3)
grid on
title('March')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23])
set(gca,'xticklabel',{[]}) 

subplot(3,4,4)
plot(x, mean_hourly_weekdays_april_profile,'LineWidth',3)
grid on
title('April')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23])
set(gca,'xticklabel',{[]}) 

subplot(3,4,5)
plot(x, mean_hourly_weekdays_may_profile,'LineWidth',3)
ylabel('kWh')
grid on
title('May')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23])
set(gca,'xticklabel',{[]}) 

subplot(3,4,6)
plot(x, mean_hourly_weekdays_june_profile,'LineWidth',3)
grid on
title('June')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23])
set(gca,'xticklabel',{[]}) 

subplot(3,4,7)
plot(x, mean_hourly_weekdays_july_profile,'LineWidth',3)
grid on
title('July')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23])
set(gca,'xticklabel',{[]}) 

subplot(3,4,8)
plot(x, mean_hourly_weekdays_august_profile,'LineWidth',3)
grid on
title('August')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([0 12 23])
set(gca,'xticklabel',{[]}) 

subplot(3,4,9)
plot(x, mean_hourly_weekdays_september_profile,'LineWidth',3)
grid on
title('September')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([12 23])
ylabel('kWh')

subplot(3,4,10)
plot(x, mean_hourly_weekdays_october_profile,'LineWidth',3)
grid on
title('October')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([12 23])

subplot(3,4,11)
plot(x, mean_hourly_weekdays_november_profile,'LineWidth',3)
grid on
title('November')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([12 23])

subplot(3,4,12)
plot(x, mean_hourly_weekdays_december_profile,'LineWidth',3)
grid on
title('December')
xlim([0 23]);
ylim([0 ceil(peak_hourly_value)])
xticks([12 23])

savefig('all_months_hourly_mean_weekdays_profile.fig')
saveas(gcf,'all_months_hourly_mean_weekdays_profile.png')
clear x

%% Create matrices for the entire year

peak_load_per_month = [max(january_month_load) max(february_month_load) max(march_month_load) max(april_month_load)...
    max(may_month_load) max(june_month_load) max(july_month_load) max(august_month_load) max(september_month_load)...
    max(october_month_load) max(november_month_load) max(december_month_load)]; 

%% Calculate the off-peak loads of the year
% This is the sum of the loads during the time slots (1-7) and (19-24).
% The preheating loads between 7-8am are included as a part of the off-peak
% loads.

offpeak_loads_scalar = (sum(sum(total_energyR(19:24,:))) + sum(sum(total_energyR(1:8,:))));
peak_loads_scalar = sum(totalenergy) - offpeak_loads_scalar;
%% Create identical matrices for the maximum life_system of the project.
% These assignments are necessary as the LIFE variables are used by the
% storage model. In this case, as we are interested in the building loads
% of one year, the first 8760 values are given to the LIFE variables.

datestampLIFE = datestamp;
datestamp2LIFE = datestamp2;
totalenergyLIFE = totalenergy;
daystampLIFE = daystamp;
index_dst1_LIFE = index_dst1;
index_dst2_LIFE = index_dst2;
datetime_dst1_LIFE = datetime_dst1;
datetime_dst2_LIFE = datetime_dst2;

%% Save variables to a file for later use
fname2 = sprintf('data_analysis_oneyear_%s.mat',ScenarioID);
save(fname2)
clear fname2
