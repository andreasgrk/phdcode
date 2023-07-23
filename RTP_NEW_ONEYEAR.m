%% Descriptor of the Functions used

% 1) allfitdist to choose and plot the best distributions for the RTP data
% 2) RTP_Distribution to plot the histogram of the empirical PDF of the RTP
% Retail data.

% Attention: RTP data only from Common Years must be imported. Leap years not supported 
% due to the incompatibility of DesignBuilder

%% INPUT FOR ELECTRICITY SPOT PRICES (ideally tab-delimited .txt file)
 % Both numbers that use commar or period as a decimal seperator are
 % supported. This is necessary as NordPool european data use the comma as
 % the decimal seperator.
 
 % Important: The timezone on the NordPool's N2EX hourly data must be
 % double checked as in 2016 and 2017, the price series start at 31/12
 % 23:00 of the previous year i.e. the first hour of the series must be replaced by the
 % last price of the data to be consistent with the UK time (CET - 1).
 
 %% Load basic data from data_analysis.m
%  load data_analysis_oneyear.mat % Use if necessary

%% Select DST Elimination
% For all cases, DST Elimination MUST take place to be consistent with
% DesignBuilder Values. Turn it off (False) only for validation/testing
% purposes.

DST_elimination = true;
 
 %% Select if input should refer to several years (life_system). 
 % If true, the RTP input of one year will be duplicated e.g. ten times for life_system = 10.
 % Setting of 'duplicate' takes place in the data_analysis script of the
 % building.
  
 %% INPUT from .txt file
 % IMPORTANT: The input RTP data must have 8760 values (for a single year).
RTP_filename = input('Enter the name for the Real Time Pricing file (.txt): ', 's');
tempdata = importdata(RTP_filename,'\t');
% if tempdata is cell, that means that data include comma as a seperator and need to be converted
if iscell(tempdata) == 1 
   numericData = [];
   for ppp = 1:size(tempdata,1)
     numericData(ppp,:) = str2double(char(strrep(tempdata(ppp,:),',','.')'));
   end
   RTP_wholesale = numericData;
else
    RTP_wholesale = tempdata; % in £/MWh
end
clear tempdata numericData ppp
  
%% No Duplication is necessary as this script is used for annual comparison purposes.
%% Fix the RTP_wholesale matrix in order to eliminate DST observation
% To achieve this, it is important for the prices to match their respective
% time period, regardless of the hour change and the following steps must
% be followed:

% a. Replace the empty time period (e.g. for 26/3/2017 1-2am) with the
%    average value of the previous and the following RTP price.

% b. Replace the value of the first time period that takes place twice
%    (e.g. 29/10/2017 1-2am) with their mean value and eliminate the second

if  DST_elimination == true
    
    RTP_wholesale = [RTP_wholesale(1:index_dst1-1); (RTP_wholesale(index_dst1-1) + RTP_wholesale(index_dst1))/2;...
    RTP_wholesale(index_dst1:index_dst2-2); (RTP_wholesale(index_dst2-1) + RTP_wholesale(index_dst2))/2;...
    RTP_wholesale(index_dst2+1:end)];
end
    
% It is assumed that the wholesale price constitutes 36.74% of the total bill (Ofgem) 
% Estimated from the Segmental Statements as of August 2018.
RTP_wholesale_kWh = RTP_wholesale/1000; % wholesale price in £ per kWh

RTP_wholesale_kWh_original = RTP_wholesale_kWh; % £ per kWh
RTP_wholesale_original = RTP_wholesale; % £ per MWh
% These two matrices will keep the original values in a matrix as some of the values may be modified 
% if they are below the minimum_wholesale_kWh_price_to_pay RTP_wholesale_original = RTP_wholesale;
% It should be pointed out that these original matrices does NOT have the
% DST elimination corrections
 
wholesale_percentage = 0.3663; % Data from Ofgem for non-domestic electricity (2017). Percentage remains the same for 2018 as the difference is less than 0.9%
minimum_wholesale_kWh_price_to_pay = 0.004; % Per kWh. If the original wholesale price is less than £0.004/kWh
% then the equivalent retail price of (£0.004 wholesale) will be paid in order to cover the
% distribution & network costs.
% For example, if the actual price is £0.001/kWh, then this will be revised to
% £0.004 and the equivalent retail price will be 0.004/wholesale_percentage.
% This means that there will be a minimum positive retail price.

RTP_retail = zeros(length(RTP_wholesale_kWh),1);
for ty = 1:length(RTP_wholesale_kWh)
    if RTP_wholesale_kWh(ty) > minimum_wholesale_kWh_price_to_pay 
       RTP_retail(ty) = RTP_wholesale_kWh(ty)/wholesale_percentage;
    elseif RTP_wholesale_kWh(ty) <-minimum_wholesale_kWh_price_to_pay
           RTP_retail(ty) = RTP_wholesale_kWh(ty)*wholesale_percentage;
    elseif (RTP_wholesale_kWh(ty) <= minimum_wholesale_kWh_price_to_pay  && RTP_wholesale_kWh(ty) >=0) ||... % zero is also included
           (RTP_wholesale_kWh(ty) >= -minimum_wholesale_kWh_price_to_pay  && RTP_wholesale_kWh(ty) <0)
        RTP_retail(ty) = abs(minimum_wholesale_kWh_price_to_pay /wholesale_percentage);
        RTP_wholesale_kWh(ty) = abs(minimum_wholesale_kWh_price_to_pay);
    end
end
clear ty
RTP_wholesale = RTP_wholesale_kWh * 1000; % To update the £/MWh wholesale values 
% in case the respective £/kWh wholesale values were replaced from the previous if-statements
RTP_retailR = reshape(RTP_retail, [24 length(RTP_retail)/24]);
RTP_retail_MWh = RTP_retail * 1000; % £/MWh retail price.
%% Plot RTP_retail against datestamp for the first year

figure9 = figure('visible','off');
plot(datestamp, RTP_retail(1:8760),'LineWidth',1)
grid on
title('Day-ahead Retail Real Time Pricing')
xlabel('Date'); ylabel('£/kWh') ;
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
savefig('RTP_retail.fig')
saveas(gcf,'RTP_retail.png')

%% Also plot wholesale price against datestamp for the first year
figure156458 = figure('visible','off');
plot(datestamp, RTP_wholesale(1:8760),'LineWidth',1)
grid on
title('Day-ahead Wholesale Real Time Pricing')
xlabel('Date'); ylabel('£/MWh') ;
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
savefig('RTP_wholesale.fig')
saveas(gcf,'RTP_wholesale.png')

%% Clear temporary variables
clearvars delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

%% CREATE THE HISTOGRAM OF THE FREQUENCY DENSITY FOR THE FIRST YEAR
% Prepare figure
if DST_elimination == true
figure859634 = figure('visible','off');
LegHandles = []; LegText = {};

% --- Plot data originally in dataset "RTP_retail data"
[CdfF,CdfX] = ecdf(RTP_retail,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(RTP_retail,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0 0.45098 0.741176],...
    'LineStyle','-', 'LineWidth',1);
xlabel('RTP Retail Electricity Tariff (£/kWh)');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'RTP_retail data';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XGrid = linspace(XLim(1),XLim(2),100);

% Adjust figure
box on;
grid on;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
set(hLegend,'Interpreter','none');

% Save
savefig('RTP_retail_empirical_PDF.fig')
saveas(gcf,'RTP_retail_empirical_PDF.png')
%Clear 
clearvars XLim XGrid LegText LegHandles hLine hLegend CdfF CdfX BinCenter BinEdge BinHeight BinInfo

figure2324543 = figure('visible','off');
RTP_distributions = allfitdist(RTP_retail(1:8760),'PDF');
xlabel('RTP Retail Price (£/kWh)')
xlim auto
ylim auto
savefig('allfitdist_RTP.fig')
saveas(gcf,'allfitdist_RTP.png')

figure059683 = figure('visible','off');
histogram(RTP_retail)
xlabel('RTP Retail Electricity Price (£/kWh)')
ylabel('Absolute Frequency')
xlim auto
ylim auto
savefig('RTP_retail_histogram_frequency.fig')
saveas(gcf,'RTP_retail_histogram_frequency.png')

end

%% Create a matrix with the minimum and the maximum RTP price for each day. Plot for the first year.

n = 24; %set the number of hourly values in a day

RTP_retail_daily_min = arrayfun(@(i) min(RTP_retail(i:i+n-1)),1:n:length(RTP_retail)-n+1)';
RTP_retail_daily_max = arrayfun(@(i) max(RTP_retail(i:i+n-1)),1:n:length(RTP_retail)-n+1)';
RTP_retail_daily_average = arrayfun(@(i) mean(RTP_retail(i:i+n-1)),1:n:length(RTP_retail)-n+1)';
RTP_retail_difference = RTP_retail_daily_max - RTP_retail_daily_min;
RTP_retail_daily_minR = RTP_retail_daily_min';
RTP_retail_daily_maxR = RTP_retail_daily_max';
RTP_retail_differenceR = RTP_retail_difference';
RTP_retail_daily_averageR = RTP_retail_daily_average';
clear n

if DST_elimination == true % This is to avoid endless loops in case the RTP prices provided are for testing/validation purposes
figure7989755 = figure('visible','off');
RTP_retail_difference_allfitdist = allfitdist(RTP_retail_difference,'PDF');
% The fit includes all the RTP data, NOT only the first year.
xlabel('Maximum - Minimum Retail Price (£/kWh)')
xlim auto
ylim auto
savefig('allfitdist_RTP_difference.fig')
saveas(gcf,'allfitdist_RTP_difference.png')  
end

% Create extra variables just in case
% Create RTP min, max and difference matrices for all the years involved.
RTP_retail_daily_maxR_FULL = zeros(24,length(RTP_retail)/24);
RTP_retail_daily_minR_FULL = zeros(24,length(RTP_retail)/24);
RTP_retail_differenceR_FULL = zeros(24,length(RTP_retail)/24);

for ff = 1:length(RTP_retail)/24
   RTP_retail_daily_maxR_FULL(:,ff) = RTP_retail_daily_maxR(ff);
   RTP_retail_daily_minR_FULL(:,ff) = RTP_retail_daily_minR(ff);
   RTP_retail_differenceR_FULL(:,ff) = RTP_retail_differenceR(ff);
end
clear ff

RTP_retail_daily_max_FULL = reshape(RTP_retail_daily_maxR_FULL,[numel(RTP_retail_daily_maxR_FULL),1]);
RTP_retail_daily_min_FULL = reshape(RTP_retail_daily_minR_FULL,[numel(RTP_retail_daily_maxR_FULL),1]);
RTP_retail_difference_FULL = reshape(RTP_retail_differenceR_FULL,[numel(RTP_retail_daily_maxR_FULL),1]);

%% Plot min, max, mean and the different(min-max) of each day for the first year
figure1555i562 = figure('visible','off');
plot(daystamp(1:365), RTP_retail_daily_max(1:365), daystamp(1:365), RTP_retail_daily_min(1:365), daystamp(1:365), RTP_retail_daily_average(1:365),'LineWidth',2)
grid on
title('Day-ahead retail max, min and average price per day')
xlabel('Date'); ylabel('Retail electricity price (£/kWh)') ;
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
xlim auto; ylim auto; legend('Maximum Price', 'Minimum Price','Mean Price','Location','northeast')
savefig('RTP_max_min.fig')
saveas(gcf,'RTP_max_min.png')

%% Plot the (maximum - minimum) price of each day, for the first year.
figure15465du = figure('visible','off');
plot(daystamp(1:365), RTP_retail_difference(1:365),'LineWidth',3)
grid on
title('Day-ahead (Maximum - Minimum) RTP Price Difference')
xlabel('Date'); ylabel('£/kWh') ;
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))
xtickangle(45)
xlim auto; ylim auto;
savefig('RTP_difference.fig')
saveas(gcf,'RTP_difference.png')

%% Save all RTP data results
% Save variables into a file for further use
fname55 = sprintf('RTP_analysis_%i_%s.mat',Year,ScenarioID);
save(fname55);
fname56 = sprintf('RTP_analysis_%i.mat',Year);
save(fname56);
clear fname55
clear fname56
fname58 = sprintf('RTP_analysis.mat');
save(fname58)
clear fname58