The first file (data_analysis_oneyear.m) imports DesignBuilder’s data generated from the building energy simulations, assigns them in variables and makes quick calculations. 
The second file (RTP_NEW_ONEYEAR.m) imports NordPool’s real-time pricing data to carry out the same objectives along with a quick statistical analysis.
The arbitrage model consists of three files, one for each operational strategy:
1a refers to Strategy E7 (exports allowed throughout the week)
1b refers to Strategy E5 (exports allowed only on working days)
2 refers to Strategy E0 (exports not allowed)

The arbitrage model is also available in its MATLAB function form where the appropriate input is needed (battery/rectifier/inverter capacities).
The function MATLAB files are indicated by 'MEGA' at the end of the filename
