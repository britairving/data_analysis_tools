# data_analysis_tools
MATLAB funcitons useful in many applications.

makefig.m                       | utility script used to create standard figure. Figure characteristics can easily be set such as position, fontsize, fontname, etc. 
makefig_subplots.m              | utility script used to create standard figure with multiple subplots. Axes positions and characterisitics very adaptable. 
standard_printfig_highrespng.m  | utility script used to save figures in a standard format. 
manual_timeSeries_QC.m          | Script to manually QC data in timeseries format. Useful for any timeseries data but created for moored oceanographic datasets.
gsw_rho_irving.m                | gsw_rho.m adapted to pass through pressure, temperature, and practical salinity (must be in specific data structure format). If uncertainties are passed through, calls gsw_errorprop_irving.m
gsw_errorprop_irving.m          | Calculate uncertainties for TEOS-10 potential temperature (pt), conservative temperature (CT), absolute salinity (SA), and density (rho) using the error propagation equations from Dai & Zhang 2018. [Dai, H., & Zhang, X. (2018). Uncertainties in climatologicalseawater density calculations.Journalof Geophysical Research: Oceans,123,2192–2212. https:././doi.org./10.1002./2017JC013427.]

