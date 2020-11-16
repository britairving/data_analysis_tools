function d = gsw_rho_irving(datastruct)
%% function gsw_rho_irving
%% gsw_rho                                in-situ density (75-term equation)
%==========================================================================
% 
% USAGE:  
%  rho = gsw_rho_irving(datastrut)
%
% DESCRIPTION:  
%
% INPUT:
%   datastruct = data structure containing fields... 
%       p    =  pressure                 [ dbar ]
%       t    =  temperature              [ ITS-90 deg C ]
%       Sp   =  practical salinty        [ ]
%       lat  =  latitude                 [decimal degrees]
%       lon  =  longitude                [decimal degrees]
%       u_p  =  pressure standard uncertainty           [ dbar ]   ** OPTIONAL **
%       u_t  =  temperature standard uncertainty        [ ITS-90 deg C ] ** OPTIONAL **
%       u_Sp =  practical salinity standard uncertainty [ ] ** OPTIONAL **
%
% OUTPUT:
%  rho  =  in-situ density                                         [ kg/m3 ]
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th November, 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
%
%
% Author:   Brita Irving <bkirving@alaska.edu>
% Created:  July 2019
%% 
if ~isstruct(datastruct)
  error('Expected data structure as input argument')
end
%% Use abbreviated name for structure
d = datastruct;

%% Calculate Absolute Salinity [g/kg] from Practical Salinity (PSS-78) [unitless]
d.SA  = gsw_SA_from_SP(d.Sp, d.p, d.lon, d.lat); % [SA, in_ocean] = gsw_SA_from_SP(Sp,p,long,lat) 

%% Calculate potential temperature (ITS-90) [deg C] from in-situ temperature (ITS-90) [deg C]
if isfield(d,'p_ref')
  p_ref = d.p_ref;
else
  p_ref = 0; % [dbar]
end
d.pt  = gsw_pt_from_t(d.SA, d.t, d.p, p_ref); % pt = gsw_pt_from_t(SA,t,p,p_ref) [ deg C ]

%% Calculate Conservative temperature [deg C] from potential temperature (ITS-90) [deg C]
d.CT  = gsw_CT_from_pt(d.SA,d.pt);            % CT = gsw_CT_from_pt(SA,pt) [ deg C ]

%% Calculate potential density of seawater  (not potential density anomaly) [ kg/m^3 ]
d.pot_rho = gsw_pot_rho_t_exact(d.SA,d.t,d.p,p_ref);

%% Calculate in-situ density [kg/m] from Absolute Salinity [g/kg], Conservative Temperature [deg C], and pressure [dbar]
d.rho = gsw_rho(d.SA, d.CT, d.p);             % rho = gsw_rho(SA,CT,p)  in-situ density   [ kg/m ]

%% Error propagation
if isfield(d,'u_Sp') && isfield(d, 'u_t')
   u = gsw_errorprop_irving(d);
   d = u; % update data structure with uncertainties
end

end %% MAIN FUNCTION gsw_rho_irving