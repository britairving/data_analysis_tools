function u = gsw_errorprop_irving(dstruct)
%% function gsw_errorprop_irving
% DESCRIPTION: 
%   Calculate uncertainties for TEOS-10 potential temperature (pt),
%   conservative temperature (CT), absolute salinity (SA), and density
%   (rho) using the error propagation equations from Dai & Zhang (2018).
%
% INPUT:
%   dstruct = data structure containing fields... 
%       p    =  pressure                 [ dbar ]
%       t    =  temperature              [ ITS-90 deg C ]
%       Sp   =  practical salinty        [ ]
%       SA   =  absolute salinity        [ g/kg ]
%       pt   =  potential temperature    [ ITS-90 deg C ]
%       CT   =  conservative temperature [ ITS-90 deg C ]
%       rho  =  in-situ density          [ kg/m^3 ]
%       u_p  =  pressure standard uncertainty     [ dbar ]   ** OPTIONAL **
%       u_t  =  temperature standard uncertainty  [ ITS-90 deg C ]
%       u_Sp =  practical salinity standard uncertainty [ ]
%
% OUTPUT:
%   u = data structure containing fields... 
%       p    =  pressure                 [ dbar ]
%       t    =  temperature              [ ITS-90 deg C ]
%       Sp   =  practical salinty        [ ]
%       SA   =  absolute salinity        [ g/kg ]
%       pt   =  potential temperature    [ ITS-90 deg C ]
%       CT   =  conservative temperature [ ITS-90 deg C ]
%       rho  =  in-situ density          [ kg/m^3 ]
%       u_p  =  pressure standard uncertainty     [ dbar ] 
%       u_t  =  temperature standard uncertainty  [ ITS-90 deg C ]
%       u_Sp =  practical salinity standard uncertainty [ ]
%       uc_pt  = combined standard uncertainty of pt    [ ITS-90 deg C ]
%       uc_CT  = combined standard uncertainty of CT    [ ITS-90 deg C ]
%       uc_SA  = combined standard uncertainty of SA    [ g/kg ]
%       uc_rho = combined standard uncertainty of rho   [ kg/m^3 ]
%
% AUTHOR:   
%   Brita Irving <bkirving@alaska.edu>
%
% REFERENCES:
%   Dai, H., & Zhang, X. (2018).Uncertainties in
%       climatologicalseawater density calculations.Journalof Geophysical
%       Research: Oceans,123,2192–2212.
%       https:././doi.org./10.1002./2017JC013427
%   IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%       seawater - 2010: Calculation and use of thermodynamic properties.  
%       Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%       UNESCO (English), 196 pp.  Available from http:././www.TEOS-10.org
%       See appendix A.20 and appendix K of this TEOS-10 Manual. 
%   Jackett, D. R., McDougall, T. J., Feistel, R., Wright, D. G., &
%       Griffies, S. M. (2006). Algorithms for density, potential
%       temperature, conservative temperature and the freezing temperature
%       of seawater. Journal Atmospheric Oceanic Technology, 23(12),
%       1709–1728. https:doi.org/10.1175/JTECH1946.1
%   Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%       polynomial expressions for the density and specifc volume of
%       seawater using the TEOS-10 standard. Ocean Modelling, 90, pp.
%       29-43.

%% Make sure data input in correct format & Define structure that contains all data and coefficients
if ~isstruct(dstruct)
  error('expected input variable to be a structure')
end
d = struct();
%% Define constants and initialize variables [Table 1. from Dai & Zhang (2018)]
d.Sp      = [];          % practical salinity
d.u_Sp    = [];          % practical salinity uncertainty
d.S35A    = 35.16504;    % [g/kg] Absolute Salinity of Reference Composition Seawater
d.ups     = d.S35A./35;  % Conversion factor of Practical Salinity to Reference Salinity
d.SR      = d.ups.*d.Sp; % [g/kg] Reference Salinity
d.SA      = [];          % Absolute Salinity (g./kg)
d.delSA   = [];          % SA - SR (g/kg) Absolute Salinity anomaly
d.t       = [];          % ITS-90 (C) In situ temperature
d.u_t     = [];          % ITS-90 (C) In situ temperature uncertainty
d.T0      = 273.15;      % K Degree Celsius zero point
d.pt      = [];          % Potential temperature (C)    [symbol = theta]
d.CT      = [];          % Conservative temperature (C) [symbol = Omega]
d.p       = [];          % Pressure (dbar)
d.u_p     = [];          % Pressure uncertainty
d.pr      = [];          % Reference pressure (dbar)
d.v       = [];          % specific volume (m^3/kg)
d.rho     = [];          % In situ density (kg/m^3) [symbol = rho]
d.C       = [];          % Conductivity (mS/cm)
d.cp0     = 3991.86795711963; % [J kg^-1 K^-1] Specific heat constant
d.u_S35A  = 0.00136;     % [g/kg] Millero et al. (2008)
d.u_delSA = 0.003;       % [g/kg] standard uncertainty of Absolute Salinity anomaly non-Baltic seas Pawlowicz et al. (2011)
d.SAu     = 40.*d.ups;   % [g/kg]
d.u_SA_forumla = 0.0024; % [g/kg]
d.CTu     = 40;          % [deg C]
d.pu      = 10.^4;       % [dbar]


%% correlation coefficients 
% approximate values taken around North Pacific/Gulf of Alaska
% correlation coeff's defined as 0 were not published in Dai & Zhang 2018
d.r_tp       = -0.8;   % r(t,p)   Figure 3  [Dai & Zhang 2018]
d.r_pt0t     = 0;      % r(d.pt0,d.t)
d.r_SAp      = 0.8;    % r(SA,p)  Figure 5  [Dai & Zhang 2018]
d.r_pt0f_pt0 = 0;      % r(pt0,f_pt0) 
d.r_SApt     = -0.8;   % r(SA,pt) Figure 9  [Dai & Zhang 2018]
d.r_pt0SA    = 0;      % r(pt0,SA) 
d.r_CTp      = -0.8;   % r(CT,p)  Figure 14 [Dai & Zhang 2018]
d.r_SACT     = -0.8;   % r(SA,CT) Figure 13 [Dai & Zhang 2018]

%% Read in input values
% Replaces the initialized variables with input arguments
fnames = fieldnames(dstruct);
for nf = 1:numel(fnames)
  if isfield(d, fnames{nf})
    d.(fnames{nf}) = dstruct.(fnames{nf});
  %elseif isfield(mtch,fnames{nf})
  %  d.(mtch.(fnames{nf})) =  dstruct.(fnames{nf});
  end
end

%% Calculate specific volume
% s, pi, and tau defined on pg#2203 Dai & Zhang 2018
d.v   = gsw_specvol(d.SA,d.CT,d.p);
d.s   = sqrt((d.SA + 24)./d.SAu);    
d.pi  = d.p./d.pu;
d.tau = d.CT./d.CTu;

%% Jackett et al. (2006) Table A1: Coefficients of the seven-term approximating polynomial pt0 for potential temperature 
% coefficients ai where i=1...7 but here named ji so do not confuse with
% a1...a29 coefficients in Dai & Zhang 2018
j1 = 8.65483913395442e-6;
j2 = -1.41636299744881e-6;
j3 = -7.38286467135737e-9;
j4 = -8.38241357039698e-6;
j5 = 2.83933368585534e-8;
j6 = 1.77803965218656e-8;
j7 = 1.71155619208233e-10;

%% Table A1: Coefficients and the Corresponding Values in Equations A.3.1-A.3.3
a1 = 61.01362420681071;
a2 = 168776.46138048015;
a3 = -22735.2785605119625;
a4 = 2574.2164453821433;
a5 = -21536.6644434977543;
a6 = 545.7340497931629;
a7 = 250.91091728474331;
a8 = 218.30489878927802;
a9 = 268.5520265845071;
a10 = 212019.028203559312;
a11 = 3734.858026725145;
a12 = 22046.7671145057618;
a13 = 465.28655623826234;
a14 = 20.6370820302376359;
a15 = 210.650848542359153;
a16 = 937.2099110620707;
a17 = 588.1802812170108;
a18 = 248.39476522971285;
a19 = 23.871557904936333;
a20 = 22.6268019854268356;
a21 = 21687.914374187449;
a22 = 246.9598888781377;
a23 = 123.59576582457964;
a24 = 248.5891069025409;
a25 = 936.3206544460336;
a26 = 2942.7827304544439;
a27 = 369.4389437509002;
a28 = 233.83664947895248;
a29 = 29.987880382780322;

%% Table K.1. Coefficients of the 75-term polynomial of Roquet et al. (2015)
% http://www.teos-10.org/pubs/TEOS-10_Manual.pdf Table K.1. Values copied
% from
% https://github.com/fabien-roquet/polyTEOS/blob/master/polyTEOS10_75t.m
% v006  was not in the polyTEOS10_75t.m script but was in the TEOS-10
% manual Table K.1.
% specific volume anomaly
v000 = 1.0769995862e-03; v100 = -3.1038981976e-04; v200 = 6.6928067038e-04; 
v300 = -8.5047933937e-04; v400 = 5.8086069943e-04; v500 = -2.1092370507e-04; 
v600 = 3.1932457305e-05; v010 = -1.5649734675e-05; v110 = 3.5009599764e-05; 
v210 = -4.3592678561e-05; v310 = 3.4532461828e-05; v410 = -1.1959409788e-05; 
v510 = 1.3864594581e-06; v020 = 2.7762106484e-05; v120 = -3.7435842344e-05; 
v220 = 3.5907822760e-05; v320 = -1.8698584187e-05; v420 = 3.8595339244e-06; 
v030 = -1.6521159259e-05; v130 = 2.4141479483e-05; v230 = -1.4353633048e-05; 
v330 = 2.2863324556e-06; v040 = 6.9111322702e-06; v140 = -8.7595873154e-06; 
v240 = 4.3703680598e-06; v050 = -8.0539615540e-07; v150 = -3.3052758900e-07; 
v060 = 2.0543094268e-07; v001 = -1.6784136540e-05; v101 = 2.4262468747e-05; 
v201 = -3.4792460974e-05; v301 = 3.7470777305e-05; v401 = -1.7322218612e-05; 
v501 = 3.0927427253e-06; v011 = 1.8505765429e-05; v111 = -9.5677088156e-06; 
v211 = 1.1100834765e-05; v311 = -9.8447117844e-06; v411 = 2.5909225260e-06; 
v021 = -1.1716606853e-05; v121 = -2.3678308361e-07; v221 = 2.9283346295e-06; 
v321 = -4.8826139200e-07; v031 = 7.9279656173e-06; v131 = -3.4558773655e-06; 
v231 = 3.1655306078e-07; v041 = -3.4102187482e-06; v141 = 1.2956717783e-06; 
v051 = 5.0736766814e-07; v002 = 3.0623833435e-06; v102 = -5.8484432984e-07; 
v202 = -4.8122251597e-06; v302 = 4.9263106998e-06; v402 = -1.7811974727e-06; 
v012 = -1.1736386731e-06; v112 = -5.5699154557e-06; v212 = 5.4620748834e-06; 
v312 = -1.3544185627e-06; v022 = 2.1305028740e-06; v122 = 3.9137387080e-07; 
v222 = -6.5731104067e-07; v032 = -4.6132540037e-07; v132 = 7.7618888092e-09; 
v042 = -6.3352916514e-08; v003 = -3.8088938393e-07; v103 = 3.6310188515e-07; 
v203 = 1.6746303780e-08; v013 = -3.6527006553e-07; v113 = -2.7295696237e-07; 
v023 = 2.8695905159e-07; v004 = 8.8302421514e-08; v104 = -1.1147125423e-07; 
v014 = 3.1454099902e-07; v005 = 4.2369007180e-09; v006 = 1.9613503930e-9;

%% Potential Temperature
% coefficients from Jackett et al. (2006) Table A1
d.pt0 = d.t + d.p.*(j1 + j2.*(d.SA.*35./d.S35A)+ j3.*d.p +j4.*d.t + j5.*(d.SA.*35./d.S35A).*d.t + j6.*d.t.^2 + j7.*d.t.*d.p);

%% Uncertainties in Absolute Salinity Calculations
% combined standard uncertainty in SA [Dai & Zhang (2018) Eqn 5]
% uc_SA = sqrt(u_SR.^2 + u_delSA.^2 + u_SA_forumla.^2);
% uc_SA = sqrt( (S35A./35.*u_Sp).^2 + (Sp./35.*u_S35A).^2 + u_delSA.^2 + u_SA_forumla.^2);
uc_SA = sqrt( (d.S35A./35.*d.u_Sp).^2 + (d.Sp./35.*d.u_S35A).^2 + d.u_delSA.^2 + d.u_SA_forumla.^2);

%% Pressure Uncertainty
if isempty(d.u_p)
    % pressure standard uncertainty  [Dai & Zhang (2018) Eqn 6]
    % u_p = sqrt(u_pref.^2 + u_pini.^2 + u_psta.^2);
    d.u_pref = 0.06 + 0.000065.*d.p; % [dbar] (Le Menn 2011)
    d.u_pini = 0.08;
    d.u_psta = 1.10;
    d.u_p = sqrt(d.u_pref.^2 + d.u_pini.^2 + d.u_psta.^2);
end

%% Uncertainties in Conservative Temperature Calculations CT = F(SA, pt);
[dCTdSA,dCTdpt,dCTdS35A] = CT_partialdiffyQ;

uc_pt0 = sqrt((dpt0dSA.*uc_SA).^2 + (dpt0dt.*d.u_t).^2 + (dpt0dp.*d.u_p).^2 +  ...
  (dpt0dS35A.*d.u_S35A).^2 + 2.*(dpt0dSA.*dpt0dt.*uc_SA.*d.u_t.*d.r_pt0t + ...
  dpt0dSA.*dpt0dp.*uc_SA.*d.u_p.*d.r_SAp+ dpt0dt.*dpt0dp.*d.u_t.*d.u_p.*d.r_tp)); % [Dai & Zhang (2018) Eqn 9]
            
uc_f_pt0 = sqrt( (df_pt0dpt0.*uc_pt0).^2 + (df_pt0dSA.*uc_SA).^2 + ...
              (df_pt0dS35A.*d.u_S35A).^2 + 2.*(df_pt0dpt0.*df_pt0dSA.*uc_pt0.*uc_SA.*d.r_pt0SA) ); % uc_SA instead of u_SA adapted from [Dai & Zhang (2018) Eqn 10]

% Uncertainty in potential temperature pt
uc_pt = sqrt( uc_pt0.^2 + uc_f_pt0.^2 + 2.*uc_pt0.*uc_f_pt0.*d.r_pt0f_pt0); % [Dai & Zhang (2018) Eqn 8]

% Uncertainty in conservative temperature CT
uc_CT = sqrt( (dCTdSA.*uc_SA).^2 + (dCTdpt.*uc_pt).^2 + (dCTdS35A.*d.u_S35A).^2 + ...
              2.*dCTdSA.*dCTdpt.*uc_SA.*uc_pt.*d.r_SApt ); % [Dai & Zhang (2018) Eqn 11]

%% Uncertainty in Density Calculations
% calculate combined standard uncertainty in density using  Dai & Zhang (2018) Eqn 13
[drhodSA, drhodCT, drhodp, drhodS35A] = rho_partialdiffyQ;
uc_rho = sqrt( (drhodSA.*uc_SA).^2   + (drhodCT.*uc_CT).^2 + (drhodp.*d.u_p).^2 + (drhodS35A.*d.u_S35A).^2 + ...
  2.*(drhodSA.*drhodCT.*uc_SA.*uc_CT.*d.r_SACT + drhodSA.*drhodp.*uc_SA.*d.u_p.*d.r_SAp + drhodCT.*drhodp.*uc_CT.*d.u_p.*d.r_CTp) );  % [Dai & Zhang (2018) Eqn 13]

%% save values to uncertainty structure u
u = dstruct();
u.u_p       = d.u_p;
u.uc_pt     = uc_pt;
u.uc_CT     = uc_CT;
u.uc_SA     = uc_SA;
u.uc_rho    = uc_rho;

%% A1. Partial Derivatives Taken by pt0 With Respect to SA, S35A, t, and p
% coefficients ji to j7 from Jackett et al. (2006) Table A1
  function val = dpt0dSA
    val = 35./d.S35A.*d.p.*(j2+j5.*d.t); % eqn A1-1
  end

  function val = dpt0dS35A
    val = (-35.*d.p.*d.SA./d.S35A.^2).*(j2 + j5.*d.t); % eqn A1-2
  end

  function val = dpt0dt
    val = 1 + d.p.*(j4 + 35./d.S35A.*j5.*d.SA + 2.*j6.*d.t + j7.*d.p); % eqn A1-3
  end

  function val = dpt0dp
    val = j1 + 35./d.S35A.*j2.*d.SA + 2.*j3.*d.p + j4.*d.t + 35./d.S35A.*j5.*d.SA.*d.t + j6.*d.t.^2 + 2.*j7.*d.p.*d.t;
  end

%% A2. Partial Derivatives Taken by f'pt0 With Respect to pt0, SA, and S35A
  function val = df_pt0dpt0
    val = -(d.cp0./(d.T0 + d.pt0).^2).*(1-0.05.*(1 - d.SA./d.S35A)).^(-1); % eqn A2-1
  end

  function val = df_pt0dSA
    val = -(0.05.*d.cp0./(d.S35A.*(d.T0 + d.pt0))).*(1-0.05.*(1 - d.SA./d.S35A)).^(-2); % eqn A2-2
  end

  function val = df_pt0dS35A
    val = (0.05.*d.SA.*d.cp0)./(d.S35A.^2.*(d.T0 + d.pt0)).*(1-0.05.*(1 - d.SA./d.S35A)).^(-2); % eqn A2-3
  end

%% A3. Partial Derivatives Taken by CT With Respect to SA, pt, and S35A
  function  [dCTdSA,dCTdpt,dCTdS35A] = CT_partialdiffyQ
    %% partial derivative of CT with respect to SA
    % dCTdSA = dCTdh0.*dh0dSA;  % eqn A3-1
    % Coefficients from Table A1 [Dai & Zhang 2018]
    x   = sqrt(d.SA./d.SAu); % x.^2 = d.SA./d.SAu;
    y   = d.pt./40;
    dCTdh0 = 1./d.cp0;
    dh0dx2 = a9 + y.*(a10 + y.*(a11 + y.*(a12 + y.*(a13 + y.*(a14 + a15.*y))))) + ...
      x.*(a16 + y.*(a17 + y.*(a18 + y.*(a19 + a20.*y))) + ...
      x.*(a21 + x.*(a22 + x.*(a23 + a24.*x)) + ...
      y.*(a25 + y.*(a26 + y.*(a27 + y.*(a28 +a29.*y))))));
    dx2dSA = 1./d.SAu;
    dh0dx  = x.^2.*(a16 + y.*(a17 + y.*(a18 + y.*(a19 + a20.*y)))) + ...
      (2.*a21.*x + 3.*a22.*x.^2 + 4.*a23.*x.^3 + 5.*a24.*x.^4) + ...
      2.*x.*y.*(a25 + y.*(a26 + y.*(a27 + y.*(a28 + a29.*y))));
    dxdSA = 1./(2.*x.*d.SAu);
    
    dh0dSA = dh0dx2.*dx2dSA + dh0dx.*dxdSA;
    
    dCTdSA = dCTdh0.*dh0dSA;  % eqn A3-1
    
    %% partial derivative of CT with respect to pt
    % dCTdpt = dCTdh0.* dh0dpt; % eqn A3-2
    dh0dy  = (a2 + 2.*a3.*y + 3.*a4.*y.^2 + 4.*a5.*y.^3 + 5.*a6.*y.^4 + 6.*a7.*y.^5 + 7.*a8.*y.^6) + ...
      x.^2.*((a10 + 2.*a11.*y + 3.*a12.*y.^2 + 4.*a13.*y.^3 + 5.*a14.*y.^4 + 6.*a15.*y.^5) + ...
      x.*((a17 + 2.*a18.*y + 3.*a19.*y.^2 + 4.*a20.*y.^3) + ...
      x.*(a25 + 2.*a26.*y + 3.*a27.*y.^2 + 4.*a28.*y.^3 + 5.*a29.*y.^4)));
    dydpt  = 1./40;
    dh0dpt = dh0dy.*dydpt;
    dCTdpt = dCTdh0.* dh0dpt; % eqn A3-2
    
    %% partial derivative of CT with respect to S35A
    % dCTdS35A = dCTdh0.*dh0dS35A; % eqn A3-3
    dx2dS35A = -d.SA./(d.S35A.*d.SAu);
    dxdS35A  = 1./(2.*x).*dx2dS35A;
    dh0dS35A = dh0dx2.*dx2dS35A + dh0dx.*dxdS35A;
    dCTdS35A = dCTdh0.*dh0dS35A; % eqn A3-3
  end % FUNCTION [dCTdSA,dCTdpt,dCTdS35A] = CT_partialdiffyQ

%% A4. Partial Derivatives Taken by rho With Respect to SA, CT, p, and S35A
  function [drhodSA, drhodCT, drhodp, drhodS35A] = rho_partialdiffyQ
    %% partial derivative of rho with respect to SA
    % drhodSA = drhodv.*dvdSA = drhodv.*dvds.*dsdSA; %eqn A4-1
    drhodv  = -1./d.v.^2;
    dsdSA   = 1./(2.*d.SAu).*sqrt(d.SAu./(d.SAu + 24));
    dvds    = v100 + 2.*v200.*d.s + 3.*v300.*d.s.^2 + 4.*v400.*d.s.^3 + 5.*v500.*d.s.^4 + 6.*v600.*d.s.^5 + ...
              d.tau.*   (v110 + 2.*v210.*d.s + 3.*v310.*d.s.^2 + 4.*v410.*d.s.^3 + 5.*v510.*d.s.^4) + ...
              d.tau.^2.* (v120 + 2.*v220.*d.s + 3.*v320.*d.s.^2 + 4.*v420.*d.s.^3) + ...
              d.tau.^3.* (v130 + 2.*v230.*d.s + 3.*v330.*d.s.^2) + ...
              d.tau.^4.* (v140 + 2.*v240.*d.s) + ...
              d.tau.^5.* (v150) + ...
              d.pi.*      (v101 + 2.*v201.*d.s + 3.*v301.*d.s.^2 + 4.*v401.*d.s.^3 + 5.*v501.*d.s.^4) + ...
              d.tau.*d.pi.*  (v111 + 2.*v211.*d.s + 3.*v311.*d.s.^2 + 4.*v411.*d.s.^3) + ...
              d.tau.^2.*d.pi.*(v121 + 2.*v221.*d.s + 3.*v321.*d.s.^2) + ...
              d.tau.^3.*d.pi.*(v131 + 2.*v231.*d.s) + ...
              d.tau.^4.*d.pi.*(v141) + ...
              d.pi.^2.*      (v102 +2.*v202.*d.s + 3.*v302.*d.s.^2 + 4.*v402.*d.s.^3) + ...
              d.tau.*d.pi.^2.*  (v112 + 2.*v212.*d.s + 3.*v312.*d.s.^2) + ...
              d.tau.^2.*d.pi.^2.*(v122 + 2.*v222.*d.s) + ...
              d.tau.^3.*d.pi.^2.*(v132) + ...
              d.pi.^3.*      (v103 + 2.*v203.*d.s) + ...
              d.tau.*d.pi.^3.*  (v113) + ...
              d.pi.^4.*      (v104);
    
    dvdSA   = dvds.*dsdSA;
    drhodSA = drhodv.*dvdSA;
    
    %% partial derivative of rho with respect to CT
    % drhodCT = drhodv.*dvdCT = drhodv.*dvdtau.*dtaudCT; % eqn A4-2
    dtaudCT = 1./d.CTu;
    dvdtau = v010 + v110.*d.s + v210.*d.s.^2 + v310.*d.s.^3 + v410.*d.s.^4 + v510.*d.s.^5 + ...
            2.*d.tau.*( v020 + v120.*d.s + v220.*d.s.^2 + v320.*d.s.^3 + v420.*d.s.^4) + ...
            3.*d.tau.^2.*(v030 + v130.*d.s + v230.*d.s.^2 + v330.*d.s.^3) + ...
            4.*d.tau.^3.*(v040 + v140.*d.s + v240.*d.s.^2) + ...
            5.*d.tau.^4.*(v050 + v150.*d.s) + ...
            6.*d.tau.^5.*(v060) + ...
            d.pi.*(v011 + v111.*d.s + v211.*d.s.^2 + v311.*d.s.^3 + v411.*d.s.^4) + ...
            2.*d.tau.*d.pi.*(v021 + v121.*d.s + v221.*d.s.^2 + v321.*d.s.^3) + ...
            3.*d.tau.^2.*d.pi.*(v031 + v131.*d.s + v231.*d.s.^2) + ...
            4.*d.tau.^3.*d.pi.*(v041 + v141.*d.s) + ...
            5.*d.tau.^4.*d.pi.*(v051) + ...
            d.pi.^2.*(v012 + v112.*d.s + v212.*d.s.^2 + v312.*d.s.^3) + ...
            2.*d.tau.*d.pi.^2.*(v022 + v122.*d.s + v222.*d.s.^2) + ...
            3.*d.tau.^2.*d.pi.^2.*(v032 + v132.*d.s) + ...
            4.*d.tau.^3.*d.pi.^2.*(v042) + ...
            d.pi.^3.*(v013 + v113.*d.s) + ...
            2.*d.tau.*d.pi.^3.*(v023) + ...
            d.pi.^4.*v014;
    
    drhodCT = drhodv.*dvdtau.*dtaudCT; % eqn A4-2
    
    %% partial derivative of rho with respect to p
    % drhodp = drhodv.*dvdp = drhodv.*dvdpi.*dpidp; % eqn A4-3
    dpidp = 1./d.pu;
    dvdpi = v001 + v101.*d.s + v201.*d.s.^2 + v301.*d.s.^3 + v401.*d.s.^4 + v501.*d.s.^5 + ...
            d.tau.*(v011 + v111.*d.s + v211.*d.s.^2 + v311.*d.s.^3 + v411.*d.s.^4) + ...
            d.tau.^2.*(v021 + v121.*d.s + v221.*d.s.^2 + v321.*d.s.^3) + ...
            d.tau.^3.*(v031 + v131.*d.s + v231.*d.s.^2) + ...
            d.tau.^4.*(v041 + v141.*d.s) + ...
            d.tau.^5.*(v051) + ...
            2.*d.pi.*(v002 + v102.*d.s + v202.*d.s.^2 + v302.*d.s.^3 + v402.*d.s.^4) + ...
            2.*d.pi.*d.tau.*(v012 + v112.*d.s + v212.*d.s.^2 + v312.*d.s.^3) + ...
            2.*d.pi.*d.tau.^2.*(v022 + v122.*d.s + v222.*d.s.^2) + ...
            2.*d.pi.*d.tau.^3.*(v032 + v132.*d.s) + ...
            2.*d.pi.*d.tau.^4.*(v042) + ...
            3.*d.pi.^2.*(v003 + v103.*d.s + v203.*d.s.^2) + ...
            3.*d.pi.^2.*d.tau.*(v013 + v113.*d.s) + ...
            3.*d.pi.^2.*d.tau.^2.*(v023) + ...
            4.*d.pi.^3.*(v004 + v104.*d.s + v014.*d.tau) + ...
            5.*d.pi.^4.*(v005) + ...
            6.*d.pi.^5.*(v006);
          
    drhodp = drhodv.*dvdpi.*dpidp; % eqn A4-3
    
    %% partial derivative of rho with respect to S35A
    % drhodS35A = drhodv.*dvdS35A = drhodv.*dvds.*dsdS35A; eqn A4-4
    dsdS35A   = -4./7.*(sqrt(d.SAu.*(d.SAu + 24)) ./ d.SAu.^2);
    drhodS35A = drhodv.*dvds.*dsdS35A; % eqn A4-4
  end % FUNCTION [drhodSA, drhodCT, drhodp, drhodS35A] = rho_partialdiffyQ

end

