function TubeV2p
% function TubeV2p
% Tube volume V -> transmural pressure pTrans, wave propagation velocity
% c0, proximal and distal zero-flow pressure pSProx, pSDist and source
% impedances ZR (prox), ZL (dist), using delayed pressures.
% Theo Arts, Maastricht University, June 29, 2018

global P

dt    = P.General.dt   ; % integration timestep
it    = round(P.t/dt)+1; % time sample counter
Len   = P.Tube.Len     ; % representative length of blood vessels
p0    = P.Tube.p0      ; % working pressure
A0    = P.Tube.A0      ; % tube cross-section at p0
Aw    = P.Tube.AWall   ; % wall cross-section
k     = P.Tube.k       ; % tube stiffness coefficient
RhoB  = P.General.RhoB ; % blood density
fWom  = 2.5/P.General.tCycleRest;% Womersley frequency in A0-flow signal
EtaB  = P.General.EtaB ; %blood viscosity

% Transmural pressure and wave velocity/impedance= fu(cross-section)
% Lo represents negative transmural pressure component for A/AWall<1.0
A    = bsxfun(@rdivide,P.Tube.V,Len); % vessel cross-section
eps  = 0.1; % clipped to minimum A/AWall value
AN   = max(bsxfun(@rdivide,A,Aw),eps); % norm to wall cross-section
A0N  = A0./Aw; % norm of reference A0 to wall cross-section
A    = bsxfun(@times,AN,Aw); % area clipped to always positive
hRPois= bsxfun(@rdivide,4*pi*EtaB*Len, A.^2);
m    = k/3-1; % stifness exponent
Lo   = max(0,1./AN-1); % if Lo>0, negative compression pressure occurs
pHiN = bsxfun(@power,bsxfun(@rdivide,AN+0.5,A0N+0.5),m);% pressure Hi
pLoN = -Lo.^2; % norm pressure component related to compression
pTrans= bsxfun(@times,p0,pHiN + pLoN); % transmural pressure
dpHi  = bsxfun(@rdivide,bsxfun(@times,pHiN,m),(AN+0.5));
dpLo  = 2*Lo./AN.^2;
dpdAN = bsxfun(@times,p0,dpHi+dpLo); 
dpdA  = dpdAN./Aw; % compliance
c0    = sqrt(0.25+dpdA.*A/RhoB); % zero flow wave velocity

%Approximation of Womersley attenuation with aWom= r Sqrt(rho w/eta)
aWom    = sqrt(A*(2*fWom*RhoB/EtaB)); % Womersley number
c0      = c0./(1+0.71./aWom+0.71./aWom.^4); % Womersley correction
Att     = fWom*6.3*sqrt(1+0.03125*aWom.^2)./(1+0.25*aWom.^2); % attenuation
idxCor  = contains(P.Tube.Name,{'LM','LAD','LCx','RCA','Cs','GCV','VAL','MCV'}); % Find coronary tubes
Att(:,idxCor) = 4*(0.85/P.General.tCycle)*Att(:,idxCor); % Necessary for stability
Z0      = RhoB * c0 ./ A; % corrected wave impedance with tube flow=0

% Flow dependency of wave velocity and impedance
% getting delayed signals for zero-flow pressure pSProx and pSDist
q     = P.Tube.q(it,:); % tube flow
vb    = q./A; % mean blood velocity
hvdc0 = 0.5*vb./c0; % blood velocity, normalized to wave velocity
% Wave velocity and wave impedance
b     = sqrt(1+hvdc0.^2);
bR    = b+hvdc0;
bL    = b-hvdc0;
cR    = sqrt(0.25+(c0.*bR).^2);% wave velocity clipped to >0.5 m/s
cL    = sqrt(0.25+(c0.*bL).^2);% wave velocity clipped to >0.5 m/s
ZR    = Z0./bR ;
ZL    = Z0./bL;
DpPois= q.*hRPois; % effect DC 

P.Tube.pSProx = P.Tube.uP(it,:).*sqrt(ZR)+DpPois; % zero flow pressure Prox
P.Tube.pSDist = P.Tube.uD(it,:).*sqrt(ZL)-DpPois; % zero flow pressure Dist
P.Tube.ZR     = ZR     ; %R-wave impedance
P.Tube.ZL     = ZL     ; %L-wave impedance
P.Tube.cR     = cR     ; %R-wave velocity
P.Tube.cL     = cL     ; %L-wave velocity
P.Tube.Att    = Att    ; %Womersley attenuation factor
P.Tube.A      = A      ; %Cross-sectional area
P.Tube.pTrans = pTrans ; %Transmural pressure
end
