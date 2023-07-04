function ChamberV2p
% function ChamberV2p
% A chamber is a cavity encapsulated in a single myocardial wall.
% Calculates: Cavity volume V ->
% myofiber stress Sf, wall tension T and cavity pressure p
% using the locally linearized T(Am) relation
% Theo Arts, Maastricht University, Oct 13, 2012

global P;
if P.Chamber.n==0 %if there is no chamber
    return
end

Chamber=P.Chamber; % Chamber structure

RhoB    = P.General.RhoB;
iCavity = Chamber.iCavity; % index of related cavity
iWall   = Chamber.iWall  ; % index of related wall
V    = max(0,P.Cavity.V(:,iCavity)); %cavity volumes
VWall= P.Wall.VWall(iWall); % wall volumes
nt   = size(V,1); % nt=number of time samples,
Vm   = max(0,V) + repmat(VWall/2,[nt,1]);% numerical safety with V<0 
Cm   = (4/3*pi./Vm).^(1/3); % mid wall curvature
Am   = (4*pi)./Cm.^2; % midwall area
Am0  = P.Wall.Am0(:,iWall); % zero tension midwall area
DADT = P.Wall.DADT(:,iWall); % wall compliance
p0   = 0.1*P.General.p0;

% Am=max(Am,Am0);% buckling with T<0
a    = 100; % buckling parameter, decrease of stiffness with buckling
DADT = DADT.*exp(a*max(0,1-Am./Am0));%Buckling reduces stiffness
T    = (Am-Am0)./DADT; % wall tension

% wall properties
P.Wall.T(:,iWall)     = T; % wall tension
P.Wall.Cm(:,iWall)    = Cm; % curvature=1/radius
P.Wall.Am(:,iWall)    = Am; % wall area
pTrans= 2*Cm.*T; % transmural pressure with effect of buckling
P.Wall.pTrans(:,iWall)= pTrans; % transmural pressure

% Anti-collapse pressure
eps   = 0.1; % lowest value of VNLo=VCavity/VWall for pressure calculation
VNLo  = max(bsxfun(@rdivide,V,VWall),eps);
dpLo  = bsxfun(@times,p0,max(0,0.5./VNLo-1).^2);% anti-collapse safety
pTrans= pTrans - dpLo;

% Cavity impedance properties, needed to make node connection
Len= 2*Vm.^(1/3); % cavity length
A  = Vm ./ Len; % cross-sectional area
% A  = ( V + 0.1*repmat(VWall,[nt,1]) ) ./Len; % cross-sectional area
% Z0 = sqrt(rhob*Len./abs(A.*DADT)); %Cavity wave impedance
Z0 = sqrt(RhoB./(abs(A).^1.5 .* DADT)); %Cavity wave impedance term
% Z0 = sqrt(4.0*RhoB./abs(A.*DADT.*Len)); %Cavity wave impedance++++
P.Cavity.A(:,iCavity) = A; % cross-sectional area for valve 
P.Cavity.Z(:,iCavity) = 0.5*Z0; % made quite small
P.Cavity.pTrans(:,iCavity)= pTrans;

end

