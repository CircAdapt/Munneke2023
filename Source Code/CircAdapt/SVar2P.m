function SVar2P
% function SVar2P
% Transfer of scaled state variables in P.SVar to P-structure (SI-units)
% Inverse of P2SVar
% Theo Arts, Maastricht University, Jun 5, 2018

global P
% Scalings
qRef= P.General.q0;
VRef= qRef*P.General.tCycle;

nCav   =P.Cavity.n; % number of cavities
nCor   =8*P.CorArtVen.n; % number of CorArtVen secret cavities (8x Myo)
nTube  =P.Tube.n; % number of cavities
nValve =P.Valve.n; % number of valves
nPatch =P.Patch.n; % number of patches with representative sarcomere
nTriSeg=P.TriSeg.n; % number of TriSeg's
% finding indices for value transfer
a=cumsum([0,1,nCav,nCor,nTube,nValve,nPatch,nPatch,nTriSeg,nTriSeg]);
iB=a(1:end-1)+1; iE=a(2:end); % successive begin and end indices

% value transfer
% Conversion of volumes, flows and distances to conveniently scaled
% variables
i=1;
P.t        =             P.SVar(:,iB(i):iE(i))  ; i=i+1;
P.Cavity.V =       exp(  P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.CorArtVen.VMyo=  exp(  P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.Tube.V   =       exp(  P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.Valve.q  = qRef* sinh( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.Patch.C  =             P.SVar(:,iB(i):iE(i))  ; i=i+1;
P.Patch.Lsi=             P.SVar(:,iB(i):iE(i))  ; i=i+1;
P.TriSeg.V = VRef* sinh( P.SVar(:,iB(i):iE(i)) ); i=i+1;
P.TriSeg.Y =       exp(  P.SVar(:,iB(i):iE(i)) ); i=i+1;

end

