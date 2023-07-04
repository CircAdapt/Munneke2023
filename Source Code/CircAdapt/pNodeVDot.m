function pNodeVDot
% function pNodeVDot
% Pressures in cavities p -> pressure in Nodes: p
%                         -> flows to cavities: VDot
% Includes collapsibility of tubes
% Theo Arts, Maastricht University, Aug 7, 2013

global P

% flow of valves to/from nodes, initialization Node.q/Y/A
nTube  =P.Tube.n;

%=== index conversion sparse matrices
Valve2NodeProx= P.Valve.Valve2NodeProx;
Valve2NodeDist= P.Valve.Valve2NodeDist;
Cavity2Node   = P.Cavity.Cavity2Node;
Tube2NodeProx = P.Tube.Tube2NodeProx;
Tube2NodeDist = P.Tube.Tube2NodeDist;

%=== Valve flow contribution to node pressure
DqValve= P.Valve.q*(-Valve2NodeProx+Valve2NodeDist); %Added flow into node
DAValve= repmat(max(P.Valve.AOpen,P.Valve.ALeak) * ...
    (Valve2NodeProx+Valve2NodeDist),size(P.t)); %Added flow area around node

%==== AV flow through ArtVen
p0Av    = P.ArtVen.p0AV;    % reference AV pressure drop
q0Av    = P.ArtVen.q0AV;    % reference AV flow
kAv     = P.ArtVen.kAV;     % exponent non-linearity AV 'resistance'
iAv2CavP= P.ArtVen.iCavity; % ArtVen ->index-> proximal cavities
iAv2CavD= iAv2CavP+1;       % ArtVen ->index-> distal cavities
iAv2NdP = P.ArtVen.iNode;   % ArtVen -> index -> arterial Node
iAv2NdD = iAv2NdP+1;        % ArtVen -> index -> venous Node
pCav    = P.Cavity.p;       % Cavity pressures
pAvCavP = pCav(:,iAv2CavP); % ArtVen: proximal cavity pressure
pAvCavD = pCav(:,iAv2CavD); % ArtVen: distal   cavity pressure
Dp      = pAvCavP-pAvCavD;  % ArtVen: AV-pressure drop between Cavities

% Correction of ArtVen-flow resistance to volume change of 
% Artven-related cavities, Resistance ~ 1/Volume^2
DpV20     = p0Av.*P.ArtVen.Len.^2./sum(P.ArtVen.A0.^-2);
VCav      = P.Cavity.V;
DpV2      = abs(Dp)./(VCav(:,iAv2CavP).^-2 + VCav(:,iAv2CavD).^-2);
qNorm     = bsxfun(@power,bsxfun(@rdivide,DpV2,DpV20),kAv).*sign(Dp);
qAv       = bsxfun(@times,qNorm,q0Av); % AV flow
P.ArtVen.q= qAv;
Facq      = P.General.FacpControl; % FacpControl= current p / pRef

% Correction of Coronary ArtVen-flow resistance to volume change of
% Coronary ArtVen-related cavities, Resistance ~ 1/Volume^2
p0AvCor   = P.CorArtVen.p0AV;       % reference coronary AV pressure drop
q0AvCor   = P.CorArtVen.q0AV;       % reference coronary AV flow
kAvCor    = P.CorArtVen.kAV;        % exponent non-linearity AV 'resistance'
iAv2CavCP = P.CorArtVen.iCavity;    % CorArtVen ->index-> proximal cavities
iAv2CavCD = iAv2CavCP+1;            % CorArtVen ->index-> distal cavities
iAv2NdCP  = P.CorArtVen.iNode;      % CorArtVen -> index -> arterial Node
iAv2NdCD  = iAv2NdCP+1;             % CorArtVen -> index -> venous Node
pCorP     = pCav(:,iAv2CavCP);      % Proximal coronary Cavity pressure
pCorD     = pCav(:,iAv2CavCD);      % Distal coronary Cavity pressure

%==== External pressure
P = calcIMP(P);
pim = P.CorArtVen.pExt;

%==== AV flow through CorArtVen
% Myocaridal pressure 
nC = P.CorArtVen.n;
pTransMyo = P.CorArtVen.pTransMyo;
pimlog    = repmat(logical([0 1 1 1 1 1 1 0]),1,nC); % only add pim to arterioles and venules
pim2      = zeros(size(pTransMyo));
pim2(:,pimlog)= pim(:,reshape((repmat(repmat([1,2,3],1,2),nC,1)+((1:nC)'-1)*3)',[1,nC*3*2])); % reshape for convenience
pMyo      = pTransMyo + pim2;
P.CorArtVen.pMyo  = pMyo;
pCor              = zeros(size(pMyo,1),10*nC); % Create total matrix of coronary pressures
pCor(:,1:10:end)  = pCorP; % starting with proximal pressure
pCor(:,10:10:end) = pCorD; % ending with distal pressure
pCor(:,repmat(logical([0 1 1 1 1 1 1 1 1 0]),1,nC))=pMyo; % intramyocardial pressures inbetween

% Myocardial volume
VCor        = zeros(size(pMyo,1),10*nC); % Create total matrix of coronary volumes
VCor(:,1:10:end)    = VCav(:,iAv2CavCP); % starting with proximal volume
VCor(:,10:10:end)   = VCav(:,iAv2CavCD); % ending with distal volume
VCor(:,repmat(logical([0 1 1 1 1 1 1 1 1 0]),1,nC))=P.CorArtVen.VMyo; % intramyocardial volumes inbetween

% Coronary flow
Len2      = P.CorArtVen.Len.^2;
if isfield(P.CorArtVen, 'facDil'), facDil = P.CorArtVen.facDil;
else,                      facDil= 1; end
A02       = ((1./facDil).*P.CorArtVen.A0).^-2;
indP      = [1,2,2,2,3,4,5,6,7,8,9]; indP2 = reshape((repmat(indP,nC,1)+((1:nC)'-1)*10)',1,11*nC);
indD      = [2,3,4,5,6,7,8,9,9,9,10];indD2 = reshape((repmat(indD,nC,1)+((1:nC)'-1)*10)',1,11*nC);
Dp        = pCor(:,indP2)-pCor(:,indD2);
A0s       = A02(indP,:)+A02(indD,:);
DpV20     = reshape(p0AvCor.*Len2./A0s,1,11*nC);
VCors     = VCor(:,indP2).^-2 + VCor(:,indD2).^-2;
DpV2      = abs(Dp)./VCors;
kAvCor2    = reshape(repmat(kAvCor',1,11)',1,nC*11);
qNormCor   = bsxfun(@power,bsxfun(@rdivide,DpV2,DpV20),kAvCor2).*sign(Dp);
q0AvCor2   = reshape(repmat(q0AvCor',1,11)',1,nC*11);
q0AvCor2(repmat(logical([0 1 1 1 1 1 1 1 1 1 0]),1,nC)) = q0AvCor2(repmat(logical([0 1 1 1 1 1 1 1 1 1 0]),1,nC)).*repmat([0.95 1 1.05],1,3*P.CorArtVen.n)/3;
qAvCor     = bsxfun(@times,qNormCor,q0AvCor2);  % AV Proximal flow
P.CorArtVen.qAr = qAvCor(:,1:11:end);
P.CorArtVen.qMyo1epi= qAvCor(:,2:11:end); P.CorArtVen.qMyo1mid= qAvCor(:,3:11:end); P.CorArtVen.qMyo1endo= qAvCor(:,4:11:end);
P.CorArtVen.qMyomepi= qAvCor(:,5:11:end); P.CorArtVen.qMyommid= qAvCor(:,6:11:end); P.CorArtVen.qMyomendo= qAvCor(:,7:11:end);
P.CorArtVen.qMyo2epi= qAvCor(:,8:11:end); P.CorArtVen.qMyo2mid= qAvCor(:,9:11:end); P.CorArtVen.qMyo2endo= qAvCor(:,10:11:end);
P.CorArtVen.qVe = qAvCor(:,11:11:end);
P.CorArtVen.q   = (1/5)*(qAvCor(:,1:11:end)+...
                        (qAvCor(:,2:11:end)+qAvCor(:,3:11:end)+qAvCor(:,4:11:end))+...
                        (qAvCor(:,5:11:end)+qAvCor(:,6:11:end)+qAvCor(:,7:11:end))+...
                        (qAvCor(:,8:11:end)+qAvCor(:,9:11:end)+qAvCor(:,10:11:end))+...
                         qAvCor(:,11:11:end));

%===== Cavity contributions to node pressure and flow
iAv2Nd    =[iAv2NdP,iAv2NdD];  % ArtVen -> index -> Node P+D
iAv2NdCor =[iAv2NdCP,iAv2NdCD];% CorArtVen -> index -> Node P+M+D
YCav  = 1./P.Cavity.Z;       % conductivity per cavity
qCav  = pCav .* YCav;        % pressure x conductivity per cavity
DYCav = YCav  * Cavity2Node; % internal Y conductance of flow source
DqCav = qCav  * Cavity2Node; % short circuit flow of node
DqCav(:,iAv2Nd)= DqCav(:,iAv2Nd) + [-qAv*Facq,qAv/Facq];
DqCav(:,iAv2NdCor)= DqCav(:,iAv2NdCor) + [-qAvCor(:,1:11:end)*Facq,qAvCor(:,11:11:end)/Facq];
%      ArtVen flow + volume control
DACav = P.Cavity.A* Cavity2Node; % flow cross-sectional area

%===== Tube contributions to node pressure and flow
YTbProx= 1./P.Tube.ZR; % conductivity per tube
YTbDist= 1./P.Tube.ZL; % conductivity per tube
pTbProx= P.Tube.pSProx; % source pressure
pTbDist= P.Tube.pSDist; % source pressure
YTbS   =[YTbProx,YTbDist] ; % source conductivity
pTbS   =[pTbProx,pTbDist] ; % source pressure
qTbS   = pTbS .* YTbS;      % source flow to node
Tb2Nd  =[Tube2NodeProx;Tube2NodeDist]; % Tube -> matrix -> Node P+D
DYTube = YTbS * Tb2Nd; % added Tube related conductivity
DqTube = qTbS * Tb2Nd; % added Tube related source flow
DATube = P.Tube.A * ( Tube2NodeProx + Tube2NodeDist ); % added area

% No waterfall
YNode=           DYCav + DYTube; % total node conductivity
qNode= DqValve + DqCav + DqTube; % total node inflow with pNode=0
ANode= DAValve + DACav + DATube; % total area, seen from node

pNode= qNode ./ YNode; % node pressure by solving NodeInflow Dq=0

% Detect waterfall conditions in Tube
iTb2NdP= P.Tube.iNodeProx;
iTb2NdD= P.Tube.iNodeDist;
pTbExt = P.Tube.p - P.Tube.pTrans; % external tube pressure
iTb2Nd =[iTb2NdP,iTb2NdD];
pTbNd  = pNode(:,iTb2Nd);
YTbNd  = YNode(:,iTb2Nd);
pTbE   =[pTbExt,pTbExt];

% Waterfall by pressure drop
dpTb = dpWf(pTbNd,pTbE,pTbS,YTbS,YTbNd); %Waterfall pressure drop
DqTb = -YTbS.*dpTb; % waterfall node inflow increment per node
pTbS = pTbS-dpTb; % waterfall pressure correction
qNode= qNode + DqTb * Tb2Nd; % waterfall node inflow correction

YTbProx= YTbS(:, 1:nTube       ); % proximal tube conductivity
YTbDist= YTbS(:,(1:nTube)+nTube); % distal tube conductivity
pTbProx= pTbS(:, 1:nTube       ); % corrected value zero-flow prox pressure
pTbDist= pTbS(:,(1:nTube)+nTube); % corrected value zero-flow dist pressure

% Waterfall in Tube
pNode  = qNode ./ YNode; % Waterfall corrected node pressures

% d/dt Cavity volume
iCav2Nd      = P.Cavity.iNode; % iNode pointed by cavity
P.Cavity.VDot= (pNode(:,iCav2Nd)-P.Cavity.p).*YCav;
P.CorArtVen.VDotMyo = zeros(size(P.CorArtVen.VMyo));
P.CorArtVen.VDotMyo(:,[1:8:end,2:8:end,3:8:end,4:8:end,5:8:end,6:8:end,7:8:end,8:8:end])= [qAvCor(:,1:11:end)-qAvCor(:,2:11:end)-qAvCor(:,3:11:end)-qAvCor(:,4:11:end),qAvCor(:,2:11:end)-qAvCor(:,5:11:end),qAvCor(:,3:11:end)-qAvCor(:,6:11:end),qAvCor(:,4:11:end)-qAvCor(:,7:11:end),qAvCor(:,5:11:end)-qAvCor(:,8:11:end),qAvCor(:,6:11:end)-qAvCor(:,9:11:end),qAvCor(:,7:11:end)-qAvCor(:,10:11:end),qAvCor(:,8:11:end)+qAvCor(:,9:11:end)+qAvCor(:,10:11:end)-qAvCor(:,11:11:end)];

% Tube inflow and outflow
pNodeProx= pNode(:,iTb2NdP); %pressure Prox-node
pNodeDist= pNode(:,iTb2NdD); %pressure Dist-node
qP=(pNodeProx-pTbProx).*YTbProx; % prox tube flow
qD=(pTbDist-pNodeDist).*YTbDist; % dist tube flow

P.Tube.qProx= qP; % Tube inflow
P.Tube.qDist= qD; % Tube outflow
P.Tube.VDot = qP-qD; % Tube volume change
P.Node.p = pNode; % Node pressure
P.Node.Y = YNode; % Total Node conductivity
P.Node.q = qNode; % zero pressure node inflow
P.Node.A = ANode; % total cross-section of connections to node

end

function dp=dpWf(pN,pE,pS,YS,YN)
% {pN,pE,pS}= pressure {Node, External, Source}
% {Y,YN}= conductivity {waterfall, Node total}
% dp= decrease of source pressure to simulate waterfall
%=====================
eps= 0.2;
PS = pS-pE;
PN = pN-pE;
OutFlow= double(PS>PN); % With inflow no waterfall
% p  = OutFlow.*max(0,min(0,PS)-PN); %Waterfall
p  = OutFlow.*max(0,min(0,PS)-PN); %Waterfall
x  = 1-YS./YN; % effect of node impedance <-> tube impedance
z = x+eps*exp(-(x/eps).^2); %safety to avoid zero division
dp= p./z;
end
