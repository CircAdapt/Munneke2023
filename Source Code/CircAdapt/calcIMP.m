function P = calcIMP(P)
% External pressure (pim) was assumed to be a summation of two mechanisms:
% (1) VE:  Myocardial stiffness, labelled as Varying Elastance
% (2) CEP: Transmission of cavity pressure into the myocardial wall, termed
%          Cavity-Induced Extracellular Pressure
% Input:  P struct
% Output: pimVE  / P.CorArtVen.pimVE  (VE component)
%         pimCEP / P.CorArtVen.pimCEP (CEP component)
%         pim    / P.CorArtVen.pExt   (total intramyocardial pressure)

% Cor2Patch matrix
Cor2Patch = P.CorArtVen.Cor2Patch;

% Distribution based on VWall
Cor2VWall = (Cor2Patch.*P.Patch.VWall) ./ sum(Cor2Patch.*P.Patch.VWall,2);

%---% Calculation of VE component
VE2Patch    = zeros(numel(P.Patch.Sf(:,3)),P.Patch.n);
VE          = zeros(numel(P.Patch.Sf(:,3)),P.CorArtVen.n); 
pimSf       = P.CorArtVen.pimSf; 
for iS = 1:numel(P.Patch.Sf(:,3))
    VEmat          = Cor2Patch*pimSf.*(P.Patch.Sf(iS,:));
    VE2Patch(iS,:) = sum(VEmat);
    VE(iS,:)       = sum(VEmat.*Cor2VWall,2);
end
pimVE = P.CorArtVen.y*VE;
P.CorArtVen.pimVE = pimVE;
P.CorArtVen.VE2Patch = VE2Patch;

%---% Calculation of CEP component
% Gather pressure components
patm    = P.Bag.p(:,strcmp(P.Bag.Name, 'Peri'));
pLV     = P.Cavity.p(:,strcmp(P.Cavity.Name, 'Lv'));
pRV     = P.Cavity.p(:,strcmp(P.Cavity.Name, 'Rv'));

% Find which patches correspond to LV, S, RV
LVpatch = startsWith(P.Patch.Name,'Lv');
Spatch  = startsWith(P.Patch.Name,'Sv');
RVpatch = startsWith(P.Patch.Name,'Rv');

% Define pressure at endocardium and epicardium
% For the LV/RV : endocardium equals LV/RV cavity pressure
%               : epicardium equals pericardial pressure
% For the S     : endocardium equals LV cavity pressure
%               : epicardium equals RV cavity pressure
CEPendo2Patch = zeros(numel(pLV),P.Patch.n);
CEPepi2Patch  = zeros(numel(pLV),P.Patch.n); 
CEPendo       = zeros(numel(pLV),P.CorArtVen.n); 
CEPepi        = zeros(numel(pLV),P.CorArtVen.n); 
for iP = 1:numel(pLV)
    % Endocardium
    CEP2endo = Cor2Patch.*(LVpatch*pLV(iP)) + Cor2Patch.*(RVpatch*pRV(iP)); %LV en RV patches
    CEP2endo = CEP2endo + Cor2Patch.*(Spatch*pLV(iP)); %Septal patches
    % Epicardium
    CEP2epi  = Cor2Patch.*(LVpatch*patm(iP)) + Cor2Patch.*(RVpatch*patm(iP)); %LV en RV patches
    CEP2epi  = CEP2epi + Cor2Patch.*(Spatch*pRV(iP));
    % CEP at endocardial and epicardial outer borders
    CEPendo2Patch(iP,:) = sum(CEP2endo);
    CEPepi2Patch(iP,:) = sum(CEP2epi);
    CEPendo(iP,:) = sum(CEP2endo.*Cor2VWall,2);
    CEPepi(iP,:) = sum(CEP2epi.*Cor2VWall,2);
end
indCorC = repmat(1:P.CorArtVen.n,3,1); indCorC = indCorC(:)';
CEPendo2= CEPendo(:,indCorC); 
CEPepi2 = CEPepi(:,indCorC); 
pimfac  = P.CorArtVen.pimfac;
pimfac  = pimfac(:)';
CEP     = pimfac.*CEPendo2 + (1-pimfac).*CEPepi2;
pimCEP  = P.CorArtVen.a*CEP;
P.CorArtVen.pimCEP = pimCEP;
P.CorArtVen.CEPendo2Patch = CEPendo2Patch;
P.CorArtVen.CEPepi2Patch = CEPepi2Patch;

%---% Total intramyocardial pressure
pim     = pimCEP + pimVE(:,indCorC);
P.CorArtVen.pExt = pim;

end