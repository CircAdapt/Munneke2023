function Wall2Patch2Sarc
% function Wall2Patch2Sarc
% After multiple patch geometry has been solved, for walls and patches
% stresses and areas are calculated for output. Needed to calculates
% derivatives of sarcomere state variables.
% Theo Arts, Maastricht University, Oct 13, 2012

global P;

for iWall = 1:P.Wall.n % wall tension T is transferred to the patches
    nPatch= P.Wall.nPatch(iWall);
    iPatch= P.Wall.iPatch(iWall)+(0:nPatch-1);
    P.Patch.T(:,iPatch)= P.Wall.T(:,repmat(iWall,[1,nPatch])); %each patch
end

Am        = P.Patch.Am0 + P.Patch.T .* P.Patch.DADT; % Patch area
P.Patch.Am= Am;
P.Patch.Ef= 0.5*log(bsxfun(@rdivide,Am,P.Patch.AmRef)); % Sarcomere strain

SarcEf2Sf % sarcomere mechanics
end

