function PatchWallA2T
% function PatchWallA2T
% Determines linear function patch area Am = fu Tension T
% Am(T)=Am0+T*DADT
% Constants Am0 and DADT are determined from Patch state variables
% and reference geometry
% Theo Arts, Maastricht University, Aug 31, 2014 // Mar 15, 2022

global P;

% Patch Am= Am0+DADT*T, provides Am0 and DADT
AmRef = P.Patch.AmRef; % midwall area for ref sarcomere length 2mu
LsRef = P.Patch.LsRef; % ref sarcomere length 2mu
VWall = P.Patch.VWall; % wall volume
LSe   = P.Patch.LenSeriesElement; % series elastic length at isometric contraction
Lsi   = P.Patch.Lsi; % unloaded sarc length= state variable

Lambda= bsxfun(@rdivide,Lsi      ,LsRef); %extension factor
Am    = bsxfun(@times  ,Lambda.^2,AmRef); %actual midwall area
Ef    = log(Lambda); %natural fiber strain

P.Patch.Ef  = Ef; % fiber strain Ef with zero length SE-element

SarcEf2Sf; % sarcomere strain->stress

Sf    = P.Patch.Sf; % sarcomere stress
DSfDEf= P.Patch.DSfDEf;

DEf   = Sf./DSfDEf; % - zero tension strain relative to Lsi
Corr  = 1+1.16*(2*LSe./Lsi+DEf); % correction factor for DADT
% Value 1.16 indicates middle of working range (2.0~isometric contraction)
DADT  = 4*Corr .* Am.^2 ./ bsxfun(@times,DSfDEf,VWall);
Am0   = Am .* exp(-2*DEf); % zero wall tension area
P.Patch.DADT= DADT; % area compliance
P.Patch.Am0 = Am0; % zero tension area

% Wall is composed of patches: Also for wall: Am(T)=Am0+DADT*T
for iWall=1:P.Wall.n
    iPatch= (P.Wall.iPatch(iWall)-1)+(1:P.Wall.nPatch(iWall));
    P.Wall.VWall(iWall) =sum(P.Patch.VWall(iPatch));
    P.Wall.Am0(:,iWall) =P.Wall.AmDead(iWall)+sum(P.Patch.Am0(:,iPatch),2);
    P.Wall.DADT(:,iWall)=sum(P.Patch.DADT(:,iPatch),2);
end
end

