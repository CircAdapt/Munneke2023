function OutDotT=SVarDot(SVarT)
% function OutDotT=SVarDot(SVarT)
% SVarT = column vector= transpose of state variables SVar vector
% OutDotT= column vector of derivatives
% State Variables -> Time derivatives of State Variables
% Ready to be used in MatLab function odeCA() to solve set of Diff Eq
% Theo Arts, Maastricht University, Jun 8, 2018
%====

global P

P.SVar= SVarT'; % store state variables SVar
SVar2P; % state variables SVar -> physiologic representation P.xx
P.tDot=ones(size(P.t)); % time derivative of time = 1

MemAlloc     ; % necessary memory allocations
ArtVenV2p    ; % arteries to veins compliant network
CorArtVenV2p ; % arteries to veins compliant network
PatchWallA2T ; % patch and wall: Am= Am0+T*DADT

% Start iter for more accurate calculation
for iter = 1:4  
    
    ChamberV2p   ; % Chamber, pTrans and Wall.Am,T
    TriSegV2p    ; % TriSeg, pTrans and Wall.Am,T
    Wall2Patch2Sarc; % filling Sarc with LsiDot,CDot

    % Recalculate patch with true Ef
    % T is known
    Am = P.Patch.Am0 + P.Patch.T.*P.Patch.DADT;
    Ef = 0.5*log(P.Patch.Am./P.Patch.AmRef);
    P.Patch.Am = Am;
    P.Patch.Ef = Ef;

    SarcEf2Sf; % sarcomere strain->stress

    Sf    = P.Patch.Sf; % sarcomere stress
    DSfDEf= P.Patch.DSfDEf;

    DEf   = Sf./DSfDEf; % - zero tension strain relative to Lsi
    DADT  = 4 .* Am.^2 ./ bsxfun(@times,DSfDEf,P.Patch.VWall);
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
    
    % To determine amount of iterations, determine error below
    % if numel(P.t) > 1
    %     if exist('YS')
    %         disp(max(abs(P.TriSeg.YS)./abs(YS)));
    %         disp(max(abs(P.TriSeg.VS)./abs(VS)));
    %         disp(max(abs(P.Wall.T)./abs(T)));
    %     end
    %     YS = P.TriSeg.YS;
    %     VS = P.TriSeg.VS;
    %     T = P.Wall.T;
    % end

% End iter
end

TubeV2p      ; % Delayed source transmural pressures and resistances of tube
pCavity      ; % Determines p in cavities by adding bag pressures
pNodeVDot    ; % Node pressures and volume derivatives VDot's
ValveqDot    ; % flow derivatives qDot
TubeDelays   ; % renders wave amplitude just at moment of reflection, delays

P2SVarDot    ; % transfer of derivatives to P-structure

OutDotT= real(P.SVarDot)'; % odeXX requires OutDotT to be a column vector
end

