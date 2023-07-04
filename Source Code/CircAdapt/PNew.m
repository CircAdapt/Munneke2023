function PNew
% function PNew
%
% Sets backbone of structure P
% Attribute names to ArtVen, TriSeg, Chambers, Valves, Tubes, 
% Walls and Patches
% Connections between elements are defined by strings stored in
% structure field P.Net
%
% Structure is built on information in P.Net
% ArtVen represents artery-peripheral resistance-vein of organ or body part
% Chamber represents a cavity enclosed by a myocardial wall (atria)
% TriSeg represents combination of two cavities with three walls
% Bag represents passive elastic bag, encapsulating part of the
% circulation, like the pericardium
% Node: named connection point
% Cavity: volume with elastic or muscular wall, connected to a node
% Wall: muscular wall can contract, composed of 1 or more patches
% Patch: contracting part of a wall, having specific mechanical properties
% Valve: valve with inertia connects proximal to distal node, may leak
% Tube: elastic wave-guiding tube connects proximal to distal node
%
% Theo Arts, Maastricht University, April 3, 2013

global P;
P=CreateP; % sets the tree structure with all necessary fields

%======= DEFINITION OF STRUCTURE BY ELEMENTS AND CONNECTIONS ==========
% Load excel sheet containing arterial tree properties
P.General.excelTopology = 'CorTerDat.xlsx';
[CorDatArt,CorTubArt] = xlsread(P.General.excelTopology,'Arterial');
[CorDatVen,CorTubVen] = xlsread(P.General.excelTopology,'Venous');
[~,CorArtVenName]  = xlsread(P.General.excelTopology,'CorArtVen');

Net.ArtVen.Name = {'Ca','Br','Ce','Fe','Pu'}; % carotid, brachial, celiac, femoral + Pulmonary 
Net.CorArtVen.Name = CorArtVenName(2:end,2)'; % coronary branches LCx (5), LAD (7), RCA (6)
Net.Chamber.Name= {'La','Ra'}; %atria
Net.TriSeg.Name = {'v'}; %ventricular unit

% Valve and Tube connections
Net.Valve.Nodes=[... % proximal and distal node of valves
    {'Cs'  ,'Ra'  };...
    {'Vc'  ,'Ra'  };...
    {'Ra'  ,'Rv'  };...
    {'Rv'  ,'PuAr'};...
    {'PuVe','La'  };...
    {'La'  ,'Lv'  };...
    {'Lv'  ,'Ao'  }];

Net.Tube.Nodes=[... % proximal and distal node of elastic tubes
    CorTubArt(2:end,2:3);...
    {'Ao'  ,'BrAr'};...
    {'BrAr','CaAr'};...
    {'BrAr','CeAr'};...
    {'CeAr','FeAr'};...
    CorTubVen(2:end,2:3);...
    {'Vc'  ,'BrVe'};...
    {'BrVe','CaVe'};...
    {'BrVe','CeVe'};...
    {'CeVe','FeVe'}];

% Bags pressurize parts of the circulation
Net.Bag.Name   ={'Peri'     ,'Thorax'}; %pericardium, thorax
Net.Bag.Chamber={{'La','Ra'},{}      }; % enclosed chambers
Net.Bag.TriSeg ={{'v'}      ,{}      }; % enclosed TriSeg
Net.Bag.ArtVen ={{}         ,{'Pu'}  }; % enclosed ArtVen
Net.Bag.CorArtVen={CorArtVenName(2:end,2)',{}}; % enclosed CorArtVen
Net.Bag.Tube   ={[CorTubArt(2:end,4);CorTubVen(2:end,4)]',{'AoBrAr','BrArCeAr','VcBrVe','BrVeCeVe'}};
Net.Bag.Bag    ={{}         ,{'Peri'}}; % pericardium inside thorax

Net.Node.Baro ={'BrAr'}; % pressure control node
Net.Wall.MultiPatch={}; % defines Wall's split in Patch's

% Transfer of structure information to P-structure
P.Net = Net;
Indexation; % mutual element relations expressed by indices

%==========================================
%============== Data filling ==============
%==========================================

% Species specific general information
P.General.q0           = 85e-6; % mean systemic flow
P.General.p0           = 12200; % mean systemic pressure
P.General.tCycle       = 0.85 ; % cycle time

P.General.RhoB         = 1050  ; % blood density
P.General.EtaB         = 0.004 ; % blood viscosity
P.General.TauAv        = 0.16  ; % AV depolarization delay
P.General.SaturationControl= 0  ;
P.General.AdaptFeedback    = 0.3; % feedback strength
P.General.nStoredBeats     = 3; % maximum number of stored beats in Ps

% Crude estimates of maximum pressures
pLv= 15000  ; % peak Lv pressure
pRv= 3800   ; % peak Rv pressure
pLa= 1400   ; % mean pLa
pRa= 560    ; % mean pRa
pPu= 1900   ; % mean pulmonary artery pressure
p0 = P.General.p0;

% General settings
P.General.tCycleRest = P.General.tCycle;
P.General.FacpControl= 1    ;
P.General.Dt         = 0.001;
P.General.ScaleVqY   = [1e-5,1e-4,1e-1];
P.General.tEnd       = 1.0;
P.General.TimeFac    = 1  ;
P.General.Fast       = 0;
P.General.AdaptFunction = 'Adapt0P';

%======= Adaptation setpoints

AdaptationParameters
SarcomereProperties
MakeHeart(pLv,pRv,pLa,pRa)
ArtVenParameters(p0,pPu,pLa,pRa)
TubeParameters
BagParameters
ValveParameters

P.t=0;
P2SVar;

save P P

end

% ================= Auxilary functions ====================
function AdaptationParameters

global P;

% ArtVen.Adapt and Tube.Adapt
Put({'ArtVen','Adapt'},'WallStress','All',[1;1]*500e3);
Put({'ArtVen','Adapt'},'vFlowMean' ,'All',[1;1]*0.17 );
Put({'ArtVen','Adapt'},'vImpact'   ,'All',[1;1]*3.0  );
Put({'CorArtVen','Adapt'},'WallStress','All',[1;1]*500e3);
Put({'CorArtVen','Adapt'},'vFlowMean' ,'All',[1;1]*0.17 );
Put({'CorArtVen','Adapt'},'vImpact'   ,'All',[1;1]*3.0  );
Put({'Tube'  ,'Adapt'},'WallStress','All',500e3);
Put({'Tube'  ,'Adapt'},'vFlowMean' ,'All',0.17 );
Put({'Tube'  ,'Adapt'},'vImpact'   ,'All',3.0  );

% Patch/Sarcomere
PatchA={};
for iP=1:P.Patch.n
    if P.Patch.Name{iP}(2)=='a'
        PatchA=[PatchA,P.Patch.Name(iP)];
    end
end
% Ventricle=default
Put({'Patch','Adapt'},'SfPasMax','All' ,6000 );
Put({'Patch','Adapt'},'SfPasAct','All' ,4800 );
Put({'Patch','Adapt'},'FacSfAct','All' ,0.61 );
Put({'Patch','Adapt'},'LsPasAct','All' ,2.1  );
Put({'Patch','Adapt'},'SfPasMax',PatchA,30000);
Put({'Patch','Adapt'},'FacSfAct',PatchA,0.35 );
end

function SarcomereProperties
global P;
PatchA={};
for iP=1:P.Patch.n
    if P.Patch.Name{iP}(2)=='a'
        PatchA=[PatchA,P.Patch.Name(iP)];
    end
end

% Ventricular: default
Put('Patch','Lsi'             ,'All' , 1.9      );
Put('Patch','C'               ,'All' , 0.001    );
Put('Patch','Depolarization'  ,'All' , [0;0]-P.General.tCycle);
Put('Patch','LsRef'           ,'All' , 2.0      );
Put('Patch','Ls0Pas'          ,'All' , 1.8      );
Put('Patch','dLsPas'          ,'All' , 0.6      );
Put('Patch','SfPas'           ,'All' , 22000    );
Put('Patch','SfPas'           ,'All' , [32711 29141 23510 22861 25431]    );
Put('Patch','Lsi0Act'         ,'All' , 1.51     );
Put('Patch','LenSeriesElement','All' , 0.04     );
Put('Patch','SfAct'           ,'All' , 84000    );% lowered P114
Put('Patch','vMax'            ,'All' , 7.0      );
Put('Patch','TR'              ,'All' , 0.25     );
Put('Patch','TD'              ,'All' , 0.25     );
Put('Patch','TimeAct'         ,'All' , 0.42     );
Put('Patch','CRest'           ,'All' , 0.0      );
Put('Patch','dT'              ,'All' , 0.0      );

%Atrial: non-default
% Put('Patch','SfPas'           ,PatchA, 50000    );
Put('Patch','SfAct'           ,PatchA, 59000    );% lowered P114
Put('Patch','TimeAct'         ,PatchA, 0.15     );
Put('Patch','TR'              ,PatchA, 0.4      );
Put('Patch','TD'              ,PatchA, 0.4      );

end
%=============================

function MakeHeart(pLv,pRv,pLa,pRa)
global P;
PatchA={'La1','Ra1'};
PatchV={'Lv1','Sv1','Rv1'};

% Geometry of heart walls
VStroke = P.General.q0*P.General.tCycle;
A0      = (580*VStroke^2)^(1/3);
Put('Patch','VWall',PatchA, VStroke*[18.968,5.5541  ].*[pLa,pRa    ]./Get('Patch','SfAct',PatchA));
Put('Patch','VWall',PatchV, VStroke*[8.4048,3.4590,11.9667].*[pLv,pLv,pRv]./Get('Patch','SfAct',PatchV));
Put('Patch','AmRef',PatchA, [0.5519,0.2839]*A0);
Put('Patch','AmRef',PatchV, [0.6564,0.3282,0.8338]*A0);

% Cavity starting values of volume state variables V
% Heart
Heart={'La','Ra','Lv','Rv'};
Put('Cavity','V','All', 0); %initialization
Put('Cavity','V',Heart, VStroke*[1.06 0.73 1.69 1.52]);
%TriSeg
Put('TriSeg','V','v',0.52*VStroke      );
Put('TriSeg','Y','v',0.80*VStroke^(1/3));
P.TriSeg.Tau = 2.5*P.General.Dt; % lowpass for V and Y pre-estimate

SplitMergePNewTube('Lv1',12);
SplitMergePNewTube('Sv1',5);
P.Net.Wall.MultiPatch={'Lv12','Sv5'}; % defines Wall's split in Patch's
end

function ArtVenParameters(p0,pPu,pLa,pRa)
global P
%---% ArtVen
Sy = {'Ca','Br','Ce','Fe'};
Pu = 'Pu';
q0 = P.General.q0;

Put('ArtVen','q0AV' ,Sy, [0.15,0.12,0.57,0.16]*.95*q0);
Put('ArtVen','q0AV' ,Pu, q0);
%   peripheral flow distribution over ArtVen vessel beds
Put('ArtVen','kAV'  ,Sy, 1);
Put('ArtVen','kAV'  ,Pu, 2);
Put('ArtVen','p0'   ,Sy,[p0 ; pRa]);
Put('ArtVen','p0'   ,Pu,[pPu; pLa]);
Put('ArtVen','k'    ,'All',[12;12]); %estimate
Put('ArtVen','Len'  ,'All',6.0*P.ArtVen.q0AV.^(1/3)); %estimate
Put('ArtVen','Len'  ,{'Br','Fe'},[0.3,0.5]); %arm and leg are long
P.ArtVen.A0=([1;1]*P.ArtVen.q0AV)./P.ArtVen.Adapt.vFlowMean;
P.ArtVen.AWall= diag([0.20;0.10])*P.ArtVen.A0;
P.ArtVen.p0AV = P.ArtVen.p0(1,:)-P.ArtVen.p0(2,:);
P.ArtVen.p0AV(1:4) = [9.3420    9.0946    9.1436    9.4357]*1e3;

% Cavities of Artven
iCavArt= P.ArtVen.iCavity;
iCavVen= iCavArt+1;
V      = P.ArtVen.A0 .* P.ArtVen.Len([1,1],:);
P.Cavity.V(iCavArt)=V(1,:);
P.Cavity.V(iCavVen)=V(2,:);

%---% CorArtVen 
% Retrieve data from Excel file
[~,CorArtVenName]  = xlsread(P.General.excelTopology,'CorArtVen');
Cor = CorArtVenName(2:end,2)';   % Coronary names
CorTer = CorArtVenName(2:end,3); % Coronary territories

% Link coronary territory to Patch
Cor2Patch = zeros(numel(CorTer),P.Patch.n); 
for i = 1:numel(CorTer)
    Cor2Patch(i,:) = strcmp(CorTer(i),P.Patch.Name)';
end
P.CorArtVen.Cor2Patch = Cor2Patch;

% Input data
VWall2Cor = sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)';
VWallT    = sum(VWall2Cor);
Put({'CorArtVen'},'q0AV',Cor,VWall2Cor./VWallT*(0.05*P.General.q0)); %Based on volume fraction of Patch
Put({'CorArtVen'},'kAV' ,Cor,1);
Put({'CorArtVen'},'p0'  ,Cor,[12200; 10600 ; 4667; 4667; 4667; 1600; 1600; 1600; 1060; 560]);%[p0 ; (p0-pRa)/4*3+pRa ; (p0-pRa)/2+pRa ; (p0-pRa)/4+pRa ; pRa]);%[p0 ; 35/7.5e-3 ; 20/7.5e-3 ; 12/7.5e-3 ; pRa]);%[p0 ; (p0-pRa)/4*3+pRa ; (p0-pRa)/2+pRa ; (p0-pRa)/4+pRa ; pRa]);
Put({'CorArtVen'},'k'   ,Cor,1.5*[16;16;16;16;16;16;16;16;16;16]);
Put({'CorArtVen'},'Len' ,Cor,6.0*P.CorArtVen.q0AV.^(1/3)); %estimate
P.CorArtVen.A0=([1.4;1.4;2.0;2.0;2.0;6.125;6.125;6.125;9.8;9.8]*P.CorArtVen.q0AV)./repmat(P.CorArtVen.Adapt.vFlowMean,5,1);
P.CorArtVen.AWall= diag([0.20;0.20;0.10;0.10;0.10;0.10;0.10;0.10;0.10;0.10])*P.CorArtVen.A0;
P.CorArtVen.p0AV = [P.CorArtVen.p0(1,:)-P.CorArtVen.p0(2,:);...
                    P.CorArtVen.p0(2,:)-P.CorArtVen.p0(3,:);...
                    P.CorArtVen.p0(2,:)-P.CorArtVen.p0(4,:);...
                    P.CorArtVen.p0(2,:)-P.CorArtVen.p0(5,:);...
                    P.CorArtVen.p0(3,:)-P.CorArtVen.p0(6,:);...
                    P.CorArtVen.p0(4,:)-P.CorArtVen.p0(7,:);...
                    P.CorArtVen.p0(5,:)-P.CorArtVen.p0(8,:);...
                    P.CorArtVen.p0(6,:)-P.CorArtVen.p0(9,:);...
                    P.CorArtVen.p0(7,:)-P.CorArtVen.p0(9,:);...
                    P.CorArtVen.p0(8,:)-P.CorArtVen.p0(9,:);...
                    P.CorArtVen.p0(9,:)-P.CorArtVen.p0(10,:)];
P.CorArtVen.FC   = 1; % Coronary flow regulation (0: off, 1: on (normal), 2: on (baseline), 3: wall mass adaptation)
P.CorArtVen.a    = 1; % Scales CEP component of IMP
P.CorArtVen.y    = 1; % Scales VE  component of IMP
P.CorArtVen.pimSf= 0.060; % Conversion factor for fiber stress to VE pressure component (~20% of max pLV)
P.CorArtVen.facq0AV = ones(1,P.CorArtVen.n); % Parameter for flow regulation, scales coronary q0AV
P.CorArtVen.facDil = ones(size(P.CorArtVen.A0)); % Parameter for flow regulation, determines amount of vasodilation

% Cavities of Artven.Cor
iCavArtCor      = P.CorArtVen.iCavity;
iCavVenCor      = iCavArtCor+1;
V      = P.CorArtVen.A0 .* P.CorArtVen.Len([1,1,1,1,1,1,1,1,1,1],:);
P.Cavity.V(iCavArtCor)  = V(1,:);
P.Cavity.V(iCavVenCor)  = V(10,:);
VMyo   = V(2:9,:); VMyo = VMyo(:)';
P.CorArtVen.VMyo      = VMyo;

end

function TubeParameters
global P
% Flow distribution calculated from ArtVen flows and shunt flows
FlowDistribution;  % mean flow distribution in ArtVen, Tubes, Valves
P.Tube.A0=abs(P.Tube.q./P.Tube.Adapt.vFlowMean); % Tube cross-section

% Sequence TubeNames:
[CorDatArt,CorTubArt] = xlsread(P.General.excelTopology,'Arterial');
[CorDatVen,CorTubVen] = xlsread(P.General.excelTopology,'Venous');
Art= [CorTubArt(2:end,4)',{'AoBrAr','BrArCaAr','BrArCeAr','CeArFeAr'}];
Ven= [CorTubVen(2:end,4)',{'VcBrVe','BrVeCaVe','BrVeCeVe','CeVeFeVe'}];
All=[Art,Ven];
p0 = P.General.p0; % arterial pressure
q0 = P.General.q0;
Aux= Get('ArtVen','p0','Br');
pRa= Aux(2); % venous pressure

% pressure distribution in tubes
Put('Tube','p0',Art,p0 );
Put('Tube','p0',Ven,pRa);
% additional tube properties
Put('Tube','k'      ,All, [CorDatArt(:,6)',8,11,11,15,CorDatVen(:,6)',8,11,11,15]);
Put('Tube','Len'    ,All, [CorDatArt(:,5)',7,12,23,15,CorDatVen(:,5)',7,12,23,15]/100);
Put('Tube','AWall'  ,All, P.Tube.A0.* ...
    (12*P.Tube.p0./P.Tube.Adapt.WallStress+P.Tube.Adapt.vImpact*0.02));

% Volume and flow
Put('Tube','V','All', P.Tube.Len.*P.Tube.A0);
Put('Tube','q','All', P.Tube.A0.*P.Tube.Adapt.vFlowMean);

% Memory allocation for Wave delay lines, needed for Tube function
Row=zeros(1 ,P.Tube.n);
P.Tube.pL = P.Tube.p0;% Left wave, non-delayed, circular storage
P.Tube.pR = P.Tube.p0;
u0=sqrt(P.Tube.p0*P.General.q0);
P.Tube.uP = u0;% Delayed and attenuated pressure signal
P.Tube.uD = u0;
P.Tube.TauL= Row+P.General.Dt;% delay
P.Tube.TauR= Row+P.General.Dt;

end

function BagParameters
global P;
Heart={'La','Ra','Lv','Rv'};
% Pericardium, Thorax Bags
P.Bag.k     = [10,10];
P.Bag.pAdapt= [200,50]; % transmural bag pressure pressure
%Estmate heart volume
VWall = sum(Get('Patch','VWall','All'));
VCav  = sum(Get('Cavity','V',Heart));
VHeart= VCav+VWall;
P.Bag.VRef= [1.4*VHeart,2.7*VHeart]; % thorax volume not really used
end

function ValveParameters

A=Get('Tube','A0','AoBrAr'); % Aortic cross-section

Put('Valve','q'    ,'All',0.0);
Put('Valve','AOpen','All',A  );
Put('Valve','ALeak','All',A*1e-6);
Put('Valve','Len'  ,'All',sqrt(A));
% Mitral and Tricuspid valve are larger
Put('Valve','AOpen',{'RaRv','LaLv'},...
    1.5* Get('Valve','AOpen',{'RaRv','LaLv'}) );
% vene-atrial orifices are inertias without valve leaflets
Put('Valve','ALeak',{'VcRa','PuVeLa'},...
    Get('Valve','AOpen',{'VcRa','PuVeLa'}) );
% L-R shunting valves are closed
Put('Valve','ALeak','CsRa',Get('Valve','AOpen','CsRa'));

% Wall: AmDead (non-contractile area)
ALv=sum(Get('Valve','AOpen',{'LvAo'  ,'LaLv'        }));
ARv=sum(Get('Valve','AOpen',{'RvPuAr','RaRv'        }));
ALa=sum(Get('Valve','AOpen',{'PuVeLa','LaLv'        }));
ARa=sum(Get('Valve','AOpen',{'VcRa'  ,'RaRv', 'CsRa'}));
Put('Wall','AmDead',{'Lv','Rv','La','Ra'},[ALv,ARv,ALa,ARa]);
end

function FlowDistribution
% calculates flow distribution through ArtVens, Valves and Tubes
% Flow ~=0 are assumed to be known. All other flows are unknown. If
% a sufficient number of flows is known, equations on steady state flow
% distribution are solved in a least squares sense

global P

nNode=P.Node.n;
nValve=P.Valve.n;
nTube=P.Tube.n;
nArtVen=P.ArtVen.n;
nCorArtVen=P.CorArtVen.n;

Ma =zeros(nNode,nArtVen);
Mc =zeros(nNode,nCorArtVen);
Mt =zeros(nNode,nTube);
Mv =zeros(nNode,nValve);

for ia=1:nArtVen
    iP=P.ArtVen.iNode(ia);
    iD=iP+1;
    Ma(iP,ia)=Ma(iP,ia)-1;
    Ma(iD,ia)=Ma(iD,ia)+1;
end
for ia=1:nCorArtVen
    iP=P.CorArtVen.iNode(ia);
    iD=iP+1;
    Mc(iP,ia)=Mc(iP,ia)-1;
    Mc(iD,ia)=Mc(iD,ia)+1;
end
for ia=1:nTube
    iP=P.Tube.iNodeProx(ia);
    iD=P.Tube.iNodeDist(ia);
    Mt(iP,ia)=Mt(iP,ia)-1;
    Mt(iD,ia)=Mt(iD,ia)+1;
end
for ia=1:nValve
    iP=P.Valve.iNodeProx(ia);
    iD=P.Valve.iNodeDist(ia);
    Mv(iP,ia)=Mv(iP,ia)-1;
    Mv(iD,ia)=Mv(iD,ia)+1;
end
M=[Ma,Mc,Mt,Mv];
Rga= 1:nArtVen; Rgc= nArtVen+(1:nCorArtVen); Rgt= nArtVen+nCorArtVen+(1:nTube); Rgv= nArtVen+nCorArtVen+nTube+(1:nValve);
nM = size(M,2);
q=[P.ArtVen.q0AV,P.CorArtVen.q0AV,P.Tube.q(1,:),P.Valve.q(1,:)];
Rg0=find(q~=0); %known flows q
Rg1=setdiff(1:nM,Rg0); % unknown flows q
q0=q(Rg0);
q1=-pinv(M(:,Rg1))*M(:,Rg0)*q0';
q(Rg1)=q1;
P.ArtVen.q0AV=q(Rga);
P.CorArtVen.q0AV=q(Rgc);
P.Tube.q=q(Rgt);
P.Valve.q=q(Rgv);
end


function P=CreateP
% Creates empty structure P with fields
P=[];

FieldsP={
    'General'
    'ArtVen'
    'CorArtVen'
    'Chamber'
    'TriSeg'
    'Valve'
    'Tube'
    'Node'
    'Cavity'
    'Wall'
    'Patch'
    'Bag'
    'Net'
    'SVar'
    'SVarDot'
    't'
    'tDot'
    };

FieldsGeneral={
    'q0'
    'p0'
    'tCycle'
    'tCycleRest'
    'DtSimulation'
    'tStart'
    'tEnd'
    'Dt'
    'ScaleVqY'
    'FacpControl'
    'TimeFac'
    'TauAv'
    'rhob'
    'AdaptFunction'
    'Fast'
    'In'
    'Out'
    };

FieldsArtVen={
    'Name'
    'n'
    'iCavity'
    'iNode'
    'k'
    'Len'
    'A0'
    'p0'
    'AWall'
    'p0AV'
    'q0AV'
    'kAV'
    'q'
    'qProx'
    'qDist'
    'Adapt'
    'ArtVen2NodeArt'
    'ArtVen2NodeVen'
    };

FieldsCorArtVen={
    'Name'
    'n'
    'iCavity'
    'iNode'
    'k'
    'Len'
    'A0'
    'p0'
    'AWall'
    'p0AV'
    'q0AV'
    'kAV'
    'q'
    'qProx'
    'qDist'
    'qAr'
    'qMA'
    'qMV'
    'qVe'
    'VDotMyo'
    'VMyo'
    'AMyo'
    'pTransMyo'
    'pMyo'
    'pExt'
    'Adapt'
    'ArtVen2NodeArt'
    'ArtVen2NodeVen'
    };

FieldsChamber={
    'Name'
    'n'
    'iCavity'
    'iWall'
    };

FieldsTriSeg={
    'Name'
    'n'
    'iCavity'
    'iWall'
    'V'
    'Y'
    'VDot'
    'YDot'
    'VS'
    'YS'
    'Tau'
    };

FieldsCavity={
    'Name'
    'n'
    'iNode'
    'V'
    'VDot'
    'A'
    'Z'
    'p'
    'pTrans'
    'Cavity2Node'
    };

FieldsWall={
    'Name'
    'n'
    'nPatch'
    'iPatch'
    'VWall'
    'Am0'
    'AmDead'
    'DADT'
    'T'
    'Cm'
    'Am'
    'pTrans'
    };

FieldsPatch={
    'Name'
    'n'
    'Lsi'
    'C'
    'LsiDot'
    'CDot'
    'Depolarization'
    'LsRef'
    'SfPas'
    'Ls0Pas'
    'dLsPas'
    'SfAct'
    'Lsi0Act'
    'LenSeriesElement'
    'vMax'
    'TimeAct'
    'TR'
    'TD'
    'CRest'
    'VWall'
    'AmRef'
    'Ef'
    'Ls'
    'SfEcm'
    'SfTit'
    'Sf'
    'DSfDEf'
    'T'
    'DADT'
    'Am0'
    'Am'
    'TauRefrac'
    'DepPath'
    'iPace'
    'dT'
    'Adapt'
    };

FieldsNode={
    'Name'
    'n'
    'iBaro'
    'q'
    'p'
    'Y'
    'A'
    };

FieldsBag={
    'Name'
    'n'
    'iCavity'
    'iWall'
    'iTube'
    'iBag'
    'OK'
    'VRef'
    'k'
    'pAdapt'
    'V'
    'pTrans'
    'p'
    'VC'
    'VW'
    'VT'
    'VB'
    };

FieldsValve={
    'Name'
    'n'
    'iNodeProx'
    'iNodeDist'
    'AOpen'
    'ALeak'
    'Len'
    'q'
    'qDot'
    'AvValves'
    'AvWalls'
    'Valve2NodeProx'
    'Valve2NodeDist'
    };

FieldsTube={
    'Name'
    'n'
    'iNodeProx'
    'iNodeDist'
    'k'
    'Len'
    'TauAtt'
    'A0'
    'p0'
    'AWall'
    'V'
    'VDot'
    'q'
    'A'
    'pTrans'
    'p'
    'pSProx'
    'pSDist'
    'qProx'
    'qDist'
    'ZR'
    'ZL'
    'pL'
    'pR'
    'uP'
    'uD'
    'TauL'
    'TauR'
    'cL'
    'cR'
    'Tube2NodeProx'
    'Tube2NodeDist'
    'Adapt'
    };

FieldsNet={
    'ArtVen'
    'Chamber'
    'TriSeg'
    'Node'
    'Valve'
    'Tube'
    'Bag'
    'Wall'
    'Depolarization'
    };

Empty10=ones(1,0);
for i=1:length(FieldsP)
    fld=FieldsP{i};
    P.(fld)=[];
    FldStr=['Fields',fld];
    if exist(FldStr)
        SubFldStr=eval(FldStr);
        for j=1:length(SubFldStr)
            P.(fld).(SubFldStr{j})=Empty10;
        end
    end
end

end

