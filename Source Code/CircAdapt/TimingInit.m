function TimingInit
% function TimingInit
% Sets patch depolarization times for upcoming beat
% Theo Arts, Maastricht University, Oct 24, 2018

global P

G        = P.General    ;
tCycleRef= G.tCycleRest ; % Reference tCycle for body size scaling
tCycle   = G.tCycle     ; % current cycle time
TimeFac  = G.TimeFac    ; % scaling of contraction time
TauAv    = G.TauAv      ; % AV-delay, controlled in Adapt0

%=== setting TimeAct [s] duration of contraction for Ls==LsStress0Act
% in atrial and ventricular wall patches
ta= 0.15*(tCycle/0.85)*TimeFac; % atrial activation duration
tv= (0.10*tCycleRef +0.40*tCycle)*TimeFac; % ventricular " "
% Find atrial patches, discrimination atria<->ventricles
iPa=[]; iPv=[];
for i=1:P.Patch.n
    if regexp(P.Patch.Name{i},'a')==2
        iPa=[iPa,i];
    else
        iPv=[iPv,i];
    end
end
P.Patch.TimeAct(iPa)=ta;
P.Patch.TimeAct(iPv)=tv;
%================ end setting TimeAct====================

%====================================================
% definition of interpatch conduction pathways
dTauAv = 0;%P.General.dTauAv;
P.Net.Depolarization.Pathway={...
    'Ra1','Ra1',tCycle;...
    'Ra1','La1',0.02 * (tCycleRef/0.85) * TimeFac;...
    'Ra1','Rv1',TauAv+dTauAv;...
    'Ra1','Lv1',TauAv+dTauAv;...
    'Rv1','Sv1',0;...
    'Sv1','Lv1',0;...
    };
P.Net.Depolarization.Pace         = 'Ra1'; % leading pacemaker patch
P.Net.Depolarization.TauInterPatch= 0; % default interpatch delay
P.Net.Depolarization.TauRefrac    = 0.25 ; % default refractory period

% === CONDUCTION PATHWAYS as defined in P.Net.Depolarization
Pathway= P.Net.Depolarization.Pathway;
nP     = length(Pathway); %number of paths
DepPath= zeros(nP,3); % initialization of [iPatchProx,iPatchDist,Tau]
for i=1:nP; % depolarization pathways copied from P.Net.Depolarization
    Path        = Pathway(i,:);
    DepPath(i,:)= [Get('Patch','Index',Path(1:2)),Path{3}];
end
% Intrawall interpatch pathways
DepP    = [];
Tau     = P.Net.Depolarization.TauInterPatch; %Default interpatch delay
for iW  = 1:P.Wall.n
    iP  = P.Wall.iPatch(iW);
    nP  = P.Wall.nPatch(iW);
    Aux = (0:nP-2)';
    dM  = [iP+Aux,(iP+1)+Aux];
    DepP=[DepP;dM;fliplr(dM)];
end
P.Patch.DepPath=[ DepPath; [DepP,repmat(Tau,[size(DepP,1),1])] ]';
InterPatchDelay = P.Patch.dT(P.Patch.DepPath(2,:)) - P.Patch.dT(P.Patch.DepPath(1,:));
P.Patch.DepPath(3,:) = P.Patch.DepPath(3,:) + InterPatchDelay;
% indexed depolarization pathways
P.Patch.TauRefrac=repmat(P.Net.Depolarization.TauRefrac,[1,P.Patch.n]);
% default refractory period
P.Patch.iPace= Get('Patch','Index',P.Net.Depolarization.Pace);
% Pacemaker patch
% === END indexation CONDUCTION PATHWAYS ================

% Shifting/preparing Tube-delayed signals
nt  = ceil(tCycle/P.General.Dt); % number of time points upcoming beat
dt  = tCycle/(nt-1); % set integer number of time steps per cycle
P.General.dt= dt; %slightly changed dt to get integer number of t-steps

t1   = P.t; % previous beat
pL1  = P.Tube.pL; % Wave signals previous beat may enter current beat
pR1  = P.Tube.pR;
uP1  = P.Tube.uP;
uD1  = P.Tube.uD;
q1   = P.Tube.q ;

t2   = (0:nt-1)'*dt; % time samples next beat
t12  = t1-t1(end)+tCycle; % time shift previous->current beat

% Delayed signals interpolated in time,relevant with change of tCycle or Dt
P.Tube.pL= interp0(t12,pL1,t2);%left  wave used as delay line
P.Tube.pR= interp0(t12,pR1,t2);%right wave used as delay line
P.Tube.uP= interp0(t12,uP1,t2);%proximal zero flow pressure
P.Tube.uD= interp0(t12,uD1,t2);%distal zero flow pressure
P.Tube.q = interp0(t12,q1 ,t2);%mean tube flow
P.Tube.dt= dt; %sampling interval Tube delay lines
% end Tube signal preparation

end

function p2=interp0(t1,p1,t2)
% converts time scale of previous beat to that of current beat
% if interpolation out of interval, boundary point is used
if numel(t1)==1
    t1=[0;t1];
    p1=p1([1,1],:);
end
t1(1)=min(0,t1(1));
p2=interp1(t1,p1,t2);
end
