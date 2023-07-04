function Adapt
%function Adapt
% Common to all adaptation procedures.
% Embeds more dedicated adaptation procedures
%    executed between beats:
% Blood pressure is controlled by change of circulating volume
% Flow is set by P.General.SaturationControl:
%     0-> target flow per ArtVen-element
% or: 1-> flow is determined by oxygen usage= AV-Dsaturation x flow
% Theo Arts, Maastricht University, April 26, 2014

global P
save PTemp P; %saves last intermediate solution

% Systemic blood flow adjusted per ArtVen-element by resistance change
% Control of systemic pressure by adjustment of circulatory blood volume
FbFac=P.General.AdaptFeedback; % Feedback constant. If FbFac==0, no control

% Finding systemic and pulmonary ArtVen's
iPu = contains(P.ArtVen.Name,'pu','IgnoreCase',true); %all indexes with 'pu' are pulmonary
iSy = ~iPu; 

%=== Estimation of TargetFlow in systemic ArtVen's for flow control
qSy   = P.ArtVen.q(:,iSy);% simulated flows in systemic ArtVen
qPu   = P.ArtVen.q(:,iPu);% simulated flows in pulmonary ArtVen
qCorAr= P.CorArtVen.qAr(:,:);% simulated flows in coronary ArtVen
q0AvSy= P.ArtVen.q0AV(:,iSy);%Get reference flow per systemic ArtVen
q0AvCor=P.CorArtVen.q0AV(:,:);%Get reference flow per coronary ArtVen

% direct flow control, ratio of flows in systemic ArtVens preserved
if isfield(P.ArtVen,'q0AVt')
    qSyTarget=q0AvSy/sum(q0AvSy)*P.ArtVen.q0AVt;
    P.ArtVen.q0AV(iSy) = qSyTarget(1:numel(q0AvSy));
else
    facq0AV = P.CorArtVen.facq0AV;
    qSyTarget=q0AvSy/sum(q0AvSy)*(P.General.q0-sum(facq0AV.*q0AvCor)); % Target flow in systemic and coronary ArtVen's
    P.ArtVen.q0AV(iSy) = qSyTarget(1:numel(q0AvSy));
end
% save Sy-Target flows in P.ArtVen.q0AV

%Control of pressure by change of circulating volume
FacpControl= mean(P.Node.p(:,P.Node.iBaro))/P.General.p0;
P.General.FacpControl=FacpControl^FbFac; % >1 -> Volume decrease

% Control of flow to P.ArtVen.q0AV for systemic ArtVen's
p0AvSy = P.ArtVen.p0AV(iSy);
kAvSy  = P.ArtVen.kAV(iSy);
FacqControlSy = exp(FbFac*(1-mean(qSy)./q0AvSy));
p0AvSyTarget = p0AvSy.*FacqControlSy.^(-1./kAvSy)/P.General.FacpControl;
P.ArtVen.p0AV(iSy) = p0AvSyTarget;

% Control of flow to P.CorArtVen for coronary ArtVen's
% 0 means no flow control
% 1 should be used for flow regulation from the reference
% 2 should be used to establish the baseline simulation (as it includes 
% coronary adaptation to the current situation)
% 3 should be used to adapt the reference simulation (in case of LBBB)
if P.CorArtVen.FC == 1 
    % Coronary Flow Regulation - Normal
    P = CorFC(P);
elseif P.CorArtVen.FC == 2 
    % Coronary Flow Regulation - For Baseline Simulation
    P = CorFCbase(P);
elseif P.CorArtVen.FC == 3 
    % Coronary Flow Regulation - For LBBB Wall Adaptation
    P = CorFCadap(P);
end

% Flow in Sys and Pu to judge steady state
FlowVec=[sum([mean(qSy),mean(qCorAr)]),sum(mean(qPu))]; % Systemic+Coronary/Pu flow
disp('Flow/q0: Sys+Cor Pu');
disp(FlowVec/P.General.q0);
if P.CorArtVen.FC == 1 || P.CorArtVen.FC == 2  || P.CorArtVen.FC == 3
    disp('Demand/Demand0:');
    disp(P.CorArtVen.SSA2Cor./P.CorArtVen.SSA2Cor0);
    disp(P.CorArtVen.VO2./P.CorArtVen.VO20);
end

% Estimate AV-delay
P.General.TauAv = 0.185*P.General.tCycle;

% ====== AdaptSpecific, different for Rest, Exercise
NoAdapt= strcmp(P.General.AdaptFunction,'Adapt');
if NoAdapt % No adaptation, regular sequence of beats
    % === Faster Steady State at rest
    VecV=P.Cavity.V;
    P.General.In =[P.General.In ;VecV(  1,:) ];
    P.General.Out=[P.General.Out;VecV(end,:) ];
    if P.General.Fast;% Put new start values in structure P
        Vec= SteadyStateP;
        i=1  ; j=P.Cavity.n; P.Cavity.V=Vec(i:j);
    end
else
    feval(P.General.AdaptFunction); % specific adapt function
end

% Judging quality of steady state
if size(P.General.Out,1)>1; % Escape if steady state is reached
    ErrVec= 1000*log( P.General.Out(end,:)./P.General.In(end,:) );
    
    StatErr = sqrt(mean(ErrVec.^2));
    
    disp(['Stationarity error= ',...
        num2str( round(sqrt(mean(ErrVec .^ 2 ))) )] );
    %=== ERROR criterium on flow stationarity
    if sqrt(mean(ErrVec .^ 2 ))<1 && (~NoAdapt || P.General.Fast)
        P.General.tEnd= P.General.tEnd-0.5*(P.General.tEnd-P.t(end));
    end
end;

% get the initial condition for next beat, P is most compact information to
% start the simulation
P2SVar; % load physiologic data in record of state variables P.SVar

end

