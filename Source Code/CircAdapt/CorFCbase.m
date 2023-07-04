function P = CorFCbase(P)
    
    % Systemic blood flow adjusted per ArtVen-element by resistance change
    % Control of systemic pressure by adjustment of circulatory blood volume
    FbFac=P.General.AdaptFeedback; % Feedback constant. If FbFac==0, no control
    
    %----% Estimation of TargetFlow in coronary ArtVen's for flow control
    iSy    = ~contains(P.ArtVen.Name,'pu','IgnoreCase',true); % {'Ca','Br','Ce','Fe'} systemic ArtVens
    Cor    = P.CorArtVen.Name;      % Coronary ArtVens
    q0AvSy = P.ArtVen.q0AV(:,iSy);  % Get reference flow per systemic ArtVen
    q0AvCor= P.CorArtVen.q0AV;      % Get reference flow per coronary ArtVen
    qCorM1P= P.CorArtVen.qMyo1epi;  % simulated epi flows in coronary ArtVen
    qCorM1M= P.CorArtVen.qMyo1mid;  % simulated mid flows in coronary ArtVen
    qCorM1N= P.CorArtVen.qMyo1endo; % simulated endo flows in coronary ArtVen
    
    % Coronary Autoregulation for baseline simulation
    %----% First, calculate the coronary flow error
    % Coronary flow error
    mqCor = [mean(P.CorArtVen.qAr      );...
             mean(P.CorArtVen.qMyo1epi );...
             mean(P.CorArtVen.qMyo1mid );...
             mean(P.CorArtVen.qMyo1endo);...
             mean(P.CorArtVen.qMyomepi );...
             mean(P.CorArtVen.qMyommid );...
             mean(P.CorArtVen.qMyomendo);...
             mean(P.CorArtVen.qMyo2epi );...
             mean(P.CorArtVen.qMyo2mid );...
             mean(P.CorArtVen.qMyo2endo);...
             mean(P.CorArtVen.qVe      )];
    if isfield(P.CorArtVen,'mqCori'), mqCori = P.CorArtVen.mqCori;
    else,                             mqCori = zeros(size(mqCor)); end
    mqCorD = abs(mqCori-mqCor);
    P.CorArtVen.mqCori = mqCor;
    
    %----% Second, Update Demand
    %      Demand is calculated by means of sarcomere stress - strain area
    %      For safety, only update demand if coronary flow error is small
    Ef      = P.Patch.Ef;                     % [-],  Fiber Strain
    Lsi0Act = P.Patch.Lsi0Act;                % [um], Zero-stress Active Length
    LsRef   = P.Patch.LsRef;                  % [um], Reference Fiber Length
    Lf      = bsxfun(@times,exp(Ef),LsRef);   % [um], Fiber Length
    Sf      = P.Patch.Sf;                     % [Pa], Total Fiber Stress 
    Sfp     = P.Patch.SfEcm+P.Patch.SfTit;    % [Pa], Passive Fiber Stress 
    Sfa     = Sf - Sfp;                       % [Pa], Active Fiber Stress
    VWall   = P.Patch.VWall;                  % [m3], Myocardial Wall Volume
    AmRef   = P.Patch.AmRef;                  % [m2], Reference Midwall Area 
    Am      = P.Patch.Am;                     % [m2], Midwall Area        
    % Determine point of maximal contraction ("end-systole")
    [~,ES]  = max(P.Patch.C);
    Efes    = diag(Ef(ES,:))';
    Sfes    = diag(Sfa(ES,:))';
    Ames    = diag(Am(ES,:))';
    Tes     = Sfes .* VWall ./ (2 * Ames);
    % Determine onset of contraction ("end-diastole")
    ED      = round(P.Patch.Depolarization(1,:)/P.General.Dt);
    ED(ED==0) = 1; % Safety
    Efed    = diag(Ef(ED,:))';
    Sfed    = diag(Sfa(ED,:))';
    % Start point
    Am0     = ((Lsi0Act./LsRef).^2).* AmRef;
    T0      = 0 .* Am0;
    Ef0     = 0.5 .* log(Am0./AmRef);
    Sf0     = 0 .* Ef0;
    % Sarcomere (active) Stress - Strain Area
    for ip = 3:P.Patch.n
        % Linear "end-systolic"-active-Tension-Area relation
        ATrA(:,ip) = linspace(Ames(ip), Am0(ip), 1000)';
        ATrT(:,ip) = linspace( Tes(ip),  T0(ip), 1000)';
        % Determine "end-systolic"-active-Stress-Strain relation   
        Efr(:,ip)  = 0.5 .* log(ATrA(:,ip) ./ AmRef(ip));
        Sfr(:,ip)  = 2 * ATrA(:,ip) .* ATrT(:,ip) ./ VWall(ip);
        % Sarcomere (active) stress-strain area (J)
        SSAr(1,ip)    = -trapz([Ef0(ip); Efes(ip);  Ef(ED(ip):ES(ip),ip); Efr(:,ip)],...
                               [Sf0(ip);  Sf0(ip); Sfa(ED(ip):ES(ip),ip); Sfr(:,ip)])*VWall(ip);
    end
    Cor2Patch = P.CorArtVen.Cor2Patch;
    SSA2Cor    = sum(SSAr.*Cor2Patch,2)';
    k1         = 4.94e-3;   % [mmol/J]
    k2         = 24.2;      % [mmol/m^3/s]
    k2cor      = k2*sum(Cor2Patch.*VWall,2)'; 
%     VO2        = k1*SSLAt2Cor/P.General.tCycle + k2;  % [mmol/m^3/s]], Oxygen consumption 
    VO2        = k1*SSA2Cor/P.General.tCycle + k2cor;  % [mmol/s]], Oxygen consumption 
    P.CorArtVen.VO2 = VO2;
    P.CorArtVen.VWall0 = sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)';
    
%     if max(mqCorD(:)) < 0.1e-9
        %-% Update Demand
        q0AVt      = VO2/sum(VO2)*(0.04*P.General.q0);
        P.CorArtVen.q0AV        = q0AVt;     % Target flow
        P.CorArtVen.facq0AV     = ones(size(1:P.CorArtVen.n)); % At baseline no flow deviation
        P.CorArtVen.VO20        = VO2;       % Oxygen consumption per coronary
        P.Patch.SSAr0           = SSAr;      % Total energy per patch
        P.CorArtVen.SSA2Cor0    = SSA2Cor;   % Total energy per coronary
        P.Patch.SSAr            = SSAr;      % Total energy per patch
        P.CorArtVen.SSA2Cor     = SSA2Cor;   % Total energy per coronary
%     else
%         P.CorArtVen.facq0AV     = ones(size(1:P.CorArtVen.n)); % At baseline no flow deviation
%         P.Patch.SSAr            = SSAr;      % Total energy per patch
%         P.CorArtVen.SSA2Cor     = SSA2Cor;   % Total energy per coronary
%     end
    
    %----% Third, CorArtVen and Tube input parameters are based on flow distribution
%     if max(mqCorD(:)) > 0.1e-9
        % Length
        P.CorArtVen.Len = 6.0*P.CorArtVen.q0AV.^(1/3); %estimate

        % Cross-section
        A0 = P.CorArtVen.A0;
        A0t=([1.4;1.4;2;2;2;6.125;6.125;6.125;9.8;9.8]*P.CorArtVen.q0AV)./repmat(P.CorArtVen.Adapt.vFlowMean,5,1);
        P.CorArtVen.A0 = A0./exp(0.3*(1-A0t./A0));
        
        % Pressure
        P.CorArtVen.p0 = repmat([12200; 10600 ; 4667; 4667; 4667; 1600; 1600; 1600; 1060; 560],1,P.CorArtVen.n);

        % Wall area
        P.CorArtVen.AWall= diag([0.20;0.20;0.10;0.10;0.10;0.10;0.10;0.10;0.10;0.10])*P.CorArtVen.A0;

        % Tube Stifness
        P.Tube.k([1:29,34:63]) = 20;
        P.Tube.k([1:29,34:63]) = 26;

        % Tube cross-section
        q0AV = P.CorArtVen.q0AV;
        P.Tube.A0(1:29)  = abs([sum(q0AV( 1:12)),...
                                sum(q0AV( 1:5 )),...
                                sum(q0AV( 1   )),...
                                sum(q0AV( 2:3 )),...
                                sum(q0AV( 2   )),...
                                sum(q0AV( 3   )),...
                                sum(q0AV( 4:5 )),...
                                sum(q0AV( 4   )),...
                                sum(q0AV( 5   )),...
                                sum(q0AV( 6:12)),...
                                sum(q0AV( 6   )),...
                                sum(q0AV( 7   )),...
                                sum(q0AV( 8:12)),...
                                sum(q0AV( 8   )),...
                                sum(q0AV( 9   )),...
                                sum(q0AV(10:12)),...
                                sum(q0AV(10   )),...
                                sum(q0AV(11   )),...
                                sum(q0AV(12   )),...
                                sum(q0AV(13:18)),...
                                sum(q0AV(18   )),...
                                sum(q0AV(13:17)),...
                                sum(q0AV(13   )),...
                                sum(q0AV(14:17)),...
                                sum(q0AV(14   )),...
                                sum(q0AV(15   )),...
                                sum(q0AV(16:17)),...
                                sum(q0AV(16   )),...
                                sum(q0AV(17   ))]./P.Tube.Adapt.vFlowMean(1:29)); % Tube cross-section
        P.Tube.A0(34:63) = abs([sum(q0AV( 1:12)),...
                                sum(q0AV( 4:5 )),...
                                sum(q0AV( 4   )),...
                                sum(q0AV( 5   )),...
                                sum([q0AV( 1:3 ), q0AV( 6:12)]),...
                                sum(q0AV( 1   )),...
                                sum(q0AV( 2:3 )),...
                                sum(q0AV( 2   )),...
                                sum(q0AV( 3   )),...
                                sum(q0AV( 6:12)),...
                                sum(q0AV( 6   )),...
                                sum(q0AV( 7   )),...
                                sum(q0AV( 8:12)),...
                                sum(q0AV( 8   )),...
                                sum(q0AV( 9   )),...
                                sum(q0AV(10:12)),...
                                sum(q0AV(10   )),...
                                sum(q0AV(11   )),...
                                sum(q0AV(12   )),...
                                sum(q0AV(13:18)),...
                                sum(q0AV(13:17)),...
                                sum(q0AV(13   )),...
                                sum(q0AV(14:17)),...
                                sum(q0AV(14   )),...
                                sum(q0AV(15   )),...
                                sum(q0AV(16:17)),...
                                sum(q0AV(16   )),...
                                sum(q0AV(17   )),...
                                sum(q0AV(18   )),...
                                sum(q0AV(18   ))]./P.Tube.Adapt.vFlowMean(34:63)); % Tube cross-section
                            
         P.Tube.AWall = P.Tube.A0.*(12*P.Tube.p0./P.Tube.Adapt.WallStress+P.Tube.Adapt.vImpact*0.02);

%     end
    
    %----% Fourth, change p0AV to account for correct resistance ratio
    %      Takes into account volume correction
    %      p0AV is equal for epi, mid and endo layers
    p0AvCor= P.CorArtVen.p0AV; p0AvCor = p0AvCor(3,:); %[mid]
    kAvCor = P.CorArtVen.kAV;
    FacqControlCor = exp(FbFac*(1-mean(qCorM1M)./(q0AvCor/3)));
    p0AvCorTarget = p0AvCor.*FacqControlCor.^(-1./kAvCor);%/P.General.FacpControl;
    % Calculate Resistance
    nC   = P.CorArtVen.n;
    p0AV = P.CorArtVen.p0AV; p0AV([2,3,4],:) = repmat(p0AvCorTarget,3,1);
    q0AV = P.CorArtVen.q0AV;
    q0AVs= repmat(q0AV',1,11)'.*[1 1/3 1/3 1/3 1/3 1/3 1/3 1/3 1/3 1/3 1]';
    q0AVs= q0AVs.*repmat([1 repmat([0.95 1 1.05],1,3) 1]',1,P.CorArtVen.n);
    Rim0 = p0AV./q0AVs;
    % Change resistance to account for intramyocardial ratio ([0.60 : 0.30 : 0.10])
    Rim0t = Rim0;
    Rim0t(6:11:end) = Rim0t(3:11:end).*(0.3/0.6);
    Rim0t(9:11:end) = Rim0t(3:11:end).*(0.1/0.6);
    p0AV([5,6,7],:) = repmat((Rim0t(6:11:end)./Rim0(6:11:end)).*p0AV(6,:),3,1);
    p0AV([8,9,10],:) = repmat((Rim0t(9:11:end)./Rim0(9:11:end)).*p0AV(9,:),3,1);
    % Calculate resistance again to account for overall ratio (0.28 : [0.60 : 0.30 : 0.10] : 0.07)
    Rim0 = p0AV./q0AVs;
    Rpp0 = Rim0(2:11:end) + Rim0(5:11:end) + Rim0(8:11:end);
    Rpm0 = Rim0(3:11:end) + Rim0(6:11:end) + Rim0(9:11:end);
    Rpn0 = Rim0(4:11:end) + Rim0(7:11:end) + Rim0(10:11:end);
    Rp0  = 1./(1./Rpp0 + 1./Rpm0 + 1./Rpn0);
    Rim0t = Rim0;
    Rim0t(1:11:end) = Rp0.*(0.28/0.65);
    Rim0t(11:11:end) = Rp0.*(0.07/0.65);
    p0AV(1,:) = (Rim0t(1:11:end)./Rim0(1:11:end)).*p0AV(1,:);
    p0AV(11,:) = (Rim0t(11:11:end)./Rim0(11:11:end)).*p0AV(11,:);
    P.CorArtVen.p0AV = p0AV;
    
    %----% Fifth, Autoregulation by means of vasodilation/contriction (facDil)
    %      Occurs mainly in the arterioles
    kAvCor  = P.CorArtVen.kAV;
    facDil0 = P.CorArtVen.facDil;
    facDil  = [facDil0(3,:) facDil0(4,:) facDil0(5,:)];
    FacqControlCor = exp(0.3*(1-mean([qCorM1P, qCorM1M, qCorM1N])./([0.95*q0AvCor q0AvCor 1.05*q0AvCor]/3)));
    facDilTarget   = facDil./FacqControlCor.^(-1./[kAvCor kAvCor kAvCor])/P.General.FacpControl;
    facDilMax      = [repmat(1.70,1,P.CorArtVen.n-1) 2.30; repmat(2.35,1,P.CorArtVen.n-1) 2.35; repmat(3,1,P.CorArtVen.n-1) 2.40];
%     P.CorArtVen.facDil(3,:) = min(3,facDilTarget(1:P.CorArtVen.n));
%     P.CorArtVen.facDil(4,:) = min(3,facDilTarget(P.CorArtVen.n+1:2*P.CorArtVen.n));
%     P.CorArtVen.facDil(5,:) = min(3,facDilTarget(2*P.CorArtVen.n+1:3*P.CorArtVen.n));
    P.CorArtVen.facDil(3,:) = min(facDilMax(1,:),facDilTarget(1:P.CorArtVen.n));
    P.CorArtVen.facDil(4,:) = min(facDilMax(2,:),facDilTarget(P.CorArtVen.n+1:2*P.CorArtVen.n));
    P.CorArtVen.facDil(5,:) = min(facDilMax(3,:),facDilTarget(2*P.CorArtVen.n+1:3*P.CorArtVen.n));
    facDil = P.CorArtVen.facDil;
%     P.CorArtVen.facDil(1,:) = 1+min(0.35,.35*((facDil(5,:)-1))/(3-1));
%     P.CorArtVen.facDil(2,:) = 1+min(0.35,.35*((facDil(5,:)-1))/(3-1));
    P.CorArtVen.facDil(1,:) = 1+min(0.40,.40*max((facDil(3:5,:)-1)./(facDilMax-1)));
    P.CorArtVen.facDil(2,:) = 1+min(0.40,.40*max((facDil(3:5,:)-1)./(facDilMax-1)));
    facDil = P.CorArtVen.facDil;
    % Dilation (facDil) change between beats
    P.CorArtVen.dfacDil = facDil0-facDil;
    % Reference facDil
    P.CorArtVen.facDil0 = P.CorArtVen.facDil;
end