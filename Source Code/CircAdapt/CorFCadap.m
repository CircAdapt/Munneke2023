function P = CorFCadap(P)
    
    %----% Estimation of TargetFlow in coronary ArtVen's for flow control
    iSy    = ~contains(P.ArtVen.Name,'pu','IgnoreCase',true); % {'Ca','Br','Ce','Fe'} systemic ArtVens
    Cor    = P.CorArtVen.Name;      % Coronary ArtVens
    q0AvSy = P.ArtVen.q0AV(:,iSy);  % Get reference flow per systemic ArtVen
    q0AvCor= P.CorArtVen.q0AV;      % Get reference flow per coronary ArtVen
    qCorM1P= P.CorArtVen.qMyo1epi;  % simulated epi flows in coronary ArtVen
    qCorM1M= P.CorArtVen.qMyo1mid;  % simulated mid flows in coronary ArtVen
    qCorM1N= P.CorArtVen.qMyo1endo; % simulated endo flows in coronary ArtVen
    
    % Coronary Autoregulation for VWall adaptation
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
    %      Demand is calculated by means of sarcomere stress - length area
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
    VO2        = k1*SSA2Cor/P.General.tCycle + k2cor;  % [mmol/s]], Oxygen consumption 
    P.CorArtVen.VO2 = VO2;
    
    if max(mqCorD(:)) < 0.1e-9
        %-% Update VWall
        VO2r   = sum(P.CorArtVen.VO2.*P.CorArtVen.Cor2Patch',2)'./P.Patch.VWall;
        VO2rn  = VO2r(3:end-1)./sum(VO2r(3:end-1));
        VO2r0  = sum(P.CorArtVen.VO20ref.*P.CorArtVen.Cor2Patch',2)'./P.Patch.VWallref;
        VO2rn0 = VO2r0(3:end-1)./sum(VO2r0(3:end-1));
        VO2d   = VO2rn-VO2rn0;
%         FacVO2ControlCor = exp(0.1*(1-VO2rn./VO2rn0));
        if isfield(P.CorArtVen,'VO2d')
%             if sum(abs(VO2d))>sum(abs(P.CorArtVen.VO2d))
%                 if isfield(P.CorArtVen,'alpha')
%                     alpha = P.CorArtVen.alpha;
%                 else
%                     alpha = VO2rn0(1)./VO2rn0(end);
%                 end
%                 alpha2 = 1.05*alpha;
%                 P.CorArtVen.alpha = alpha2;
%                 facalpha = ones(size(VO2rn0)); facalpha(1:12) = alpha2;
%                 VO2rn0f = ones(size(VO2rn0)).*facalpha./sum(ones(size(VO2rn0)).*facalpha).*sum(VO2rn0);
% %                 VO2rn0f2 = VO2rn0f./VO2rn0;
%                 FacVO2ControlCor = exp(0.1*(1-VO2rn./(VO2rn0f)));
%             else
                FacVO2ControlCor = exp(0.1*(1-VO2rn./(VO2rn0)));
%             end
% %             FacVO2ControlCor(abs(VO2d)>abs(P.CorArtVen.VO2d)) = 1;
% %             P.CorArtVen.VO2d(abs(VO2d)<abs(P.CorArtVen.VO2d)) = VO2d(abs(VO2d)<abs(P.CorArtVen.VO2d));
        else
            FacVO2ControlCor = exp(0.1*(1-VO2rn./(VO2rn0)));
        end
        P.CorArtVen.VO2d = VO2d;
        VWallTarget = VWall(3:end-1)./FacVO2ControlCor;
        VWallfmax = 0.25;
        VWallTarget = min((1+VWallfmax)*P.Patch.VWallref(3:end-1),max((1-VWallfmax)*P.Patch.VWallref(3:end-1),VWallTarget));
        P.Patch.VWall(3:end-1) = VWallTarget;
        
        %-% Update Demand
        q0AV       = P.CorArtVen.q0AV;
        VO20       = P.CorArtVen.VO20;
        q0AVt      = VO2.*q0AV./VO20;
        facq0AV    = q0AVt./q0AV;
        qSyTarget  = q0AvSy/sum(q0AvSy)*(P.General.q0-sum(facq0AV.*q0AV)); % Target flow in systemic and coronary ArtVen's
        P.ArtVen.q0AV(iSy)      = qSyTarget;
        P.CorArtVen.facq0AV     = facq0AV;
        P.Patch.SSAr0           = SSAr;      % Total energy per patch
        P.CorArtVen.SSA2Cor0    = SSA2Cor;   % Total energy per coronary
    else
        P.Patch.SSAr            = SSAr;      % Total energy per patch
        P.CorArtVen.SSA2Cor     = SSA2Cor;   % Total energy per coronary
    end
    
    %----% Third, Autoregulation by means of vasodilation/contriction (facDil)
    %      Occurs mainly in the arterioles
    kAvCor  = P.CorArtVen.kAV;
    facq0AV = P.CorArtVen.facq0AV;
    q0AvCor = q0AvCor.*facq0AV;
    facDil0 = P.CorArtVen.facDil;
    facDil  = [facDil0(3,:) facDil0(4,:) facDil0(5,:)];
    FacqControlCor = exp(0.3*(1-mean([qCorM1P, qCorM1M, qCorM1N])./([0.95*q0AvCor q0AvCor 1.05*q0AvCor]/3)));
    facDilTarget = facDil./FacqControlCor.^(-1./[kAvCor kAvCor kAvCor])/P.General.FacpControl;
    facDilMax      = [repmat(1.70,1,P.CorArtVen.n-1) 2.30; repmat(2.35,1,P.CorArtVen.n-1) 2.35; repmat(3,1,P.CorArtVen.n-1) 2.40];
%     P.CorArtVen.facDil(3,:) = min(3,facDilTarget(1:P.CorArtVen.n));
%     P.CorArtVen.facDil(4,:) = min(3,facDilTarget(P.CorArtVen.n+1:2*P.CorArtVen.n));
%     P.CorArtVen.facDil(5,:) = min(3,facDilTarget(2*P.CorArtVen.n+1:3*P.CorArtVen.n));
    P.CorArtVen.facDil(3,:) = min(facDilMax(1,:),facDilTarget(1:P.CorArtVen.n));
    P.CorArtVen.facDil(4,:) = min(facDilMax(2,:),facDilTarget(P.CorArtVen.n+1:2*P.CorArtVen.n));
    P.CorArtVen.facDil(5,:) = min(facDilMax(3,:),facDilTarget(2*P.CorArtVen.n+1:3*P.CorArtVen.n));
    facDil = P.CorArtVen.facDil;
%     P.CorArtVen.facDil(1,:) = 1+min(0.35,.35*max(([facDil(3,:); facDil(4,:); facDil(5,:)]-1))/(3-1));
%     P.CorArtVen.facDil(2,:) = 1+min(0.35,.35*max(([facDil(3,:); facDil(4,:); facDil(5,:)]-1))/(3-1));
    P.CorArtVen.facDil(1,:) = 1+min(0.40,.40*max((facDil(3:5,:)-1)./(facDilMax-1)));
    P.CorArtVen.facDil(2,:) = 1+min(0.40,.40*max((facDil(3:5,:)-1)./(facDilMax-1)));
    facDil = P.CorArtVen.facDil;
    % Dilation (facDil) change between beats
    P.CorArtVen.dfacDil = facDil0-facDil;
end