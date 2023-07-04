%% Reconstruct Figures

%% Figure 2 - External work and potential energy
% Gather data
% load('PRef.mat')
load('PRef_acuteLBBB.mat')
% load('PRef_chronicLBBB.mat')
[Ef, Sf, Sfa, Efr, Sfr, ED, SSA] = calcSSA(P);
[EDSSR_Ef,EDSSR_Sf] = calcEDSSR(P);

%----% Create figure
f=figure;
ititle = {{'Septum','(Early-activated)'},{'Left ventricular free wall','(Late-activated)'}};
ipatch = [19, 5]; % Septum, LV free wall
% Subplots
for is = 1:2
    subplot(1,2,is); hold on;

    % Active end-systolic stress strain relation
    plot(Efr(:,ipatch(is)),Sfr(:,ipatch(is))*1e-3,'LineWidth',3,'Color',[0 0 0]);

    % Passive end-diastolic stress strain relation
    plot(EDSSR_Ef,EDSSR_Sf(:,ipatch(is))*1e-3,'LineWidth',3,'Color',[0 0 0]);

    % Stress-strain loop
    plot(Ef(:,ipatch(is)),Sf(:,ipatch(is))*1e-3,'LineWidth',3,'Color',[0 0 0]);
        
    % End-systolic point
    plot(Efr(1,ipatch(is)),Sfr(1,ipatch(is))*1e-3,'p','Color',[ 47  85 151]/255,'MarkerFaceColor',[ 47  85 151]/255,'MarkerSize',16);

    % End-diastolic point
    plot(Ef(ED(ipatch(is)),ipatch(is)),P.Patch.Sf(ED(ipatch(is)),ipatch(is))*1e-3,'x','Color',[ 47  85 151]/255,'LineWidth',3,'MarkerSize',16);

    % Subplot design
    xlim([-0.3 0.2]); xticks(-0.4:0.1:0.2); xlabel('Strain (-)')
    ylim([0 60]); ylabel('Stress (kPa)');
    title(ititle{is});
    set(gca,'LineWidth',2,'FontSize',18);
end
% Figure design
set(gcf,'Color','white');
set(f,'Position',[0.2042    0.3420    1.3312    0.4200]*1e3);

% Second x axis in sarcomere length
f=figure;
subplot(1,2,1);
plot(Ef,Sf*1e-3);
xlim([-0.3 0.2]); xticks(-0.4:0.1:0.2); xlabel('Strain (-)');
ylim([0 80]); ylabel('Stress (kPa)');
set(gca,'LineWidth',2,'FontSize',18);
subplot(1,2,2);
semilogx(exp(Ef)*2,Sf*1e-3);
xlim([exp(-0.3)*2 exp(0.2)*2]); xlabel('Sarcomere Length (\mum)');
xticks(1.5:0.1:2.5); xticklabels({'','1.6','','1.8','','2.0','','2.2','','2.4'});
ylim([0 80]); ylabel('Stress (kPa)');
set(gca,'LineWidth',2,'FontSize',18);
% Figure design
set(gcf,'Color','white');
set(f,'Position',[0.2042    0.3420    1.3312    0.4200]*1e3);

%% Figure 3 - Resting arterial flow velocity
%----% Color scheme
cRCA = [  0 114 189]/255;
cLAD = [237 126  50]/255;
cLCx = [158   0   0]/255;

%----% Load data
% REFERENCE (synchronous)
load('PRef.mat');
% Valve opening/closing
[nVAvalveOp,nVAvalveCl,nAVvalveOp,nAVvalveCl] = ValveEvents(P,1,'L');
% Coronary territory mass
MCor  = 1.055e6*sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)'; % [g]
% Gather Arterial Flow Velocity Data
qVe   = P.CorArtVen.qVe*1e6;   
qArT  = Get('Tube','q',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^3 / s]
qArT  = qArT(end-length(qVe)+1:end,:);                % [m^3 / s]
qArT  = qArT*60e6;%./[sum(MCor(1:5)) sum(MCor(6:12)) sum(MCor(13:end))]; % [mL/min/g]
AArT  = Get('Tube','A',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^2]
vArT  = qArT./AArT;                                   % [m / s]

% ACUTE LBBB (asynchronous)
load('PRef_acuteLBBB.mat');
% Valve opening/closing
[nVAvalveOpa,nVAvalveCla,nAVvalveOpa,nAVvalveCla] = ValveEvents(P,1,'L');
% Coronary territory mass
MCor  = 1.055e6*sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)'; % [g]
% Gather Arterial Flow Velocity Data
qVe   = P.CorArtVen.qVe*1e6;   
qArTa = Get('Tube','q',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^3 / s]
qArTa = qArTa(end-length(qVe)+1:end,:);               % [m^3 / s]
qArTa = qArTa*60e6;%./[sum(MCor(1:5)) sum(MCor(6:12)) sum(MCor(13:end))]; % [mL/min/g]
AArT  = Get('Tube','A',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^2]
vArTa = qArTa./AArT;                                  % [m / s]

% CHRONIC LBBB (asynchronous)
load('PRef_chronicLBBB.mat'); 
% Coronary territory mass
MCor  = 1.055e6*sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)'; % [g]
% Gather Arterial Flow Velocity Data
qVe   = P.CorArtVen.qVe*1e6;   
qArTav= Get('Tube','q',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^3 / s]
qArTav= qArTav(end-length(qVe)+1:end,:);              % [m^3 / s]
qArTav= qArTav*60e6;%./[sum(MCor(1:5)) sum(MCor(6:12)) sum(MCor(13:end))]; % [mL/min/g]
AArT  = Get('Tube','A',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^2]
vArTav= qArTav./AArT;                                 % [m / s]

%----% Information per Coronary
FigNum  = [311 312 313];            % Figure number
cCor    = [cLCx; cLAD; cRCA];       % Coronary color
% yylim   = [-0.25 1.50];             % YLimit
yylim   = [-100  600];               % YLimit

%----% Resting Arterial Flow Velocity Figures
for iAr = 1:size(qArT,2)
    % Figure
    figure(FigNum(iAr)); 
    s1=subplot(1,2,1); hold on;
    % Ejection area in grey
    area([nVAvalveOp-nAVvalveCl nVAvalveCl-nAVvalveCl]*P.General.tCycle/length(qVe),[yylim(2) yylim(2)],...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.10,'BaseValue',yylim(1))
    % Valve opening/closing lines
%     plot(repmat(nAVvalveCl-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nVAvalveCl-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nAVvalveOp-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
%     plot(repmat(nVAvalveOp-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
    % Horizontal line at zero
    plot([0 P.General.tCycle],[0 0],'k:','LineWidth',1)
    % Coronary Artery Flow Velocity
%     plot([P.t(1:end)-P.t(1)],[vArT(nAVvalveCl:end,iAr) ;vArT(1:nAVvalveCl-1,iAr)],'-','LineWidth',3,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[qArT(nAVvalveCl:end,iAr) ;qArT(1:nAVvalveCl-1,iAr)],'-','LineWidth',3,'Color',cCor(iAr,:));
    % Subplot Design
    xlim([0 P.General.tCycle]);
    ylim(yylim);
    set(gca,'XTick',[(nVAvalveCl-nAVvalveCl)*P.General.tCycle/length(qVe) P.General.tCycle ],'XTickLabel',{'ES','ED'});
    set(gca,'XTick',[],'XTickLabel',{},'XColor','none');
%     set(gca,'YTick',-0.25:0.125:1.50,'YTickLabel',{'-0.25','','0','','0.25','','0.50','','0.75','','1.00','','1.25','','1.50'});
%     set(gca,'YTick',-150:75:750,'YTickLabel',{'-150','','0','','150','','300','','450','','600','','750'});
    set(gca,'YTick',yylim(1):50:yylim(2),'YTickLabel',{'-100','','0','','100','','200','','300','','400','','500','','600'});
    set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','box','off','YGrid','on','GridAlpha',0.05);
%     xlabel('Time (s)');
%     ylabel('Flow velocity (m/s)');
    ylabel('Flow (ml/min)');
    set(s1,'Position',[0.1575    0.1267    0.3484    0.7983]);
    
    s2=subplot(1,2,2); hold on;
    % Ejection area in grey
    area([nVAvalveOpa-nAVvalveCla nVAvalveCla-nAVvalveCla]*P.General.tCycle/length(qVe),[yylim(2) yylim(2)],...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.10,'BaseValue',yylim(1))
    % Valve opening/closing lines
%     plot(repmat(nAVvalveCla-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nVAvalveCla-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nAVvalveOpa-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
%     plot(repmat(nVAvalveOpa-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
    % Horizontal line at zero
    plot([0 P.General.tCycle],[0 0],'k:','LineWidth',1)
    % Coronary Artery Flow Velocity
%     plot([P.t(1:end)-P.t(1)],[vArTa(nAVvalveCla:end,iAr);vArTa(1:nAVvalveCla-1,iAr)],'-','LineWidth',3,'Color',cCor(iAr,:));
%     plot([P.t(1:end)-P.t(1)],[vArTav(nAVvalveCla:end,iAr);vArTav(1:nAVvalveCla-1,iAr)],'-','LineWidth',1.5,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[qArTa(nAVvalveCla:end,iAr);qArTa(1:nAVvalveCla-1,iAr)],'-','LineWidth',3,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[qArTav(nAVvalveCla:end,iAr);qArTav(1:nAVvalveCla-1,iAr)],'-','LineWidth',1.5,'Color',cCor(iAr,:));
    % Subplot Design
    xlim([0 P.General.tCycle]);
    ylim(yylim);
    set(gca,'XTick',[(nVAvalveCla-nAVvalveCla)*P.General.tCycle/length(qVe) P.General.tCycle],'XTickLabel',{'ES','ED'});
    set(gca,'XTick',[],'XTickLabel',{},'XColor','none');
%     set(gca,'YTick',-0.25:0.125:1.50,'YTickLabel',{''});
%     set(gca,'YTick',-150:75:750,'YTickLabel',{''});
    set(gca,'YTick',yylim(1):50:yylim(2),'YTickLabel',{''},'YColor','none');
    set(gcf,'Color','white')
    set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','box','off','YGrid','on','GridAlpha',0.05);
    set(s2,'Position',[0.5478    0.1267    0.3484    0.7983]);
    
    % Figure Design
    set(figure(FigNum(iAr)),'Position',[588.2000  169.8000  636.8000  574.4000]);
    set(gcf,'Color','white')
end

%% Figure 3 - Hyperemic arterial flow velocity
%----% Color scheme
cRCA = [  0 114 189]/255;
cLAD = [237 126  50]/255;
cLCx = [158   0   0]/255;

%----% Load data
% REFERENCE (synchronous)
load('PRef_Hyp.mat');
% Valve opening/closing
[nVAvalveOp,nVAvalveCl,nAVvalveOp,nAVvalveCl] = ValveEvents(P,1,'L');
% Coronary territory mass
MCor  = 1.055e6*sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)'; % [g]
% Gather Arterial Flow Velocity Data
qVe   = P.CorArtVen.qVe*1e6;   
qArT  = Get('Tube','q',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^3 / s]
qArT  = qArT(end-length(qVe)+1:end,:);                % [m^3 / s]
qArT  = qArT*60e6;%./[sum(MCor(1:5)) sum(MCor(6:12)) sum(MCor(13:end))]; % [mL/min/g]
AArT  = Get('Tube','A',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^2]
vArT  = qArT./AArT;                                   % [m / s]

% ACUTE LBBB (asynchronous)
load('PRef_acuteLBBB_Hyp.mat');
% Valve opening/closing
[nVAvalveOpa,nVAvalveCla,nAVvalveOpa,nAVvalveCla] = ValveEvents(P,1,'L');
% Coronary territory mass
MCor  = 1.055e6*sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)'; % [g]
% Gather Arterial Flow Velocity Data
qVe   = P.CorArtVen.qVe*1e6;   
qArTa = Get('Tube','q',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^3 / s]
qArTa = qArTa(end-length(qVe)+1:end,:);               % [m^3 / s]
qArTa = qArTa*60e6;%./[sum(MCor(1:5)) sum(MCor(6:12)) sum(MCor(13:end))]; % [mL/min/g]
AArT  = Get('Tube','A',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^2]
vArTa = qArTa./AArT;                                  % [m / s]

% CHRONIC LBBB (asynchronous)
load('PRef_chronicLBBB_Hyp.mat');
% Coronary territory mass
MCor  = 1.055e6*sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)'; % [g]
% Gather Arterial Flow Velocity Data
qVe   = P.CorArtVen.qVe*1e6;   
qArTav= Get('Tube','q',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^3 / s]
qArTav= qArTav(end-length(qVe)+1:end,:);              % [m^3 / s]
qArTav= qArTav*60e6;%./[sum(MCor(1:5)) sum(MCor(6:12)) sum(MCor(13:end))]; % [mL/min/g]
AArT  = Get('Tube','A',{'LMLCx1','LMLAD1','AoRCA1'}); % [m^2]
vArTav= qArTav./AArT;                                 % [m / s]

%----% Information per Coronary
FigNum  = [321 322 323];            % Figure number
cCor    = [cLCx; cLAD; cRCA];       % Coronary color
yylim   = [-100  600];               % YLimit

%----% Hyperemic Arterial Flow Velocity Figures
for iAr = 1:size(qArT,2)
    % Figure
    figure(FigNum(iAr));  
    s1=subplot(1,2,1); hold on;
    % Ejection area in grey
    area([nVAvalveOp-nAVvalveCl nVAvalveCl-nAVvalveCl]*P.General.tCycle/length(qVe),[yylim(2) yylim(2)],...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.10,'BaseValue',yylim(1))
    % Valve opening/closing lines
%     plot(repmat(nAVvalveCl-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nVAvalveCl-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nAVvalveOp-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
%     plot(repmat(nVAvalveOp-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
    % Horizontal line at zero
    plot([0 P.General.tCycle],[0 0],'k:','LineWidth',1)
    % Coronary Artery Flow Velocity
%     plot([P.t(1:end)-P.t(1)],[vArT(nAVvalveCl:end,iAr) ;vArT(1:nAVvalveCl-1,iAr)],'-','LineWidth',3,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[qArT(nAVvalveCl:end,iAr) ;qArT(1:nAVvalveCl-1,iAr)],'-','LineWidth',3,'Color',cCor(iAr,:));
    % Figure design
    xlim([0 P.General.tCycle]);
    ylim(yylim);
    set(gca,'XTick',[(nVAvalveCl-nAVvalveCl)*P.General.tCycle/length(qVe) P.General.tCycle],'XTickLabel',{'ES','ED'});
    set(gca,'XTick',[],'XTickLabel',{},'XColor','none');
%     set(gca,'YTick',-0.25:0.125:1.50,'YTickLabel',{'-0.25','','0','','0.25','','0.50','','0.75','','1.00','','1.25','','1.50'});
%     set(gca,'YTick',-150:75:750,'YTickLabel',{'-150','','0','','150','','300','','450','','600','','750'});
    set(gca,'YTick',yylim(1):50:yylim(2),'YTickLabel',{'-100','','0','','100','','200','','300','','400','','500','','600'});
    set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','box','off','YGrid','on','GridAlpha',0.05);
    xlabel('Time (s)');
%     ylabel('Flow velocity (m/s)');
    ylabel('Flow (ml/min)');
    set(s1,'Position',[0.1575    0.1267    0.3484    0.7983]);
    
    s2=subplot(1,2,2); hold on;
    % Ejection area in grey
    area([nVAvalveOpa-nAVvalveCla nVAvalveCla-nAVvalveCla]*P.General.tCycle/length(qVe),[yylim(2) yylim(2)],...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.10,'BaseValue',yylim(1))
    % Valve opening/closing lines
%     plot(repmat(nAVvalveCla-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nVAvalveCla-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nAVvalveOpa-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
%     plot(repmat(nVAvalveOpa-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
    % Horizontal line at zero
    plot([0 P.General.tCycle],[0 0],'k:','LineWidth',1)
    % Coronary Artery Flow Velocity
%     plot([P.t(1:end)-P.t(1)],[vArTa(nAVvalveCla:end,iAr) ;vArTa(1:nAVvalveCla-1,iAr) ],'-','LineWidth',3,'Color',cCor(iAr,:));
%     plot([P.t(1:end)-P.t(1)],[vArTav(nAVvalveCla:end,iAr);vArTav(1:nAVvalveCla-1,iAr)],'-','LineWidth',1.5,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[qArTa(nAVvalveCla:end,iAr) ;qArTa(1:nAVvalveCla-1,iAr) ],'-','LineWidth',3,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[qArTav(nAVvalveCla:end,iAr);qArTav(1:nAVvalveCla-1,iAr)],'-','LineWidth',1.5,'Color',cCor(iAr,:));
    % Subplot design
    xlim([0 P.General.tCycle]);
    ylim(yylim);
    set(gca,'XTick',[(nVAvalveCla-nAVvalveCla)*P.General.tCycle/length(qVe) P.General.tCycle],'XTickLabel',{'ES','ED'});
    set(gca,'XTick',[],'XTickLabel',{},'XColor','none');
%     set(gca,'YTick',-0.25:0.125:1.50,'YTickLabel',{''});
%     set(gca,'YTick',-150:75:750,'YTickLabel',{''});
    set(gca,'YTick',yylim(1):50:yylim(2),'YTickLabel',{''},'YColor','none');
    set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','box','off','YGrid','on','GridAlpha',0.05);
    set(s2,'Position',[0.5478    0.1267    0.3484    0.7983]);
    
    % Figure Design
    set(figure(FigNum(iAr)),'Position',[588.2000  169.8000  636.8000  574.4000]);
    set(gcf,'Color','white')
end

%% Figure 3 - Average intramyocardial pressure
%----% Color scheme
cRCA = [  0 114 189]/255;
cLAD = [237 126  50]/255;
cLCx = [158   0   0]/255;

%----% Load data
% REFERENCE (synchronous)
load('PRef.mat');
% Valve opening/closing
[nVAvalveOp,nVAvalveCl,nAVvalveOp,nAVvalveCl] = ValveEvents(P,1,'L');
% Gather Intramyocardial Pressure Data
qVe   = P.CorArtVen.qVe*1e6;   
Pw    = sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)';   % [m3]
CEPn  = P.CorArtVen.pimCEP(:,2:3:end)*7.5e-3;           % [mmHg]
VE    = P.CorArtVen.pimVE*7.5e-3;                       % [mmHg]  
CEPc  = [sum((P.CorArtVen.pimCEP(:, 2:3:15 )*7.5e-3.*Pw( 1:5  )),2)/sum(Pw( 1:5  )),...
         sum((P.CorArtVen.pimCEP(:,17:3:35 )*7.5e-3.*Pw( 6:12 )),2)/sum(Pw( 6:12 )),...
         sum((P.CorArtVen.pimCEP(:,38:3:end)*7.5e-3.*Pw(13:end)),2)/sum(Pw(13:end))]; % [mmHg]
VEc   = [sum((P.CorArtVen.pimVE( :, 1:5    )*7.5e-3.*Pw( 1:5  )),2)/sum(Pw( 1:5  )),...
         sum((P.CorArtVen.pimVE( :, 6:12   )*7.5e-3.*Pw( 6:12 )),2)/sum(Pw( 6:12 )),...
         sum((P.CorArtVen.pimVE( :,13:end  )*7.5e-3.*Pw(13:end)),2)/sum(Pw(13:end))]; % [mmHg]
VEmax = [max(P.CorArtVen.pimVE(:, 1:5  )'*7.5e-3)',...
         max(P.CorArtVen.pimVE(:, 6:12 )'*7.5e-3)',...
         max(P.CorArtVen.pimVE(:,13:end)'*7.5e-3)'];
VEmin = [min(P.CorArtVen.pimVE(:, 1:5  )'*7.5e-3)',...
         min(P.CorArtVen.pimVE(:, 6:12 )'*7.5e-3)',...
         min(P.CorArtVen.pimVE(:,13:end)'*7.5e-3)'];
     
% ACUTE LBBB (asynchronous)
load('PRef_acuteLBBB.mat');
% Valve opening/closing
[nVAvalveOpa,nVAvalveCla,nAVvalveOpa,nAVvalveCla] = ValveEvents(P,1,'L');
% Gather Intramyocardial Pressure Data
qVe   = P.CorArtVen.qVe*1e6;   
Pw    = sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)';   % [m3]
CEPn  = P.CorArtVen.pimCEP(:,3:3:end)*7.5e-3;           % [mmHg]
VE    = P.CorArtVen.pimVE*7.5e-3;                       % [mmHg]  
CEPca = [sum((P.CorArtVen.pimCEP(:, 2:3:15 )*7.5e-3.*Pw( 1:5  )),2)/sum(Pw( 1:5  )),...
         sum((P.CorArtVen.pimCEP(:,17:3:35 )*7.5e-3.*Pw( 6:12 )),2)/sum(Pw( 6:12 )),...
         sum((P.CorArtVen.pimCEP(:,38:3:end)*7.5e-3.*Pw(13:end)),2)/sum(Pw(13:end))]; % [mmHg]
VEca  = [sum((P.CorArtVen.pimVE( :, 1:5    )*7.5e-3.*Pw( 1:5  )),2)/sum(Pw( 1:5  )),...
         sum((P.CorArtVen.pimVE( :, 6:12   )*7.5e-3.*Pw( 6:12 )),2)/sum(Pw( 6:12 )),...
         sum((P.CorArtVen.pimVE( :,13:end  )*7.5e-3.*Pw(13:end)),2)/sum(Pw(13:end))]; % [mmHg]
VEmaxa= [max(P.CorArtVen.pimVE(:, 1:5  )'*7.5e-3)',...
         max(P.CorArtVen.pimVE(:, 6:12 )'*7.5e-3)',...
         max(P.CorArtVen.pimVE(:,13:end)'*7.5e-3)'];
VEmina= [min(P.CorArtVen.pimVE(:, 1:5  )'*7.5e-3)',...
         min(P.CorArtVen.pimVE(:, 6:12 )'*7.5e-3)',...
         min(P.CorArtVen.pimVE(:,13:end)'*7.5e-3)'];
     
% CHRONIC LBBB (asynchronous)   
load('PRef_chronicLBBB.mat');
% Gather Intramyocardial Pressure Data
qVe   = P.CorArtVen.qVe*1e6;   
Pw    = sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)';   % [m3]
CEPn  = P.CorArtVen.pimCEP(:,3:3:end)*7.5e-3;           % [mmHg]
VE    = P.CorArtVen.pimVE*7.5e-3;                       % [mmHg]  
CEPcav= [sum((P.CorArtVen.pimCEP(:, 2:3:15 )*7.5e-3.*Pw( 1:5  )),2)/sum(Pw( 1:5  )),...
         sum((P.CorArtVen.pimCEP(:,17:3:35 )*7.5e-3.*Pw( 6:12 )),2)/sum(Pw( 6:12 )),...
         sum((P.CorArtVen.pimCEP(:,38:3:end)*7.5e-3.*Pw(13:end)),2)/sum(Pw(13:end))]; % [mmHg]
VEcav = [sum((P.CorArtVen.pimVE( :, 1:5    )*7.5e-3.*Pw( 1:5  )),2)/sum(Pw( 1:5  )),...
         sum((P.CorArtVen.pimVE( :, 6:12   )*7.5e-3.*Pw( 6:12 )),2)/sum(Pw( 6:12 )),...
         sum((P.CorArtVen.pimVE( :,13:end  )*7.5e-3.*Pw(13:end)),2)/sum(Pw(13:end))]; % [mmHg]
VEmaxav=[max(P.CorArtVen.pimVE(:, 1:5  )'*7.5e-3)',...
         max(P.CorArtVen.pimVE(:, 6:12 )'*7.5e-3)',...
         max(P.CorArtVen.pimVE(:,13:end)'*7.5e-3)'];
VEminav=[min(P.CorArtVen.pimVE(:, 1:5  )'*7.5e-3)',...
         min(P.CorArtVen.pimVE(:, 6:12 )'*7.5e-3)',...
         min(P.CorArtVen.pimVE(:,13:end)'*7.5e-3)'];

%----% Information per Coronary
FigNum  = [321 322 323];            % Figure number
cCor    = [cLCx; cLAD; cRCA];       % Coronary color
yylim   = [-20 120];                % YLimit

%----% Average Intramyocardial Pressure Figures
for iAr = 1:size(qArT,2)
    % Figure
    figure(FigNum(iAr)); 
    s1=subplot(1,2,1);
    yyaxis left; set(gca,'YColor','none'); ylim(yylim); set(gca,'YTick',-20:10:120,'YTickLabel',{''});
    yyaxis right; set(gca,'YColor',[0 0 0]);
    hold on;
    % Ejection area in grey
    area([nVAvalveOp-nAVvalveCl nVAvalveCl-nAVvalveCl]*P.General.tCycle/length(qVe),[yylim(2) yylim(2)],...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.10,'BaseValue',yylim(1))
    % Valve opening/closing lines
%     plot(repmat(nAVvalveCl-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nVAvalveCl-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nAVvalveOp-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
%     plot(repmat(nVAvalveOp-nAVvalveCl,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
    % Horizontal line at zero
    plot([0 P.General.tCycle],[0 0],'k:','LineWidth',1)
    % Intramyocardial pressure - CEP
    plot([P.t(1:end)-P.t(1)],[CEPc(nAVvalveCl:end,iAr) ;CEPc(1:nAVvalveCl-1,iAr) ],':','LineWidth',3,'Color',cCor(iAr,:));
    % Intramyocardial pressure - VE
    plot([P.t(1:end)-P.t(1)],[VEc(nAVvalveCl:end,iAr) ;VEc(1:nAVvalveCl-1,iAr) ],'-.','LineWidth',3,'Color',cCor(iAr,:));
    % Intramyocardial pressure - average
    plot([P.t(1:end)-P.t(1)],[CEPc(nAVvalveCl:end,iAr)+VEc(nAVvalveCl:end,iAr) ;CEPc(1:nAVvalveCl-1,iAr)+VEc(1:nAVvalveCl-1,iAr)],'-','LineWidth',3,'Color',cCor(iAr,:));
    % Subplot design
    xlim([0 P.General.tCycle]);
    ylim(yylim);
    set(gca,'XTick',[(nVAvalveCl-nAVvalveCl)*P.General.tCycle/length(qVe) P.General.tCycle],'XTickLabel',{'ES','ED'});
    set(gca,'XTick',[],'XTickLabel',{},'XColor','none');
    set(gca,'YTick',-20:10:120,'YTickLabel',{''},'YColor','none');
    set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','box','off','YGrid','on','GridAlpha',0.05);
    set(s1,'Position',[0.0575    0.1267    0.3484    0.7983]);
    
    s2=subplot(1,2,2);
    yyaxis left; set(gca,'YColor','none'); ylim(yylim); set(gca,'YTick',-20:10:120,'YTickLabel',{''});
    yyaxis right; set(gca,'YColor',[0 0 0]);
    hold on;
    % Ejection area in grey
    area([nVAvalveOpa-nAVvalveCla nVAvalveCla-nAVvalveCla]*P.General.tCycle/length(qVe),[yylim(2) yylim(2)],...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.10,'BaseValue',yylim(1))
    % Valve opening/closing lines
%     plot(repmat(nAVvalveCla-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nVAvalveCla-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k-','LineWidth',1)
%     plot(repmat(nAVvalveOpa-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
%     plot(repmat(nVAvalveOpa-nAVvalveCla,1,2)*P.General.tCycle/length(qVe),yylim,'k--','LineWidth',1)
    % Horizontal line at zero
    plot([0 P.General.tCycle],[0 0],'k:','LineWidth',1)
    % Intramyocardial pressure - CEP
    plot([P.t(1:end)-P.t(1)],[CEPca(nAVvalveCla:end,iAr)  ;CEPca(1:nAVvalveCla-1,iAr) ],':','LineWidth',3  ,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[CEPcav(nAVvalveCla:end,iAr) ;CEPcav(1:nAVvalveCla-1,iAr)],':','LineWidth',1.5,'Color',cCor(iAr,:));
    % Intramyocardial pressure - VE
    plot([P.t(1:end)-P.t(1)],[VEca(nAVvalveCla:end,iAr)   ;VEca(1:nAVvalveCla-1,iAr)  ],'-.','LineWidth',3  ,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[VEcav(nAVvalveCla:end,iAr)  ;VEcav(1:nAVvalveCla-1,iAr) ],'-.','LineWidth',1.5,'Color',cCor(iAr,:));
    % Intramyocardial pressure - average
    plot([P.t(1:end)-P.t(1)],[CEPca(nAVvalveCla:end,iAr)+VEca(nAVvalveCla:end,iAr)  ;CEPca(1:nAVvalveCla-1,iAr)+VEca(1:nAVvalveCla-1,iAr) ],'-','LineWidth',3  ,'Color',cCor(iAr,:));
    plot([P.t(1:end)-P.t(1)],[CEPcav(nAVvalveCla:end,iAr)+VEcav(nAVvalveCla:end,iAr) ;CEPcav(1:nAVvalveCla-1,iAr)+VEcav(1:nAVvalveCla-1,iAr)],'-','LineWidth',1.5,'Color',cCor(iAr,:));
    % Subplot design
    xlim([0 P.General.tCycle]);
    ylim(yylim);
    set(gca,'XTick',[(nVAvalveCla-nAVvalveCla)*P.General.tCycle/length(qVe) P.General.tCycle],'XTickLabel',{'ES','ED'});
    set(gca,'XTick',[],'XTickLabel',{},'XColor','none');
    set(gca,'YTick',-20:10:120,'YTickLabel',{'-20','','0','','20','','40','','60','','80','','100','','120'});
    set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','box','off','YGrid','on','GridAlpha',0.05);
    xlabel('Time (s)');
    ylabel('Intramyocardial pressure (mmHg)');
    set(s2,'Position',[0.4478    0.1267    0.3484    0.7983]);
    
    % Figure Design
    set(gcf,'Color','white')
    set(figure(FigNum(iAr)),'Position',[588.2000  169.8000  636.8000  574.4000]);
end

%% Figure 4 - Activation delay
%----% Add subfolder AHABullseye to the path
addpath('AHABullseye');

%----% Load P-struct
PName = {'PRef.mat';...
         'PRef_acuteLBBB.mat';...
         'PRef_chronicLBBB.mat'};
iTitle = {{'REFERENCE'   ,'Activation Delay'},...
          {'ACUTE LBBB'  ,'Activation Delay'},...
          {'CHRONIC LBBB','Activation Delay'}};
for iPName = 1:3
    load(PName{iPName});

    %----% Create figure
    f=figure;
    Lab = 0; % 0 if labels off, 1 if labels on

    %----% Load RV data
    % Bullseye data
    [theta,rho] = getWedgeBorder(360, 2*360, 0, 2);
    wedgeHandle = polar(gca,theta,rho);
    XData = get(wedgeHandle,'XData');
    YData = get(wedgeHandle,'YData');
    XData=XData-1;

    % Color data
    cRV = Get('Patch','dT','Rv1');

    % Plot
    close(f);
    figure; hold on; axis equal;
    fillBullseyeRV(cRV,0,2,0,360);
    plot(XData,YData,'k');

    %----% Load LV data
    % Bullseye data
    c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);

    % Color data
    cLVac = Get('Patch','dT',{'Lv12'})*1e3+0.01;
    cLVa  = Get('Patch','dT',{'Lv11','Lv9','Sv5','Lv10'})*1e3;
    cLVm  = Get('Patch','dT',{'Lv8','Lv5','Sv3','Sv4','Lv6','Lv7'})*1e3;
    cLVb  = Get('Patch','dT',{'Lv4','Lv1','Sv1','Sv2','Lv2','Lv3'})*1e3;
    fillBullseye(0,0.0,0.5,0,360);          % to set color range
    fillBullseye(70,0.0,0.5,0,360);         % to set color range
    fillBullseye(cLVac,0.0,0.5,0,360);      % apical cap
    fillBullseye(cLVa ,0.5,1.0,-45,315);    % apex
    fillBullseye(cLVm ,1.0,1.5,0,360);      % mid
    fillBullseye(cLVb ,1.5,2.0,0,360);      % base
    uistack(c,'top');

    %----% Labels
    if Lab == 1
        dec = 0; % Amount of decimals
        % RV 
        ang = (0)*(pi/180); [x,y] = pol2cart(ang, -2.475);
        t = text(x,y,string(round(cRV,dec))); set(t,'HorizontalAlignment','Center');
        % LV Apical cap
        ang = (60*0)*(pi/180); [x,y] = pol2cart(ang, 0);
        t = text(x,y,string(round(cLVac,dec))); set(t,'HorizontalAlignment','Center');
        % LV Apex
        labels = string(round(cLVa,dec))';
        for k = 1:length(labels)
            ang = (90*(k-1))*(pi/180); [x,y] = pol2cart(ang, 0.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Mid
        labels = string(round(cLVm,dec))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.25);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Base
        labels = string(round(cLVb,dec));
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
    end

    %----% Define colormap and add colorbar
    mymap = flip([linspace(  1, 29, 250)'/255,linspace(  1, 62, 250)'/255,linspace(  1,129, 250)'/255;...
             linspace( 29,  9, 250)'/255,linspace( 62,141, 250)'/255,linspace(129,138, 250)'/255;...
             linspace(  9, 86, 250)'/255,linspace(141,196, 250)'/255,linspace(138, 60, 250)'/255;...
             linspace( 86,170, 250)'/255,linspace(196,213, 250)'/255,linspace( 60, 33, 250)'/255;...
             linspace(170,253, 250)'/255,linspace(213,231, 250)'/255,linspace( 33, 37, 250)'/255;...
             linspace(253,255, 250)'/255,linspace(231,255, 250)'/255,linspace( 37,255, 250)'/255]);
    colormap(mymap);
    b = colorbar;

    %----% Figure design
    set(b,'Ticks',[0 25:10:65],'TickLabels',{'0','25','35','45','55','65 ms'},'Location','east');
    set(b, 'Position',[0.7664    0.24    0.0381    0.7644-2*(0.24-0.1359)]);
    title(iTitle{iPName});
    set(gcf,'Color','white');
    ax=gca; ax.XColor = 'none'; ax.YColor = 'none';
end

%% Figure 4 - Wall thickness
%----% Add subfolder AHABullseye to the path
addpath('AHABullseye');

%----% Load P-struct
PName = {'PRef.mat';...
         'PRef_acuteLBBB.mat';...
         'PRef_chronicLBBB.mat'};
iTitle = {{'REFERENCE'   ,'Wall Thickness'},...
          {'ACUTE LBBB'  ,'Wall Thickness'},...
          {'CHRONIC LBBB','Wall Thickness'}};
for iPName = 1:3
    load(PName{iPName});

    %----% Create figure
    f=figure;
    Lab = 0; % 0 if labels off, 1 if labels on

    %----% RV
    % Bullseye data
    [theta,rho] = getWedgeBorder(360, 2*360, 0, 2);
    wedgeHandle = polar(gca,theta,rho);
    XData = get(wedgeHandle,'XData');
    YData = get(wedgeHandle,'YData');
    XData=XData-1;

    % Color data
    cRV = 1e3*(Get('Patch','VWall','Rv1'))./mean((Get('Patch','Am',{'Rv1'})));

    % Plot
    close(f);
    figure; hold on; axis equal;
    title('Wall thickness')
    fillBullseyeRV(cRV,0,2,0,360);
    plot(XData,YData,'k');

    %----% LV
    % Bullseye data
    c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);

    % Color data
    ED    = round(P.Patch.Depolarization(1,end)/P.General.Dt);
    Amac  = Get('Patch','Am',{'Lv12'});
    Ama   = Get('Patch','Am',{'Lv11','Lv9','Sv5','Lv10'});
    Amm   = Get('Patch','Am',{'Lv8','Lv5','Sv3','Sv4','Lv6','Lv7'});
    Amb   = Get('Patch','Am',{'Lv4','Lv1','Sv1','Sv2','Lv2','Lv3'});
    cLVac = 1e3*(Get('Patch','VWall',{'Lv12'}))./Amac(ED,:);
    cLVa  = 1e3*(Get('Patch','VWall',{'Lv11','Lv9','Sv5','Lv10'}))./Ama(ED,:);
    cLVm  = 1e3*(Get('Patch','VWall',{'Lv8','Lv5','Sv3','Sv4','Lv6','Lv7'}))./Amm(ED,:);
    cLVb  = 1e3*(Get('Patch','VWall',{'Lv4','Lv1','Sv1','Sv2','Lv2','Lv3'}))./Amb(ED,:);
    fillBullseye(4,0.0,0.5,0,360);          % to set color range
    fillBullseye(16.5,0.0,0.5,0,360);         % to set color range
    fillBullseye(cLVac,0.0,0.5,0,360);      % apical cap
    fillBullseye(cLVa ,0.5,1.0,-45,315);    % apex
    fillBullseye(cLVm ,1.0,1.5,0,360);      % mid
    fillBullseye(cLVb ,1.5,2.0,0,360);      % base
    uistack(c,'top');

    %----% Labels
    if Lab == 1
        dec = 1; % Amount of decimals
        % RV 
        ang = (0)*(pi/180); [x,y] = pol2cart(ang, -2.475);
        t = text(x,y,string(num2str(cRV','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apical cap
        ang = (60*0)*(pi/180); [x,y] = pol2cart(ang, 0);
        t = text(x,y,string(num2str(cLVac','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apex
        labels = string(num2str(cLVa','%.1f'))';
        for k = 1:length(labels)
            ang = (90*(k-1))*(pi/180); [x,y] = pol2cart(ang, 0.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Mid
        labels = string(num2str(cLVm','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.25);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Base
        labels = string(num2str(cLVb','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
    end
    
    %----% Define colormap and add colorbar
    mymap = [linspace(  1, 29, 250)'/255,linspace(  1, 62, 250)'/255,linspace(  1,129, 250)'/255;...
             linspace( 29,  9, 250)'/255,linspace( 62,141, 250)'/255,linspace(129,138, 250)'/255;...
             linspace(  9, 86, 250)'/255,linspace(141,196, 250)'/255,linspace(138, 60, 250)'/255;...
             linspace( 86,170, 250)'/255,linspace(196,213, 250)'/255,linspace( 60, 33, 250)'/255;...
             linspace(170,253, 250)'/255,linspace(213,231, 250)'/255,linspace( 33, 37, 250)'/255;...
             linspace(253,255, 250)'/255,linspace(231,255, 250)'/255,linspace( 37,255, 250)'/255];
    colormap(mymap);
    b = colorbar;
    
    %----% Figure design
    set(b,'Ticks',[4:1.25:16.5],'TickLabels',{'4.0','','6.5','','9.0','','11.5','','14.0','','16.5 mm'},'Location','east');
    set(b, 'Position',[0.7664    0.24    0.0381    0.7644-2*(0.24-0.1359)]);
    title(iTitle{iPName});
    set(gcf,'Color','white');
    ax=gca; ax.XColor = 'none'; ax.YColor = 'none';
end

%% Figure 4 - Hyperemic flow per tissue weight
%----% Add subfolder AHABullseye to the path
addpath('AHABullseye');

%----% Load P-struct
PName = {'PRef_Hyp.mat';...
         'PRef_acuteLBBB_Hyp.mat';...
         'PRef_chronicLBBB_Hyp.mat'};
iTitle = {{'REFERENCE'   ,'Hyperemic Flow'},...
          {'ACUTE LBBB'  ,'Hyperemic Flow'},...
          {'CHRONIC LBBB','Hyperemic Flow'}};
for iPName = 1:3
    load(PName{iPName});

    %----% Create figure
    f=figure;
    Lab = 0;  % 0 if labels off, 1 if labels on

    %----% RV
    % Bullseye data
    [theta,rho] = getWedgeBorder(360, 2*360, 0, 2);
    wedgeHandle = polar(gca,theta,rho);
    XData = get(wedgeHandle,'XData');
    YData = get(wedgeHandle,'YData');
    XData=XData-1;

    % Color data
    cRV = 60e6*mean(Get('CorArtVen','q','RCA'))./(1.055e6*Get('Patch','VWall','Rv1')); % [ml/min/g]

    % Plot
    close(f);
    figure; hold on; axis equal;
    title('Demand')
    fillBullseyeRV(cRV,0,2,0,360);
    plot(XData,YData,'k');

    %----% LV
    % Bullseye data
    c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);

    % Color data
    cLVac = (60e6*mean(Get('CorArtVen','q',P.CorArtVen.Name(12)))./(1.055e6*Get('Patch','VWall',{'Lv12'}))); % [ml/min/g]
    cLVa  = (60e6*mean(Get('CorArtVen','q',P.CorArtVen.Name([3,10,11,17])))./(1.055e6*Get('Patch','VWall',{'Lv11','Lv9','Sv5','Lv10'}))); % [ml/min/g]
    cLVm  = (60e6*mean(Get('CorArtVen','q',P.CorArtVen.Name([2,8,9,14,16,4])))./(1.055e6*Get('Patch','VWall',{'Lv8','Lv5','Sv3','Sv4','Lv6','Lv7'}))); % [ml/min/g]
    cLVb  = (60e6*mean(Get('CorArtVen','q',P.CorArtVen.Name([1,6,7,13,15,5])))./(1.055e6*Get('Patch','VWall',{'Lv4','Lv1','Sv1','Sv2','Lv2','Lv3'}))); % [ml/min/g]
    fillBullseye(2.25  ,0.0,0.5,0,360);      % to set color range
    fillBullseye(4.75  ,0.0,0.5,0,360);      % to set color range
    fillBullseye(cLVac,0.0,0.5,0,360);      % apical cap
    fillBullseye(cLVa ,0.5,1.0,-45,315);    % apex
    fillBullseye(cLVm ,1.0,1.5,0,360);      % mid
    fillBullseye(cLVb ,1.5,2.0,0,360);      % base
    uistack(c,'top');

    %----% Labels
    if Lab == 1
        dec = 1; % Amount of decimals
        % RV 
        ang = (0)*(pi/180); [x,y] = pol2cart(ang, -2.475);
        t = text(x,y,string(num2str(cRV','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apical cap
        ang = (60*0)*(pi/180); [x,y] = pol2cart(ang, 0);
        t = text(x,y,string(num2str(cLVac','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apex
        labels = string(num2str(cLVa','%.1f'))';
        for k = 1:length(labels)
            ang = (90*(k-1))*(pi/180); [x,y] = pol2cart(ang, 0.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Mid
        labels = string(num2str(cLVm','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.25);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Base
        labels = string(num2str(cLVb','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
    end
    
    %----% Define colormap and add colorbar
    mymap = [linspace(  1, 29, 250)'/255,linspace(  1, 62, 250)'/255,linspace(  1,129, 250)'/255;...
             linspace( 29,  9, 250)'/255,linspace( 62,141, 250)'/255,linspace(129,138, 250)'/255;...
             linspace(  9, 86, 250)'/255,linspace(141,196, 250)'/255,linspace(138, 60, 250)'/255;...
             linspace( 86,170, 250)'/255,linspace(196,213, 250)'/255,linspace( 60, 33, 250)'/255;...
             linspace(170,253, 250)'/255,linspace(213,231, 250)'/255,linspace( 33, 37, 250)'/255;...
             linspace(253,255, 250)'/255,linspace(231,255, 250)'/255,linspace( 37,255, 250)'/255];
    colormap(mymap);
    b = colorbar;
    
    %----% Figure design
%     set(b,'Ticks',[2.5:0.5:7.5],'TickLabels',{'2.5','','3.5','','4.5','','5.5','','6.5','','7.5 \mumol/min/g'},'Location','east');
    set(b,'Ticks',[2.25:0.25:4.75],'TickLabels',{'2.25','','2.75','','3.25','','3.75','','4.25','','4.75 ml/min/g'},'Location','east');
%     set(b,'Ticks',linspace(0.35,1.35,11)*4.71,'TickLabels',{'0.35','','0.55','','0.75','','0.95','','1.15','','1.35 ml/min/g'},'Location','east'); % Label for resting flow
    set(b, 'Position',[0.7664    0.24    0.0381    0.7644-2*(0.24-0.1359)]);
    title(iTitle{iPName});
    set(gcf,'Color','white');
    ax=gca; ax.XColor = 'none'; ax.YColor = 'none';
end

%% Figure 4 - Myocardial oxygen demand per tissue weight
%----% Add subfolder AHABullseye to the path
addpath('AHABullseye');

%----% Load P-struct
PName = {'PRef.mat';...
         'PRef_acuteLBBB.mat';...
         'PRef_chronicLBBB.mat'};
iTitle = {{'REFERENCE'   ,'Myocardial Oxygen Demand'},...
          {'ACUTE LBBB'  ,'Myocardial Oxygen Demand'},...
          {'CHRONIC LBBB','Myocardial Oxygen Demand'}};
for iPName = 1:3
    load(PName{iPName});

    %----% Create figure
    f=figure;
    Lab = 0;  % 0 if labels off, 1 if labels on

    %----% RV
    % Bullseye data
    [theta,rho] = getWedgeBorder(360, 2*360, 0, 2);
    wedgeHandle = polar(gca,theta,rho);
    XData = get(wedgeHandle,'XData');
    YData = get(wedgeHandle,'YData');
    XData=XData-1;

    % Color data
    cRV = 1e3*60*Get('CorArtVen','VO2','RCA')./(1.055e6*Get('Patch','VWall','Rv1')); % [umol/min/g]

    % Plot
    close(f);
    figure; hold on; axis equal;
    title('Demand')
    fillBullseyeRV(cRV,0,2,0,360);
    plot(XData,YData,'k');

    %----% LV
    % Bullseye data
    c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);

    % Color data
    cLVac = (60*1e3*Get('CorArtVen','VO2',P.CorArtVen.Name(12))./(1.055e6*Get('Patch','VWall',{'Lv12'}))); % [umol/min/g]
    cLVa  = (60*1e3*Get('CorArtVen','VO2',P.CorArtVen.Name([3,10,11,17]))./(1.055e6*Get('Patch','VWall',{'Lv11','Lv9','Sv5','Lv10'}))); % [umol/min/g]
    cLVm  = (60*1e3*Get('CorArtVen','VO2',P.CorArtVen.Name([2,8,9,14,16,4]))./(1.055e6*Get('Patch','VWall',{'Lv8','Lv5','Sv3','Sv4','Lv6','Lv7'}))); % [umol/min/g]
    cLVb  = (60*1e3*Get('CorArtVen','VO2',P.CorArtVen.Name([1,6,7,13,15,5]))./(1.055e6*Get('Patch','VWall',{'Lv4','Lv1','Sv1','Sv2','Lv2','Lv3'}))); % [umol/min/g]
    fillBullseye(1.5  ,0.0,0.5,0,360);      % to set color range
    fillBullseye(6.5  ,0.0,0.5,0,360);      % to set color range
    fillBullseye(cLVac,0.0,0.5,0,360);      % apical cap
    fillBullseye(cLVa ,0.5,1.0,-45,315);    % apex
    fillBullseye(cLVm ,1.0,1.5,0,360);      % mid
    fillBullseye(cLVb ,1.5,2.0,0,360);      % base
    uistack(c,'top');

    %----% Labels
    if Lab == 1
        dec = 1; % Amount of decimals
        % RV 
        ang = (0)*(pi/180); [x,y] = pol2cart(ang, -2.475);
        t = text(x,y,string(num2str(cRV','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apical cap
        ang = (60*0)*(pi/180); [x,y] = pol2cart(ang, 0);
        t = text(x,y,string(num2str(cLVac','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apex
        labels = string(num2str(cLVa','%.1f'))';
        for k = 1:length(labels)
            ang = (90*(k-1))*(pi/180); [x,y] = pol2cart(ang, 0.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Mid
        labels = string(num2str(cLVm','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.25);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Base
        labels = string(num2str(cLVb','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
    end
    
    %----% Define colormap and add colorbar
    mymap = [linspace(  1, 29, 250)'/255,linspace(  1, 62, 250)'/255,linspace(  1,129, 250)'/255;...
             linspace( 29,  9, 250)'/255,linspace( 62,141, 250)'/255,linspace(129,138, 250)'/255;...
             linspace(  9, 86, 250)'/255,linspace(141,196, 250)'/255,linspace(138, 60, 250)'/255;...
             linspace( 86,170, 250)'/255,linspace(196,213, 250)'/255,linspace( 60, 33, 250)'/255;...
             linspace(170,253, 250)'/255,linspace(213,231, 250)'/255,linspace( 33, 37, 250)'/255;...
             linspace(253,255, 250)'/255,linspace(231,255, 250)'/255,linspace( 37,255, 250)'/255];
    colormap(mymap);
    b = colorbar;
    
    %----% Figure design
    set(b,'Ticks',[1.5:0.5:6.5],'TickLabels',{'1.5','','2.5','','3.5','','4.5','','5.5','','6.5 \mumol/min/g'},'Location','east');
%     set(b,'Ticks',linspace(0.35,1.35,11)*4.71,'TickLabels',{'0.35','','0.55','','0.75','','0.95','','1.15','','1.35 ml/min/g'},'Location','east'); % Label for resting flow
    set(b, 'Position',[0.7664    0.24    0.0381    0.7644-2*(0.24-0.1359)]);
    title(iTitle{iPName});
    set(gcf,'Color','white');
    ax=gca; ax.XColor = 'none'; ax.YColor = 'none';
end

%% Figure 4 - Myocardial flow reserve
%----% Add subfolder AHABullseye to the path
addpath('AHABullseye');

%----% Load P-struct
PName = {'PRef.mat';...
         'PRef_acuteLBBB.mat';...
         'PRef_chronicLBBB.mat'};
iTitle = {{'REFERENCE'   ,'Myocardial Flow Reserve'},...
          {'ACUTE LBBB'  ,'Myocardial Flow Reserve'},...
          {'CHRONIC LBBB','Myocardial Flow Reserve'}};
for iPName = 1:3
    [~, Pname,~] = fileparts(PName{iPName});
    files = {strjoin({Pname,    '.mat'},'');...
             strjoin({Pname,'_Hyp.mat'},'')};
    
    %----% Gather data
    load(files{1})
    qRest = mean(P.CorArtVen.q);
    load(files{2})
    qHyp = mean(P.CorArtVen.q);
    P.CorArtVen.MFR = qHyp./qRest;

    %----% Create figure
    f=figure;
    Lab = 0;  % 0 if labels off, 1 if labels on
    
    %----% RV
    % Bullseye data
    [theta,rho] = getWedgeBorder(360, 2*360, 0, 2);
    wedgeHandle = polar(gca,theta,rho);
    XData = get(wedgeHandle,'XData');
    YData = get(wedgeHandle,'YData');
    XData=XData-1;

    % Color data
    cRV = Get('CorArtVen','MFR',P.CorArtVen.Name(18));

    % Plot
    close(f);
    figure; hold on; axis equal;
    title('MFR')
    fillBullseyeRV(cRV,0,2,0,360);
    plot(XData,YData,'k');

    %----% LV
    % Bullseye data
    c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);

    % Color data
    cLVac = (Get('CorArtVen','MFR',P.CorArtVen.Name(12)));
    cLVa  = (Get('CorArtVen','MFR',P.CorArtVen.Name([3,10,11,17])));
    cLVm  = (Get('CorArtVen','MFR',P.CorArtVen.Name([2,8,9,14,16,4])));
    cLVb  = (Get('CorArtVen','MFR',P.CorArtVen.Name([1,6,7,13,15,5])));
    fillBullseye(1.5    ,0.0,0.5,0,360);      % to set color range
    fillBullseye(6.5    ,0.0,0.5,0,360);      % to set color range
    fillBullseye(cLVac,0.0,0.5,0,360);      % apical cap
    fillBullseye(cLVa ,0.5,1.0,-45,315);    % apex
    fillBullseye(cLVm ,1.0,1.5,0,360);      % mid
    fillBullseye(cLVb ,1.5,2.0,0,360);      % base
    uistack(c,'top');

    %----% Labels
    if Lab == 1
        dec = 1; % Amount of decimals
        % RV 
        ang = (0)*(pi/180); [x,y] = pol2cart(ang, -2.475);
        t = text(x,y,string(num2str(cRV','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apical cap
        ang = (60*0)*(pi/180); [x,y] = pol2cart(ang, 0);
        t = text(x,y,string(num2str(cLVac','%.1f'))'); set(t,'HorizontalAlignment','Center');
        % LV Apex
        labels = string(num2str(cLVa','%.1f'))';
        for k = 1:length(labels)
            ang = (90*(k-1))*(pi/180); [x,y] = pol2cart(ang, 0.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Mid
        labels = string(num2str(cLVm','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.25);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
        % LV Base
        labels = string(num2str(cLVb','%.1f'))';
        for k = 1:length(labels)
            ang = (60*(k-1)+30)*(pi/180); [x,y] = pol2cart(ang, 1.75);
            t = text(x,y,labels{k}); set(t,'HorizontalAlignment','Center');
        end
    end

    %----% Define colormap and add colorbar
    mymap = [linspace(  1, 29, 250)'/255,linspace(  1, 62, 250)'/255,linspace(  1,129, 250)'/255;...
             linspace( 29,  9, 250)'/255,linspace( 62,141, 250)'/255,linspace(129,138, 250)'/255;...
             linspace(  9, 86, 250)'/255,linspace(141,196, 250)'/255,linspace(138, 60, 250)'/255;...
             linspace( 86,170, 250)'/255,linspace(196,213, 250)'/255,linspace( 60, 33, 250)'/255;...
             linspace(170,253, 250)'/255,linspace(213,231, 250)'/255,linspace( 33, 37, 250)'/255;...
             linspace(253,255, 250)'/255,linspace(231,255, 250)'/255,linspace( 37,255, 250)'/255];
    colormap(mymap);
    b = colorbar;
    
    %----% Figure design
    set(b,'Ticks',[1.5:0.5:6.5],'TickLabels',{'1.5','','2.5','','3.5','','4.5','','5.5','','6.5'},'Location','east');
    set(b, 'Position',[0.7664    0.24    0.0381    0.7644-2*(0.24-0.1359)]);
    title(iTitle{iPName});
    set(gcf,'Color','white');
    ax=gca; ax.XColor = 'none'; ax.YColor = 'none';
end

%% Graphical Abstract - SSA loops
% External work and potential energy
% Gather data
PName = {'PRef.mat';...
         'PRef_acuteLBBB.mat';...
         'PRef_chronicLBBB.mat'};
     
f=figure;

for ip = 1:length(PName)
    load(PName{ip});

    [Ef, Sf, Sfa, Efr, Sfr, ED, SSA] = calcSSA(P);
    [EDSSR_Ef,EDSSR_Sf] = calcEDSSR(P);

    %----% Create figure
    ititle = {{'Septum','(Early-activated)'},{'Left ventricular free wall','(Late-activated)'}};
    ipatch = [18, 10]; % Septum, LV free wall
    % Subplots
    for is = 1:2
        subplot(3,2,(ip-1)*2+is); hold on;

        % Active end-systolic stress strain relation
        plot(Efr(:,ipatch(is)),Sfr(:,ipatch(is))*1e-3,'LineWidth',3,'Color',[0 0 0]);

        % Passive end-diastolic stress strain relation
        [~,ind] = min(abs(EDSSR_Ef-Ef(ED(ipatch(is)),ipatch(is))));
        plot(EDSSR_Ef(1:ind),EDSSR_Sf(1:ind,ipatch(is))*1e-3,'LineWidth',3,'Color',[0 0 0]);

        % Stress-strain loop
        plot(Ef(:,ipatch(is)),Sf(:,ipatch(is))*1e-3,'LineWidth',3,'Color',[0 0 0]);

    %     % End-systolic point
    %     plot(Efr(1,ipatch(is)),Sfr(1,ipatch(is))*1e-3,'p','Color',[ 47  85 151]/255,'MarkerFaceColor',[ 47  85 151]/255,'MarkerSize',16);
    % 
    %     % End-diastolic point
    %     plot(Ef(ED(ipatch(is)),ipatch(is)),P.Patch.Sf(ED(ipatch(is)),ipatch(is))*1e-3,'x','Color',[ 47  85 151]/255,'LineWidth',3,'MarkerSize',16);

        % Subplot design
        xlim([-0.3 0.2]); xticks([]); %xlabel('Strain (-)')
        ylim([0 70]); yticks([]); %ylabel('Stress (kPa)');
        if ip == 1
            title(ititle{is});
        end
        set(gca,'LineWidth',2,'FontSize',18);
    end
end
% Figure design
set(gcf,'Color','white');
set(f,'Position',[204.2000   51.0000  557.8000  910.0000]);

%% Graphical Abstract - Bars
% External work and potential energy
% Gather data
PName = {'PRef.mat';...
         'PRef_acuteLBBB.mat';...
         'PRef_chronicLBBB.mat'};
     
f=figure;
for ip = 1:length(PName)
    load(PName{ip});

    VO2        = P.CorArtVen.VO2;       % [mmol/s]]    Oxygen consumption
    SSA2Cor    = P.CorArtVen.SSA2Cor;   % [J]          Stress-Strain Area (SSA)
    k1         = 4.94e-3;               % [mmol/J]
    k2         = 24.2;                  % [mmol/m^3/s]
    k2VO2      = k2*sum(P.CorArtVen.Cor2Patch.*P.Patch.VWall,2)'; % [mmol/s]
    for ipa = 1:P.Patch.n
        EW(ipa) = -trapz([P.Patch.Ef(:,ipa)],[P.Patch.Sf(:,ipa)])*P.Patch.VWall(ipa);
    end
    
    SSAVO2     = k1*SSA2Cor/P.General.tCycle;                   % [mmol/s] SSA part of the total oxygen consumption
    EWVO2      = k1*sum(EW.*P.CorArtVen.Cor2Patch,2)'/P.General.tCycle; % [mmol/s] EW  part of the total oxygen consumption (Exteral Work)
    PEVO2      = SSAVO2 - EWVO2;                                % [mmol/s] PE  part of the total oxygen consumption (Potential Energy)
    BMVO2      = k2VO2;                                         % [mmol/s] BM  part of the total oxygen consumption (Basal Metabolism)
    
    VWall2Cor  = sum(P.Patch.VWall.*P.CorArtVen.Cor2Patch,2)';
    EWVO2      = 1e3*60*EWVO2./(1.055e6*VWall2Cor);                 % [umol/min/g] EW  part of the total oxygen consumption (Exteral Work)
    PEVO2      = 1e3*60*PEVO2./(1.055e6*VWall2Cor);                 % [umol/min/g] PE  part of the total oxygen consumption (Potential Energy)
    BMVO2      = 1e3*60*BMVO2./(1.055e6*VWall2Cor);                 % [umol/min/g] BM  part of the total oxygen consumption (Basal Metabolism)
        
    %----% Create figure
    ititle = {{'Septum','(Early-activated)'},{'Left ventricular free wall','(Late-activated)'}};
    icor = [14, 2]; % Septum, LV free wall coronary territory
    ipatch = [18, 10]; % Septum, LV free wall

    subplot(1,3,ip); hold on;
    bar([BMVO2(icor);...
         EWVO2(icor);...
         PEVO2(icor)]','stacked');

    % Subplot design
        xlim([0 3]); xticks([1 2]); %xlabel('Strain (-)')
        ylim([0 6]); %yticks([]); %ylabel('Stress (kPa)');
%     if ip == 1
%         title(ititle{is});
%     end
    set(gca,'LineWidth',2,'FontSize',18);
    
end
% Figure design
set(gcf,'Color','white');
set(f,'Position',[459.4000  234.6000  557.6000  236.0000]);

%% Functions
%------%
function [Ef, Sf, Sfa, Efr, Sfr, ED, SSA] = calcSSA(P)
    % Init
    Lsi0Act = P.Patch.Lsi0Act;
    LsRef   = P.Patch.LsRef;
    AmRef   = P.Patch.AmRef;
    VWall   = P.Patch.VWall;
    Ef      = P.Patch.Ef;
    Sf      = P.Patch.Sf;
    Sfa     = P.Patch.Sf - P.Patch.SfEcm - P.Patch.SfTit;

    % Start point
    Am0     = ((Lsi0Act./LsRef).^2).* AmRef;
    T0      = 0 .* Am0;
    Ef0     = 0.5 .* log(Am0./AmRef);
    Sf0     = 0 .* Ef0;

    % End-diatole point
    ED      = ceil(P.Patch.Depolarization(1,:)./P.General.Dt);
    ED(ED==0) = 1;
    Efed    = diag(P.Patch.Ef(ED,:))';
    Sfed    = diag(P.Patch.Sf(ED,:)-P.Patch.SfEcm(ED,:)-P.Patch.SfTit(ED,:))';

    % End-systole point
    [~,ES]  = max(P.Patch.C);
    Efes    = diag(P.Patch.Ef(ES,:))';
    Sfes    = diag(P.Patch.Sf(ES,:)-P.Patch.SfEcm(ES,:)-P.Patch.SfTit(ES,:))';
    Ames    = diag(P.Patch.Am(ES,:))';
%     Tes     = diag(P.Patch.T(ES,:))';
    Tes     = Sfes .* VWall ./ (2 * Ames);

    % End-systole to start point
    for ip = 1:P.Patch.n
        ATrA(:,ip) = linspace(Ames(ip), Am0(ip), 1000)';
        ATrT(:,ip) = linspace( Tes(ip),  T0(ip), 1000)';
    end
    Efr     = 0.5 .* log(ATrA ./ AmRef);
    Sfr     = 2 * ATrA .* ATrT ./ VWall;

    % SSA
    for ip = 1:P.Patch.n
        SSA(ip) = -trapz([Ef0(ip); Efes(ip);  Ef(ED(ip):ES(ip),ip); Efr(:,ip)],...
                         [Sf0(ip);  Sf0(ip); Sfa(ED(ip):ES(ip),ip); Sfr(:,ip)]);
    end
end
%------%

%------%
function [EDSSR_Ef,EDSSR_Sf] = calcEDSSR(P)
    % Retrieve nonlinear passive stress-strain behavior
    nSteps = 1000;
    Efss   = linspace(-0.2, 0.2, nSteps);
    for iss = 1:nSteps
        if isfield(P.Patch, 'k1'), k1 = Get('Patch','k1','Lv1');
        else,                      k1= 10; end
        LsRef = Get('Patch','LsRef','All');
        SfAct = Get('Patch','SfAct','All');
        Ls0Pas= Get('Patch','Ls0Pas','All');
        dLsPas= Get('Patch','dLsPas','All');
        SfPas = Get('Patch','SfPas','All');
        k2    = 0.01; 
        kk3   = 2*(LsRef./dLsPas);
        LfP   = bsxfun(@times,exp(Efss(iss)), LsRef./Ls0Pas);
        yEcm  = LfP.^k1;
        SfEcm = bsxfun(@times, yEcm-1, 0.0349*SfPas);% Ls=Ls0Pas, zero stress
        y     = bsxfun(@power,LfP,kk3);
        SfTit = bsxfun(@times, y-1, k2*SfAct);% titin is softer than ecm, proportional with SfAct
        SfPasT(iss,:) = SfTit + SfEcm;
    end
    EDSSR_Ef = Efss;
    EDSSR_Sf = SfPasT;
end
