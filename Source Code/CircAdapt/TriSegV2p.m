function TriSegV2p
% function TriSegV2p
% TriSeg is a 3-wall structure (Left,Septal,Right) with 2 cavities (R,L)
% Calculates: cavity volumes V -> dimensions of the 'double bubble',
% myofiber stress Sf, wall tension T and cavity pressures p
% VS and YS repesent septal volume displacement and junction radius.
% State variables P.TriSeg V and Y represent initial estimates to calculate
% VS and YS accurately.
% Theo Arts, Maastricht University, Oct 13, 2012

global P;
if P.TriSeg.n==0 %if there is no chamber
    return
end

TriSeg  = P.TriSeg;
Wall    = P.Wall;
Tau     = P.General.dt; % lowpass for VS,YS 1st estimate
RhoB    = P.General.RhoB; % blood density
iCavity = TriSeg.iCavity; %related cavities
% VLT     = TriSeg.VL   ; % left cavity volume
% VRT     = TriSeg.VR   ; % right cavity volume
n       = TriSeg.n    ; % number of TriSeg's
iWall   = TriSeg.iWall; %related walls
VT      = TriSeg.V    ; % init septal cap volume, to be solved
YT      = TriSeg.Y    ; % init radius junction circle, to be solved
Am0W    = Wall.Am0    ; %zero stress wall area
DADTW   = Wall.DADT   ; %wall compliance
VWallW  = Wall.VWall  ; % 3 wall volumes

for i=1:n %for all TriSeg's
    iC   = iCavity(:,i)+(0:1); % 2 cavities
    iW   = iWall(:,i)  +(0:2); % 3 walls
    Am0  = Am0W(:,iW) ; %zero stress wall area
    DADT = DADTW(:,iW); %wall compliance
    VWall= VWallW(iW) ; % 3 wall volumes
%     nt   = size(Am0,1); %number of time points
    VWL  = mean(VWall(1:2)); % enclosed wall volume 1st cavity
    VWR  = mean(VWall(2:3)); % enclosed wall volume 2nd cavity
    VLT  = P.Cavity.V(:,iC(1)); % left cavity volume
    VRT  = P.Cavity.V(:,iC(2)); % right cavity volume
    VL   = VLT(:,i)+VWL; % midwall enclosed LV-volume
    VR   = VRT(:,i)+VWR; % midwall enclosed RV-volume
    VS= TriSeg.V(:,i); % septal rightward volume shift
    YS= TriSeg.Y(:,i); % radius LV-RV junction
%     VCav = P.Cavity.V(:,iC);
%     V    = max(0,bsxfun(@plus,VCav,[VWL,VWR])); %2 midwall volumes(t)

    VS0=VS; % storage initial value septal cap volume
    YS0=YS; % storage initial value junction radius
    
    % Calculation iteration VS and YS increment for better equilibrium 
    for j=1:2 % 1-2 iterations appear to be sufficient
        [dVS,dYS]=DvDp(VS,YS,VL,VR,Am0,DADT);
        VS=VS+dVS; % Improved estimate of TriSeg geometry
        YS=YS+dYS; % only one iteration has been applied
    end  
    
    % Limitation of time-derivatives to min and max value
    % for better numerical safety
    AL= Am0(:,1); %zero stress wall area L
    AS= Am0(:,2); %zero stress wall area S
    AR= Am0(:,3); %zero stress wall area R
    F1= 1.08; % maximum allowed stretch from zero load state
    F3= F1^3;
    yLo= sqrt(AS./(1+AS./min(AL,AR))/pi)/F1;% assume sphere AL+AS
    yHi= F1*sqrt(AS./(1+AS./(max(AL,AR)*F1^2))/pi);
    vLo= -F3*(2*pi/3)*(AL/(2*pi)).^1.5+VL; % given wall area -> max volume
    vHi=  F3*(2*pi/3)*(AR/(2*pi)).^1.5-VR; % given wall area -> max volume
    YS= max(yLo,min(yHi,YS)); % avoid too large steps, numerical safety
    VS= max(vLo,min(vHi,VS));
        
    % final TriSeg geometry requires some recalculations
    VM   = [VS-VL,VS,VS+VR]; %1st estimate cap-volumes
    ARef = pi*YS.^2; % normalization reference area
    V    = bsxfun(@rdivide,VM  ,ARef.*YS); % normalized capvolumes
    A0   = bsxfun(@rdivide,Am0 ,ARef    ); % normalized zero stress area
    dTdA = bsxfun(@rdivide,ARef,DADT); % area normalized wall stiffnes
    % Auxilary normalized variables
    X    = VdPi2X(V); % normalized cap-height
    XX   = X.^2;
    RR   = XX+1; % normalized wall area
    C    = 2*X./RR; % normalized wall curvature
    
    % Revert normalization
    T    = max(0,(RR-A0).*dTdA);% [N/m] wall tension
    Am   = bsxfun(@times,ARef,RR); % wall area
    Cm   = bsxfun(@rdivide,C,YS); % wall curvature
    pTrans = 2*Cm.*T; % transmural pressure with effect of buckling
    
    % Cavity impedance properties, needed to make node connection
    Vm  = [VL+VWL,VR+VWR];
    Len = 2*Vm.^(1/3); % cavity length
    A   = Vm ./ Len; % cross-sectional area
    dTdA1= bsxfun(@rdivide,dTdA(:,[1,3]),ARef);    
    Z0  = 0.5*sqrt(RhoB .* dTdA1 ./abs(A).^1.5); %0.5*Cavity wave impedance
    
    % State variable derivatives
    TriSeg.VDot(:,i)= (VS-VS0)/Tau; % V serves as initial estimate
    TriSeg.YDot(:,i)= (YS-YS0)/Tau; % Y serves as initial estimate
   
    % writing geometric data TriSeg to be used for output, not for solving
    TriSeg.VS(:,i)= VS; % final solution Rv-Sv-Lv
    TriSeg.YS(:,i)= YS; % volume of septal cap, only used for output

    % writing wall data
    Wall.Am(:,iW)    = Am ; % wall area
    Wall.Cm(:,iW)    = Cm ; % wall curvature
    Wall.T(:,iW)     = T  ; % wall tension
    Wall.pTrans(:,iW)= bsxfun(@minus,pTrans,mean(pTrans,2)); % pTrans
    
%     TriSeg.AL(:,iTr)= A(:,1) ; % cross-sectional area LV
%     TriSeg.AR(:,iTr)= A(:,2) ; % cross-sectional area RV
%     TriSeg.ZL(:,iTr)= Z0(:,1); % source impedance LV
%     TriSeg.ZR(:,iTr)= Z0(:,2); % source impedance RV
    P.Cavity.A(:,iC)= A;% [A(:,1) A(:,2)]; % cross-sectional area
    P.Cavity.Z(:,iC)= Z0;% [Z0(:,1) Z0(:,2)]; % source impedance
    
    % anti-collapse counter pressure for numerical safety
    p0  = 0.2*P.General.p0;
    eps = 0.1;
    VNLo= max(eps,bsxfun(@rdivide,[VL,VR],[VWL,VWR])-1);
    dpLo= bsxfun(@times,p0,max(0,0.3./VNLo-1).^2);% anti-collapse safety
%     TriSeg.pTransL(:,iTr)= -pTrans(:,1)-dpLo(:,1); % LV transmural pressure
%     TriSeg.pTransR(:,iTr)=  pTrans(:,3)-dpLo(:,2); % RV transmural pressure
    P.Cavity.pTrans(:,iC)= [-pTrans(:,1)-dpLo(:,1),pTrans(:,3)-dpLo(:,2)];
    
end

P.Wall = Wall;
P.TriSeg=TriSeg;
end

%=== Additional function ============ 
function [dV,dY]=DvDp(VS,YS,VL,VR,Am0,DADT)
% input: VS=septal cap volume, YS= junction radius
% [VL, VR]= [left,right] midwall enclosed volumes
% Aw0,DADT= [zero load wall area, compliance wall] for L- S- R- walls
VM   = [VS-VL,VS,VS+VR]; % cap volumes LSR
ARef = pi*YS.^2; % reference area for normalization
V    = bsxfun(@rdivide,VM  ,ARef.*YS); %normalized cap volumes
A0   = bsxfun(@rdivide,Am0 ,ARef    ); %normalized wall areas
dTdA = bsxfun(@rdivide,ARef,DADT); % normalized wall compliance

X    = VdPi2X(V); % solves 3rd order polynomial analytically
XX   = X.^2;
RR   = XX+1; %wall area (normalized to junction area (pi y^2)
RR2  = RR.^2;

dA = RR-A0; % wall area stretch increment
T  = max(0,dA.*dTdA);% buckle possibility

Av =  4*X  ./RR; % partial derivative dA/dVS, normalized
Ap = (1-XX)./RR; % partial derivative dA/dYS, normalized

Avv= 8*Ap./RR2; % 2nd order partial derivatives
Avp= -2*Av./RR2;
App= 3*V.*Av./RR2;

Fv= sum(T.*Av,2); %2*Tx, dE/dVS ~force
Fp= sum(T.*Ap,2); %-Ty, dE/dYS ~force

Fvv=sum(bsxfun(@times,dTdA,Av.^2 )+T.*Avv,2); % 2nd derivative,'stiffness'
Fvp=sum(bsxfun(@times,dTdA,Av.*Ap)+T.*Avp,2); % 2nd derivative,'stiffness' 
Fpp=sum(bsxfun(@times,dTdA,Ap.^2 )+T.*App,2); % 2nd derivative,'stiffness'

Det=Fvv.*Fpp-Fvp.*Fvp; % determinant
dv= -( Fv.*Fpp-Fp.*Fvp)./Det; %solution of 2 parameters by matrix inversion
dp= -(-Fv.*Fvp+Fp.*Fvv)./Det;
dV=  YS.^3 .* dv; % reversion of normalization
dY= 0.5*YS .* dp;
end

function X=VdPi2X(VdPi)
% Analytical solution of X as a function of cap volume V
% VdPi = normalized Vcap/(pi*Rjunction^3)
% assumption: cap junction radius is normalized to 1
% Solving 3rd order polynomial
% Mathematica: Solve[x^3 + 3x  - 2V == 0, x]
% Q= (V + Sqrt(V^2 + 1))^(1/3);  x= Q - 1/Q;
V = 3*VdPi;
Q = (V + sqrt(V.^2+1)).^(1/3);
X= (Q - 1./Q);
end

% function TriSegV2p
% % function TriSegV2p
% % TriSeg is a 3-wall structure (Left,Septal,Right) with 2 cavities (R,L)
% % Calculates: cavity volumes V -> dimensions of the 'double bubble',
% % myofiber stress Sf, wall tension T and cavity pressures p
% % VS and YS repesent septal volume displacement and junction radius.
% % State variables P.TriSeg V and Y represent initial estimates to calculate
% % VS and YS accurately.
% % Theo Arts, Maastricht University, Oct 13, 2012
% 
% global P;
% if P.TriSeg.n==0 %if there is no chamber
%     return
% end
% 
% TriSeg=P.TriSeg;
% 
% n  = TriSeg.n; % number of TriSeg's
% iCavity = TriSeg.iCavity; %related cavities
% iWall   = TriSeg.iWall  ; %related walls
% RhoB    = P.General.RhoB; % blood density
% 
% for i=1:n %for all TriSeg's
%     iC   = iCavity(:,i)+(0:1); % 2 cavities
%     iW   = iWall(:,i)  +(0:2); % 3 walls
%     Am0  = P.Wall.Am0(:,iW); %zero stress wall area
%     DADT = P.Wall.DADT(:,iW); %wall compliance
%     nt   = size(Am0,1); %number of time points
%     VWall= P.Wall.VWall(iW); % 3 wall volumes
%     VWL  = mean(VWall(1:2)); % enclosed wall volume 1st cavity
%     VWR  = mean(VWall(2:3)); % enclosed wall volume 2nd cavity
%     VCav = P.Cavity.V(:,iC);
%     V    = max(0,bsxfun(@plus,VCav,[VWL,VWR])); %2 midwall volumes(t)
% 
%     VS= TriSeg.V(:,i);
%     YS= TriSeg.Y(:,i); % 1st estimate of [V,Y]
%     YRef= mean(YS); VRef= YRef^3;
%     dvR = 0.02;     dyR = dvR;
%     dv=dvR*VRef;
%     dy=dyR*YRef; % increments for d[Txy]/dVY
% 
%     as=sin(pi/12); ac=cos(pi/12); bs=sin(pi/4);
%     % coefficients to for a equilateral triangle in 2D
%     [Tx0a,Ty0a]=VY2Txy(VS-bs*dv,YS+bs*dy,Am0,DADT,V); % tension Txy =fu(VY)
%     [TxVa,TyVa]=VY2Txy(VS+ac*dv,YS+as*dy,Am0,DADT,V); % partial V derivative
%     [TxYa,TyYa]=VY2Txy(VS-as*dv,YS-ac*dy,Am0,DADT,V); % partial Y derivative
%     % matrix of coefficxients to determine partial derivatives
%     DX=[-bs  bs
%          ac  as
%         -as -ac];
%     Ddvy= pinv(DX')*diag([1/dv,1/dy]);
%     DTx= [Tx0a,TxVa,TxYa]*Ddvy; % [dTx/dx,dTx/dy]
%     DTy= [Ty0a,TyVa,TyYa]*Ddvy; % [dTy/dx,dTy/dy]
% 
%     DET= DTx(:,1).*DTy(:,2)-DTx(:,2).*DTy(:,1); % determinant
%     Tx0=mean([Tx0a,TxVa,TxYa],2);
%     Ty0=mean([Ty0a,TyVa,TyYa],2);
%     dV=(-DTy(:,2).*Tx0+DTx(:,2).*Ty0)./DET;
%     dY=(+DTy(:,1).*Tx0-DTx(:,1).*Ty0)./DET;
%     VS=VS+dV;
%     YS=YS+dY;
%     Tau=TriSeg.Tau;
%     TriSeg.VDot(:,i)= dV/Tau; % keep track of solution
%     TriSeg.YDot(:,i)= dY/Tau;
% 
%     for j=1:1 % extra iterations for solution of TriSeg geometry
%         [Tx0,Ty0]=VY2Txy(VS,YS,Am0,DADT,V);
%         dV=(-DTy(:,2).*Tx0+DTx(:,2).*Ty0)./DET;
%         dY=(+DTy(:,1).*Tx0-DTx(:,1).*Ty0)./DET;
%         VS=VS+dV;
%         YS=YS+dY;
%     end
%     
%     % writing geometric data TriSeg
%     TriSeg.YS(:,i)=YS; % final solution Rv-Sv-Lv junction radius
%     TriSeg.VS(:,i)=VS; % volume of septal cap
%     
%     % writing wall data
%     [Tx0,Ty0,Am,Cm,T]=VY2Txy(VS,YS,Am0,DADT,V); % solution TriSeg
%     P.Wall.Am(:,iW)= max(Am0,Am); % wall area
%     P.Wall.Cm(:,iW)= Cm; % wall curvature
%     P.Wall.T(:,iW) = max(0,T) ; % wall tension
%     pTrans= 2*Cm.*T; % transmural pressure
%     P.Wall.pTrans(:,iW)=pTrans;
%      
%     % Cavity impedance properties, needed to make node connection +++++++
%     % efficienter maken
%     Vw= repmat([VWL,VWR],[nt,1]);
%     Vm= V + Vw;
%     Len= 2*Vm.^(1/3);
%     A  = ( V + 0.1*Vw ) ./Len;
%     Z0 = sqrt(RhoB ./ abs(A.*DADT(:,[1,3]).*Len)); %Cavity wave impedance
%     P.Cavity.A(:,iC) = A; % cross-sectional area for valve inflow and outflow pressure
%     P.Cavity.Z(:,iC) = Z0;
%     P.Cavity.pTrans(:,iC)= [-pTrans(:,1),pTrans(:,3)];
%     
% end
% 
% P.TriSeg=TriSeg;
% end
% function [Tx,Ty,Am,Cm,Tm]=VY2Txy1(VS,YS,Am0,DADT,VLR)
% % 1st order approximation of TriSeg according to J. Lumens et al.
% % For a wall with zero-stress area Am0 and compliance DADT
% % cap volume VS, junction radius YS and cavity volumes VLR=[VL,VR]
% % Result: Summed axial and radial tension components [Tx,Ty] on junction
% % for 3 walls: Am= midwall area, Cm= midwall curvature, Tm= wall tension
% ARef=YS.^2;
% VRef=YS.^3;
% VS  = VS./VRef;
% Am0 = bsxfun(@rdivide,Am0 ,ARef);
% DADT= bsxfun(@rdivide,DADT,ARef);
% VLR = bsxfun(@rdivide,VLR ,VRef);
% Vm  = [VLR,VS] * [...
%     -1     0     0
%      0     0     1
%      1     1     1];
% % Solving 3rd order polynomial
% % Mathematica: Solve[x^3 + 3x  - 2V == 0, x]
% % Q= (V + Sqrt(V^2 + 1))^(1/3);  x= Q - 1/Q;
% SignVm= sign(Vm); Vm=abs(Vm);
% V     = (3/pi)*Vm;
% Q     = (V + sqrt(V.^2 + 1)).^(1/3);
% Xm    = SignVm .* ( Q - 1./Q );
% 
% %calculate midwall area Am and curvature Cm=1/rm
% X2    = Xm.^2;
% R2    = X2+1;
% Am    = pi*R2; % midwall cap area, buckling with T<0%+++++
% Cm    = 2*Xm./R2; % midwall cap curvature
% 
% % calculation of tension T and components Tx, Ty
% Tm  = (Am-Am0)./DADT;
% Sin = Cm;
% Cos = (1-X2)./R2;
% Txi = Cos.*Tm; %
% Tyi = Sin.*Tm;
% Tx  = sum(Txi,2); % normalized axial tension component
% Ty  = sum(Tyi,2); % normalized radial tension component
% Am  = bsxfun(@times  ,Am,ARef);
% Cm  = bsxfun(@rdivide,Cm,YS  );
% end
% 
% function [Tx,Ty,Am,Cm,Tm]=VY2Txy(VS,YS,Am0,DADT,VLR)
% % 1st order approximation of TriSeg according to J. Lumens et al.
% % For a wall with zero-stress area Am0 and compliance DADT
% % cap volume VS, junction radius YS and cavity volumes VLR=[VL,VR]
% % Result: Summed axial and radial tension components [Tx,Ty] on junction
% % for 3 walls: Am= midwall area, Cm= midwall curvature, Tm= wall tension
% Vm= [VLR,VS]* [...
%     -1     0     0
%      0     0     1
%      1     1     1];
% Ym=[YS,YS,YS];
% 
% % Solving 3rd order polynomial
% % Mathematica: Solve[x^3 + 3y^2x  - 2V == 0, x]
% % Q= (V + Sqrt(V^2 + y^6))^(1/3);  x= Q - y^2/Q;
% SignVm= sign(Vm); Vm=abs(Vm);
% V     = (3/pi)*Vm;
% Q     = (V + sqrt(V.^2 + Ym.^6)).^(1/3);
% Xm    = SignVm .* ( Q - Ym.^2 ./ Q );
% 
% %calculate midwall area Am and curvature Cm=1/rm
% X2    = Xm.^2; Y2= Ym.^2;
% R2    = X2+Y2;
% % Am    = pi*R2; % midwall cap area
% Am    = pi*R2; % midwall cap area, buckling with T<0
% Cm    = 2*Xm./R2; % midwall cap curvature
% 
% % calculation of tension T and components Tx, Ty
% Tm=(Am-Am0)./DADT;
% Sin= Ym.*Cm;
% Cos= (Y2-X2)./R2;
% Txi = Cos.*Tm; %
% Tyi = Sin.*Tm;
% TRef=sqrt(sum(Tm.^2,2)); 
% Tx= sum(Txi,2)./TRef; % axial tension component
% Ty= sum(Tyi,2)./TRef; % radial tension component
% 
% end
% 
% function [Tx,Ty,Am,Cm,Tm]=VY2Txy1(VS,YS,Am0,DADT,VLR)
% % 1st order approximation of TriSeg according to J. Lumens et al.
% % For a wall with zero-stress area Am0 and compliance DADT
% % cap volume VS, junction radius YS and cavity volumes VLR=[VL,VR]
% % Result: Summed axial and radial tension components [Tx,Ty] on junction
% % for 3 walls: Am= midwall area, Cm= midwall curvature, Tm= wall tension
% ARef=YS.^2;
% VRef=YS.^3;
% VS  = VS./VRef;
% Am0 = bsxfun(@rdivide,Am0 ,ARef);
% DADT= bsxfun(@rdivide,DADT,ARef);
% VLR = bsxfun(@rdivide,VLR ,VRef);
% Vm  = [VLR,VS] * [...
%     -1     0     0
%      0     0     1
%      1     1     1];
% % Solving 3rd order polynomial
% % Mathematica: Solve[x^3 + 3x  - 2V == 0, x]
% % Q= (V + Sqrt(V^2 + 1))^(1/3);  x= Q - 1/Q;
% SignVm= sign(Vm); Vm=abs(Vm);
% V     = (3/pi)*Vm;
% Q     = (V + sqrt(V.^2 + 1)).^(1/3);
% Xm    = SignVm .* ( Q - 1./Q );
% 
% %calculate midwall area Am and curvature Cm=1/rm
% X2    = Xm.^2;
% R2    = X2+1;
% Am    = pi*R2; % midwall cap area, buckling with T<0%+++++
% Cm    = 2*Xm./R2; % midwall cap curvature
% 
% % calculation of tension T and components Tx, Ty
% Tm  = (Am-Am0)./DADT;
% Sin = Cm;
% Cos = (1-X2)./R2;
% Txi = Cos.*Tm; %
% Tyi = Sin.*Tm;
% Tx  = sum(Txi,2); % normalized axial tension component
% Ty  = sum(Tyi,2); % normalized radial tension component
% Am  = bsxfun(@times  ,Am,ARef);
% Cm  = bsxfun(@rdivide,Cm,YS  );
% end
% 
% function [Tx,Ty,Am,Cm,Tm]=VY2Txy(VS,YS,Am0,DADT,VLR)
% % 1st order approximation of TriSeg according to J. Lumens et al.
% % For a wall with zero-stress area Am0 and compliance DADT
% % cap volume VS, junction radius YS and cavity volumes VLR=[VL,VR]
% % Result: Summed axial and radial tension components [Tx,Ty] on junction
% % for 3 walls: Am= midwall area, Cm= midwall curvature, Tm= wall tension
% Vm= [VLR,VS]* [...
%     -1     0     0
%      0     0     1
%      1     1     1];
% Ym=[YS,YS,YS];
% 
% % Solving 3rd order polynomial
% % Mathematica: Solve[x^3 + 3y^2x  - 2V == 0, x]
% % Q= (V + Sqrt(V^2 + y^6))^(1/3);  x= Q - y^2/Q;
% SignVm= sign(Vm); Vm=abs(Vm);
% V     = (3/pi)*Vm;
% Q     = (V + sqrt(V.^2 + Ym.^6)).^(1/3);
% Xm    = SignVm .* ( Q - Ym.^2 ./ Q );
% 
% %calculate midwall area Am and curvature Cm=1/rm
% X2    = Xm.^2; Y2= Ym.^2;
% R2    = X2+Y2;
% % Am    = pi*R2; % midwall cap area
% Am    = pi*R2; % midwall cap area, buckling with T<0
% Cm    = 2*Xm./R2; % midwall cap curvature
% 
% % calculation of tension T and components Tx, Ty
% Tm=(Am-Am0)./DADT;
% Sin= Ym.*Cm;
% Cos= (Y2-X2)./R2;
% Txi = Cos.*Tm; %
% Tyi = Sin.*Tm;
% TRef=sqrt(sum(Tm.^2,2)); 
% Tx= sum(Txi,2)./TRef; % axial tension component
% Ty= sum(Tyi,2)./TRef; % radial tension component
% 
% end
