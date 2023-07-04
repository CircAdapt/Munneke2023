function CorArtVenV2p
%function CorArtVenV2p
% Art/Ven hemodynamics with peripheral resistance in between
% Flow waves enter Art and Ven.
% volume V-> transmural pressure pTrans(V) and wave impedance Z(V)
% Theo Arts, Maastricht University, Oct 12, 2012

global P;

CorArtVen = P.CorArtVen; % CorArtVen structure
nCor      = CorArtVen.n;
RhoB   = P.General.RhoB;
% indices referring to related cavities and walls
iCavity= [CorArtVen.iCavity; CorArtVen.iCavity+1]; iCavity=iCavity(:);
iCavArt= CorArtVen.iCavity;
iCavVen= CorArtVen.iCavity+1;
iCorMyoArt  = 1:8:nCor*8;
iCorMyoArt2P= 2:8:nCor*8;
iCorMyoArt2M= 3:8:nCor*8;
iCorMyoArt2N= 4:8:nCor*8;
iCorMyoVen2P= 5:8:nCor*8;
iCorMyoVen2M= 6:8:nCor*8;
iCorMyoVen2N= 7:8:nCor*8;
iCorMyoVen  = 8:8:nCor*8;

%--- Arterial and venous resistance and compliance
VArt   = P.Cavity.V(:,iCavArt); % arterial cavity volume (state variable)
VVen   = P.Cavity.V(:,iCavVen); % venous cavity volume (state variable)
VMyoArt= P.CorArtVen.VMyo(:,iCorMyoArt); % secret myocardial arteriole cavity volume (state variable)
VMyoArt2P= P.CorArtVen.VMyo(:,iCorMyoArt2P); % secret myocardial venule cavity volume (state variable)
VMyoArt2M= P.CorArtVen.VMyo(:,iCorMyoArt2M); % secret myocardial venule cavity volume (state variable)
VMyoArt2N= P.CorArtVen.VMyo(:,iCorMyoArt2N); % secret myocardial venule cavity volume (state variable)
VMyoVen2P= P.CorArtVen.VMyo(:,iCorMyoVen2P); % secret myocardial arteriole cavity volume (state variable)
VMyoVen2M= P.CorArtVen.VMyo(:,iCorMyoVen2M); % secret myocardial arteriole cavity volume (state variable)
VMyoVen2N= P.CorArtVen.VMyo(:,iCorMyoVen2N); % secret myocardial arteriole cavity volume (state variable)
VMyoVen= P.CorArtVen.VMyo(:,iCorMyoVen); % secret myocardial venule cavity volume (state variable)
V      = [VArt VMyoArt VMyoArt2P VMyoArt2M VMyoArt2N VMyoVen2P VMyoVen2M VMyoVen2N VMyoVen VVen]; 
% arterial myo venous length are equal
Len    = CorArtVen.Len([1,1,1,1,1,1,1,1,1,1],:)'; Len=Len(:)'; % repesentative length of blood vessels
p0     = CorArtVen.p0'; p0 = p0(:)'; % working pressure
if isfield(P.CorArtVen, 'facDil'), facDil = P.CorArtVen.facDil;
else,                              facDil = 1; end
A0     = (facDil.*CorArtVen.A0)'; A0 = A0(:)'; 
Aw     = CorArtVen.AWall'; Aw = Aw(:)'; % Wall cross-section
k      = CorArtVen.k'; k = k(:)'; % stiffness parameter of fibers in vessel wall

p0Row  = reshape(p0,[1,numel(p0)]);
A0Row  = reshape(A0,[1,numel(A0)]);
AwRow  = reshape(Aw,[1,numel(A0)]);
kRow   = reshape(k ,[1,numel(A0)]);
A      = max(1e-10,bsxfun(@rdivide,V,Len)); % vessel cross-section, safety A>0

% Transmural pressure and wave velocity/impedance= fu(cross-section)
% Lo represents negative transmural pressure component for A/AWall<1.0
eps  = 0.1; % clipped to minimum A/AWall value
AN   = max(bsxfun(@rdivide,A,AwRow),eps); % norm to wall cross-section
A0N  = A0Row./AwRow; % norm of reference A0 to wall cross-section
m    = kRow/3-1; % stifness exponent
Lo   = max(0,1./AN-1); % if Lo>0, negative compression pressure occurs
pHiN = bsxfun(@power,bsxfun(@rdivide,AN+0.5,A0N+0.5),m);% pressure Hi
dpHi = bsxfun(@times,pHiN,m); % diff pHiN2 to AN
dpLo = -(1./A0).*(A./A0).^(.7*(Aw./(A0))-1.4).*(.7*(Aw./(A0))-.4);
pTrans=bsxfun(@times,p0Row,pHiN+min(0,(-(A./(A0)).^(.7*(Aw./(A0))-.4))+1)); % transmural pressure, new collapsible tube law
dpdA  = bsxfun(@times,p0,dpHi./Aw+dpLo); 
c0    = sqrt(0.25+dpdA.*A/RhoB); % zero flow wave velocity
Z0    = RhoB * c0 ./ A; % wave impedance with flow=0

P.Cavity.pTrans(:,[iCavArt,iCavVen])= [pTrans(:,1:nCor) ,pTrans(:,end-nCor+1:end)]; % transmural pressure
P.Cavity.Z(:,[iCavArt,iCavVen])     = [    Z0(:,1:nCor) ,    Z0(:,end-nCor+1:end)]; % wave impedance
P.Cavity.A(:,[iCavArt,iCavVen])     = [     A(:,1:nCor) ,     A(:,end-nCor+1:end)]; % cross-section

P.CorArtVen.pTransMyo = zeros(size(pTrans(:,nCor+1:9*nCor)));
P.CorArtVen.pTransMyo(:,[iCorMyoArt,iCorMyoArt2P,iCorMyoArt2M,iCorMyoArt2N,iCorMyoVen2P,iCorMyoVen2M,iCorMyoVen2N,iCorMyoVen])= [pTrans(:,nCor+1:2*nCor), pTrans(:,2*nCor+1:3*nCor), pTrans(:,3*nCor+1:4*nCor), pTrans(:,4*nCor+1:5*nCor), pTrans(:,5*nCor+1:6*nCor), pTrans(:,6*nCor+1:7*nCor), pTrans(:,7*nCor+1:8*nCor), pTrans(:,8*nCor+1:9*nCor)]; % transmural pressure
P.CorArtVen.AMyo = zeros(size(A(:,nCor+1:9*nCor)));
P.CorArtVen.AMyo(:,[iCorMyoArt,iCorMyoArt2P,iCorMyoArt2M,iCorMyoArt2N,iCorMyoVen2P,iCorMyoVen2M,iCorMyoVen2N,iCorMyoVen])     = [     A(:,nCor+1:2*nCor) ,     A(:,2*nCor+1:3*nCor),      A(:,3*nCor+1:4*nCor),      A(:,4*nCor+1:5*nCor),      A(:,5*nCor+1:6*nCor),      A(:,6*nCor+1:7*nCor),      A(:,7*nCor+1:8*nCor),      A(:,8*nCor+1:9*nCor)]; % cross-section

end
