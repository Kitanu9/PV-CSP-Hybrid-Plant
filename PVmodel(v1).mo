model FirstPV 
  parameter Integer c = 25;
  parameter Real pi = Modelica.Constants.pi;
  parameter Real d = pi/180; 
  parameter Real S = 1367;
  parameter Real phi = -20*d;
  parameter Real day = 60;
  parameter Real beta = 45*d; 
  parameter Real e= Modelica.Constants.e;
  
  parameter Real Voc = 41.3;
  parameter Real Isc = 11.37;
  parameter Real T = 25;
  parameter Real k = Modelica.Constants.k;
  parameter Real q = Modelica.Constants.q; 
  parameter Real a = -3.47;
  parameter Real b = -0.0594;
  parameter Real Ta = 25;
  parameter Real WS = 10;
  parameter Real Tcnd = 3;
  parameter Real TCV = 0.0023;
  parameter Real TCI = 0.0005;
  
  Real hourtime[c] = 6:(12/(c-1)):18;
  Real delta;
  Real omega[c];
  Real costz[c];
  Real theta[c];
  Real Ib[c];
  Real Ib_beta[c];
  Real Id[c];
  Real Id_beta[c];
  Real AM[c];
  Real alpha[c];
  
  Real Voc_X[c];
  Real Isc_X[c];
  Real X[c];
  Real Gpoa[c];
  Real Tm[c];
  Real Tcell[c];
  Real VocTcell[c];
  Real IscTcell[c];
   
equation
  delta = 23.45*sin((360/365)*(284+day));
  omega = d.*(15).*(hourtime .- 12);
  costz = sin(phi)*sin(delta) .+ cos(phi)*cos(delta).*cos(omega);
  theta = acos(costz); 
  alpha = (pi/2) .- theta;
  AM = (1)./cos(theta);
  Ib = S.*(0.7.^(AM.^0.678));
  Ib_beta = Ib.*sin(beta .+ alpha);
  Id = Ib.*0.1;
  Id_beta = Id.*(1 + cos(beta))./(1 .+ sin(alpha));
  Gpoa = Ib_beta + Id_beta;
  X = Gpoa./1000;
  Voc_X = Voc .+ (((k*T)/q).*log(X));
  Isc_X = Isc.*X;
  Tm = Gpoa .* exp(a+b*WS) .+ Ta;
  Tcell = Tm .+ X.*Tcnd;
  VocTcell = Voc_X.*(1 .- TCV.*(Tcell .- 25));
  IscTcell = Isc_X.*(1 .+ TCI.*(Tcell .- 25));
  
end FirstPV;
