fC = (1-(C/Cmax)) * pow((1-C),2);

rho = C * rhoif + ( 1 - C ) * rhoif;

etaN = 0.75 * pow(C/Cmax,2) * pow(1-C/Cmax,-2);
etaS = pow(1-C/Cmax,-2);
etaP = etaS - 1.;

lambda2 = 0.81 * (C/Cmax) + 0.66;
lambda3 = -0.0088 * (C/Cmax) + 0.54;

gradU = fvc::grad(U);
E = gradU + gradU.T();
tensor I(1.,0,0,0,1.,0,0,0,1.); 
Q = (1.+lambda2 + lambda3)/3.*I;
gammaShear = sqrt(2*(E && E));

sigmaP = 2*eta0*etaP*E - eta0*etaN*gammaShear*Q;
sigmaF = 2.*eta0*E;

if (abs(rhoif.value()-rhoip.value())<0.01) 
{
    jflux =  (2./9.*(pow(a,2)*fC)/(eta0)*fvc::div(sigmaP));
}
else
{
    jflux = fC * C * vst*100.;   //+ 1000.*(2./9.*(pow(a,2)*fC)/(eta0)*fvc::div(sigmaP));
}

// volVectorField jg((fC*C*vst));
//volVectorField j((2./9.*(pow(a,2)*fC)/(eta0)*fvc::div(sigmaP)));