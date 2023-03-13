Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoip("rhoip",transportProperties);
dimensionedScalar rhoif("rhoif",transportProperties);
dimensionedScalar eta0("eta0",transportProperties);
dimensionedScalar Cmax("Cmax",transportProperties);
dimensionedVector g("g",transportProperties);
dimensionedScalar a("a",transportProperties);
dimensionedScalar vstdt("vstdt",transportProperties);


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<<"Reading field C\n"<<endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField fC
(
    IOobject
    (
        "fC",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ((1-(C/Cmax)) * pow((1-C),2))
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rhoi
(
    IOobject
    (
        "rhoi",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(rhoip - rhoif)
);


volVectorField vst
(
    IOobject
    (
        "vst",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ((2./9.*a*a*rhoi/eta0)*g)
);


Info<< "Reading field jflux\n" << endl;
volVectorField jflux
(
    IOobject
    (
        "jflux",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (C * rhoif + ( 1. - C ) * rhoif)
);

volScalarField etaS
(
    IOobject
    (
        "etaS",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (pow(1.-C/Cmax,-2))
);

volScalarField etaN
(
    IOobject
    (
        "etaN",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (0.75 * pow(C/Cmax,2) * pow(1.-C/Cmax,-2))
);

volScalarField etaP
(
    IOobject
    (
        "etaP",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (etaS - 1.)
);

volScalarField lambda2
(
    IOobject
    (
        "lambda2",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (0.81 * (C/Cmax) + 0.66)
);

volScalarField lambda3
(
    IOobject
    (
        "lambda3",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (-0.0088 * (C/Cmax) + 0.54)
);

volTensorField gradU
(
    IOobject
    (
        "gradU",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (fvc::grad(U))
);

volTensorField E
(
    IOobject
    (
        "E",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (gradU + gradU.T())
);

volScalarField gammaShear
(
    IOobject
    (
        "gammaShear",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (sqrt(2*(E && E)))
);

volTensorField Q
(
    IOobject
    (
        "Q",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (1.+lambda2 + lambda3)/3.*tensor(1.,0,0,0,1.,0,0,0,1.)
);

volTensorField sigmaF
(
    IOobject
    (
        "sigmaF",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (2.*eta0*E)
);

volTensorField sigmaP
(
    IOobject
    (
        "sigmaP",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (2*eta0*etaP*E - eta0*etaN*gammaShear*Q)
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

