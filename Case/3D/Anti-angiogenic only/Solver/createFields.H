Info << "Reading transportProperties\n" << endl;

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

Info << "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Cf\n" << endl;
volScalarField Cf
(
    IOobject
    (
        "Cf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Caa\n" << endl;
volScalarField Caa
(
    IOobject
    (
        "Caa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field pvf\n" << endl;
volScalarField pvf
(
    IOobject
    (
        "pvf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field por\n" << endl;
volScalarField por
(
    IOobject
    (
        "por",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field ktrans\n" << endl;
volScalarField ktrans
(
    IOobject
    (
        "ktrans",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field hc\n" << endl;
volScalarField hc
(
    IOobject
    (
        "hc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field sapuv\n" << endl;
volScalarField sapuv
(
    IOobject
    (
        "sapuv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field lfc\n" << endl;
volScalarField lfc
(
    IOobject
    (
        "lfc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field lp\n" << endl;
volScalarField lp
(
    IOobject
    (
        "lp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Peff\n" << endl;
volScalarField Peff
(
    IOobject
    (
        "Peff",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Daa\n" << endl;
volScalarField Daa
(
    IOobject
    (
        "Daa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field beta\n" << endl;
volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field gamma\n" << endl;
volScalarField gamma
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Vcell\n" << endl;
volScalarField Vcell
(
    IOobject
    (
        "Vcell",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Vtissue\n" << endl;
volScalarField Vtissue
(
    IOobject
    (
        "Vtissue",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field infSite\n" << endl;
volScalarField infSite
(
    IOobject
    (
        "infSite",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field w\n" << endl;
volScalarField w
(
    IOobject
    (
        "w",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field kb\n" << endl;
volScalarField kb
(
    IOobject
    (
        "kb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Ustar\n" << endl;
volVectorField Ustar
(
    IOobject
    (
        "Ustar",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field Dstar\n" << endl;
volScalarField Dstar
(
    IOobject
    (
        "Dstar",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field kstar\n" << endl;
volScalarField kstar
(
    IOobject
    (
        "kstar",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/*IOdictionary cpDict
(
    IOobject
    (
        "cpDict",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);



autoPtr<Function1<scalar> > cpValue(Function1<scalar>::New("cpValue", cpDict));

dimensionedScalar Cp("Cp", dimensionSet(0, -3, 0, 0, 1, 0, 0), 0.0);*/


dimensionedScalar ktransavg
(
    transportProperties.lookup("ktransavg")
);

dimensionedScalar kaae
(
    transportProperties.lookup("kaae")
);

dimensionedScalar kak
(
    transportProperties.lookup("kak")
);

dimensionedScalar Q
(
    transportProperties.lookup("Q")
);

dimensionedScalar C_infSite
(
    transportProperties.lookup("C_infSite")
);

dimensionedScalar Vinf
(
    transportProperties.lookup("Vinf")
);

dimensionedScalar Df
(
    transportProperties.lookup("Df")
);

dimensionedScalar ke
(
    transportProperties.lookup("ke")
);

dimensionedScalar Cinf
(
    transportProperties.lookup("Cinf")
);

dimensionedScalar Pics_ecs
(
    transportProperties.lookup("Pics_ecs")
);

dimensionedScalar Pcm_ecs
(
    transportProperties.lookup("Pcm_ecs")
);

dimensionedScalar Kecs
(
    transportProperties.lookup("Kecs")
);

dimensionedScalar Kics
(
    transportProperties.lookup("Kics")
);

dimensionedScalar P
(
    transportProperties.lookup("P")
);

#   include "createPhi.H" 

