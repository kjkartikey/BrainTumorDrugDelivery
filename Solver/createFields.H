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


/*Info << "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);*/


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


/*Info<< "Reading field J\n" << endl;
    volScalarField J
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

  Info<< "Reading field I\n" << endl;
    volScalarField I
    (
        IOobject
        (
            "I",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

   Info<< "Reading field L\n" << endl;
    volScalarField L
    (
        IOobject
        (
            "L",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );*/

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


Info << "Reading field D\n" << endl;
volScalarField D
(
    IOobject
    (
        "D",
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

/*Info << "Reading field Z\n" << endl;
volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);*/

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

IOdictionary cpDict
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

dimensionedScalar Cp("Cp", dimensionSet(0, -3, 0, 0, 1, 0, 0), 0.0);


dimensionedScalar ktransavg
(
    transportProperties.lookup("ktransavg")
);

dimensionedScalar Daa
(
    transportProperties.lookup("Daa")
);

dimensionedScalar kaae
(
    transportProperties.lookup("kaae")
);

dimensionedScalar kak
(
    transportProperties.lookup("kak")
);

dimensionedScalar alpha
(
    transportProperties.lookup("alpha")
);

dimensionedScalar beta
(
    transportProperties.lookup("beta")
);

dimensionedScalar gamma
(
    transportProperties.lookup("gamma")
);
#   include "createPhi.H" 

