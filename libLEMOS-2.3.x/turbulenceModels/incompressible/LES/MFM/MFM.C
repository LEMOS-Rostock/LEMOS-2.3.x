#include "MFM.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
 {
namespace incompressible
 {
namespace LESModels
 {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MFM, 0);
addToRunTimeSelectionTable(LESModel, MFM, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MFM::MFM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
  LESModel(modelName, U, phi, transport, turbulenceModelName),
  thermoDict
       (
           IOobject
           (
               "thermophysicalProperties",
               U.mesh().time().constant(),
               U.mesh(),
               IOobject::MUST_READ,
               IOobject::NO_WRITE
           )
       ),
  filterPtr(LESfilter::New(U.mesh(), coeffDict(), "filter")),
  filter(filterPtr()),
  filter2Ptr(LESfilter::New(U.mesh(), coeffDict(), "filter2")),
  filter2(filter2Ptr()),
  viscLengthScale_
  (
     IOobject
     (
      "viscLengthScale",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      0.0,
      "zeroGradient"
  ),
  Ureynolds_
  (
     IOobject
     (
      "UsgsUsgs",
      runTime_.timeName(),
      mesh_,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_
      //dimensionedTensor("UU", dimensionSet(0,2,-2,0,0,0,0), tensor(0,0,0,0,0,0,0,0,0)),
      //"zeroGradient"
  ),
  uDelta_
  (
     IOobject
     (
      "uDelta",
      runTime_.timeName(),
      mesh_,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_
   ),
  sgsDissipation_
  (
     IOobject
     (
      "sgsDissipation",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedScalar("UU", dimensionSet(0,2,-3,0,0,0,0), 0.0),
      "zeroGradient"
  ),
  N_
  (
     IOobject
     (
      "N",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedScalar("N", dimless, 10.0),
      zeroGradientFvPatchScalarField::typeName
  ),

 Re_
  (
     IOobject
     (
      "Re",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedScalar("Re", dimless, 10.0),
      zeroGradientFvPatchScalarField::typeName
  ),
  Cai_
  (
     IOobject
     (
      "Cai",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedScalar("Cai_", dimless, 0.0),
      zeroGradientFvPatchScalarField::typeName
   )
{
  updateSubGridScaleFields(dev(symm(fvc::grad(U)))); 

  printCoeffs();
  
  Info << "MFM created" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void MFM::updateSubGridScaleFields(const volSymmTensorField& S)
{
   updateCai(S);
   updateN(S);
}


// Reynolds number based on the norm of strain rate tensor
void MFM::ReS(const volSymmTensorField& S) 
{
   Re_ = mag(S)*sqr(delta())/nu();
   max(Re_,Re_,1.0);
   Info << "updating Reynolds number ReS: " << average(Re_) << endl;

   Re_.correctBoundaryConditions();
}

// Reynolds number based on the local velocity
void MFM::ReU(const volSymmTensorField& S)
{
   Re_ = mag(U())*(delta())/nu();
   max(Re_,Re_,1.0);
   Info << "updating Reynolds number ReU: " << average(Re_) << endl;
   Re_.correctBoundaryConditions();
}
       

// Anisotropy factor, based on Reynolds number 
void MFM::updateCai(const volSymmTensorField& S) 
{
   ReS(S);
   //ReU(S);
   Cai_= 1.0-pow(Re_,-3.0/16.0);
   max(Cai_,Cai_,0.0);
 
   Info << "updating anisotropy factor: " << average(Cai_) << endl;
   Cai_.correctBoundaryConditions();
}


scalar MFM::Csgs() const
{
   return readScalar(thermoDict.lookup("Csgs"));
}

tmp<volSymmTensorField> MFM::B() const
{
  return  symm(-(F1()*U()*uDelta_+ F1()*uDelta_*U()) - (sqr(F1())*uDelta_*uDelta_)) ;
}

// Scaling factor 2^(-2/3*N)*sqrt(2^(4/3*N)-1) 
tmp<volScalarField> MFM::F1() const
{
  //return  Csgs()*Cai_*pow((1-pow(alpha,-4.0/3.0)),-0.5)*pow(2.0,-2.0/3.0*N_)*sqrt(pow(2.0,4.0/3.0*N_)-1.0);
  return  Csgs()*pow((1-pow(alpha,-4.0/3.0)),-0.5)*pow(2.0,-2.0/3.0*N_)*sqrt(pow(2.0,4.0/3.0*N_)-1.0);
}

tmp<volSymmTensorField> MFM::devReff() const
{
  return  - 2*nu()*dev(symm(fvc::grad(U())));
}



tmp<volScalarField> MFM::epsilon() const
{
  return
  (
     -(F1()*U()*uDelta_+ F1()*uDelta_*U() + sqr(F1())*uDelta_*uDelta_)
     &&
     dev(symm(fvc::grad(U())))
  );
}


tmp<volScalarField> MFM::k() const
{
   return tmp<volScalarField>
   (
      new volScalarField
      (
         IOobject
         (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
         mesh_,
         dimensionedScalar("k", nu()().dimensions(), 0.0)
      )
   );
}



tmp<fvVectorMatrix> MFM::divDevReff(volVectorField& U) const
{
  volVectorField uDelta("uDelta", U-filter2(U));

  return
  (
     - fvm::laplacian(nu(), U, "laplacian(nu,U)")
     + filter( fvc::div((F1()*U*uDelta)+ (F1()*uDelta*U) + (sqr(F1())*uDelta*uDelta), "div(MFM)")) 
  );
}

tmp<fvVectorMatrix> MFM::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

/*
tmp<volVectorField> MFM::F(const volScalarField &f) const
{

  return
    (
      *LeoPhi_[f.name()] - (*turbulentDiffusivity_[f.name()]) * fvc::grad(f)
     );

}

*/
tmp<volScalarField> MFM::molecularDiffusivityCoeff(word name) const
{
  return
    (
     laminarDiffusivity_[name] 
     );
}

tmp<volScalarField> MFM::turbulentDiffusivityCoeff(word name) const
{
  return
    (
     *turbulentDiffusivity_[name]
     );
}

tmp<volVectorField> MFM::Feff(const volScalarField &f) const
{
  return
    (
      *LeoPhi_[f.name()] - ( laminarDiffusivity_[f.name()] + *turbulentDiffusivity_[f.name()] ) * fvc::grad(f)
     );
}


tmp<fvScalarMatrix> MFM::divFeff(volScalarField &f) const
{
  return
     (
      fvm::laplacian((laminarDiffusivity_[f.name()] + *turbulentDiffusivity_[f.name()]), f, "laplacian(Deff,F)")
      - fvc::div((*LeoPhi_[f.name()]))
      );
}

tmp<fvScalarMatrix> MFM::divFsgs(volScalarField &f) const
{
  return
     (
      fvm::laplacian((*turbulentDiffusivity_[f.name()]), f, "laplacian(Deff,F)")
      - fvc::div((*LeoPhi_[f.name()]))
      );
}


void MFM::correct(const tmp<volTensorField>& gradU)
{
  LESModel::correct(gradU);
  uDelta_ = U()-filter2(U());
  sgsDissipation_ = epsilon();
  viscLengthScale_ = F1();
  Ureynolds_ = B();

  /// strain rate tensor
  updateSubGridScaleFields(dev(symm(gradU)));
  //Info << "average B: " << average(F1()) << endl;
}

// Number of cascade steps log2(delta/Kolmogorov)
void MFM::updateN(const volSymmTensorField& S)
{

  Info << "calculating number of cascade steps" << endl;

//  if(viscLengthScale_.headerOk())
  //{
  //  // N_ = (log(delta()/viscLengthScale_)/log(2.0));
 // }
 // else
  {
       ReU(S);
     // ReS(S);

       scalar cv = readScalar(thermoDict.lookup("cv"));


       scalarField tmp = cv*pow(Re_,0.75);
       max(tmp,tmp,1.0);

       N_.internalField() = log(tmp)/log(2.0);
       N_.correctBoundaryConditions();
       //Info << max(N_) << endl;
       //Info << min(N_) << endl;

     Info << "updating number of cascade steps: " << average(N_) << endl;
  }
/*
   volScalarField cellDist
   (
                    IOobject
                    (
                        "N",
                        U().time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    0.0,
                    zeroGradientFvPatchScalarField::typeName
                );
  cellDist.internalField()=Cai_*Csgs();
  cellDist.write();
*/
}


void MFM::registerScalarField(volScalarField &f, scalar molecularDiffusivityCoeff) 
{

	word name = f.name();

	Info << "register ScalarField " <<  name << endl;

	registeredScalarFields_.insert
	(
		name,
		f
	);

        laminarDiffusivity_.insert
        (
          name,
          volScalarField(
          	IOobject
                (
                        "D_"+name,
                        f.time().timeName(),
                        f.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                ),
                f.mesh(),
                dimensionedScalar("D_"+name,dimensionSet(0,2,-1,0,0,0,0), molecularDiffusivityCoeff)
          )
        );

        turbulentDiffusivity_.insert
        (
          name,
          new volScalarField (
                IOobject
                (
                        "Dt_"+name,
                        f.time().timeName(),
                        f.mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                ),
                f.mesh(),
                dimensionedScalar("Dt_"+name,dimensionSet(0,2,-1,0,0,0,0), 0.0)
          )
        );


        LeoPhi_.insert
        (
                name,
                new volVectorField
                (
                        IOobject
                        (
                                "leonardScalarFlux_" + name,
                                runTime_.timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedVector("zero",  dimensionSet(0,1,-1,0,0,0,0) * f.dimensions() , pTraits<vector>::zero)
                )
        );



}

//Abschluss
bool MFM::read()
{
    if (LESModel::read())
    {
        filter.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 } // End namespace LESmodels
 }
 } // End namespace Foam

// ************************************************************************* //
