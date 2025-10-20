/*---------------------------------------------------------------------------*\
  createFields.H
  Robust initialization of fields for converging cylindrical shock
  using Friedlander profile and stiffened gas EOS.
  Automatically reads gamma, pInf (a), and Cv from phaseProperties.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <cmath>

// ------------------------------------------------------------------------- //
// Read/create fields

Info<< "Reading/creating fields\n" << endl;

volScalarField p
(
    IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
);

volScalarField rho1
(
    IOobject("rho1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
);

volScalarField rho2
(
    IOobject("rho2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
);

volScalarField alpha1
(
    IOobject("alpha1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
);

volScalarField T
(
    IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
);

volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
);

Info<< "Fields ready.\n" << endl;

// ------------------------------------------------------------------------- //
// Read stiffened-gas parameters from phaseProperties

IOdictionary phaseProps
(
    IOobject("phaseProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
);

Info<< "Reading fluid (liquid) phase parameters from phaseProperties\n" << endl;

// --- Read liquid (fluid) phase properties
const dictionary& fluidDict = phaseProps.subDict("fluid");

scalar gamma = 4.4;
scalar pInf  = 6.0e8;
scalar C_heat = 4285;

if (fluidDict.found("equationOfState"))
{
    const dictionary& eos = fluidDict.subDict("equationOfState");
    gamma = eos.lookupOrDefault<scalar>("gamma", gamma);
    pInf  = eos.lookupOrDefault<scalar>("a", pInf);
}
if (fluidDict.found("thermodynamics"))
{
    const dictionary& thermo = fluidDict.subDict("thermodynamics");
    C_heat = thermo.lookupOrDefault<scalar>("Cv", C_heat);
}

Info<< "  gamma   = " << gamma << nl
    << "  pInf    = " << pInf  << nl
    << "  Cv      = " << C_heat << nl << endl;

// ------------------------------------------------------------------------- //
// Shock and ambient parameters

scalar r0      = 2.29e-03;
scalar a0      = 0.1e-03;
int n          = 16;

scalar Ps      = 815e6;
scalar t_star  = 0.9e-6;
scalar u_shock = 2200;

scalar cs      = 1500;
scalar rho_amb = 1000;

scalar P0      = 1.01e5;        // ambient pressure
scalar P2      = Ps;            // shock peak

scalar r_arr   = 15.5e-03;
scalar r_shock = r_arr - cs*t_star;

scalar x_off = 0;
scalar y_off = 0;

const scalar SMALL    = 1e-12;
const scalar U_CAP    = 1.0e5;
const scalar pAmbient = 1.01e5;

// ------------------------------------------------------------------------- //
// Correct stiffened-gas post-shock density

scalar rho_post = rho_amb *
    ((gamma+1)*(P2 + pInf) + (gamma-1)*(P0 + pInf)) /
    ((gamma-1)*(P2 + pInf) + (gamma+1)*(P0 + pInf));

// Post-shock particle velocity
scalar U_post = Foam::sqrt( (P2 - P0) * (1.0/rho_amb - 1.0/rho_post) );

Info<< "Post-shock rho: " << rho_post << ", U_post: " << U_post << endl;

// ------------------------------------------------------------------------- //
// Initialize fields with debug protections

label count_rprime  = 0;
label count_rmag    = 0;
label count_pfried  = 0;
label count_rhoadbt = 0;

forAll(p, cellI)
{
    scalar xx = mesh.C()[cellI].x();
    scalar yy = mesh.C()[cellI].y();

    // Interface perturbation
    scalar azim  = Foam::atan2(yy, xx);
    scalar pertR = r0 - a0*Foam::cos(n*azim);
    scalar r2    = xx*xx + yy*yy;
    alpha1[cellI] = (r2 <= pertR*pertR) ? 0.0 : 1.0;

    // Position relative to shock center
    scalar x_shift = xx - x_off;
    scalar y_shift = yy - y_off;
    scalar r_mag   = Foam::sqrt(x_shift*x_shift + y_shift*y_shift);

    // --- Converging shock: r_prime positive inside shock front ---
    scalar r_prime = r_mag - r_shock;      // positive for r < r_shock
    scalar tEff    = r_prime / Foam::max(u_shock, SMALL);

    // Friedlander pressure
    scalar p_friedlander = pAmbient;
    if (tEff >= 0.0)
    {
        scalar arg       = tEff / Foam::max(t_star, SMALL);
        scalar friedTerm = (1.0 - arg) * Foam::exp(-arg);
        p_friedlander    = Ps * friedTerm;
    }

    // Debug and clamp
    if (r_prime < 1e-8 && count_rprime < 10)
    {
        Info << "Warning: r_prime very small at cell " << cellI << ": " << r_prime << endl;
        count_rprime++;
    }
    if (r_mag < 1e-8 && count_rmag < 10)
    {
        Info << "Warning: r_mag very small at cell " << cellI << ": " << r_mag << endl;
        count_rmag++;
    }
    if (p_friedlander < pAmbient && count_pfried < 10)
    {
        Info << "Warning: p_friedlander below ambient at cell " << cellI << ": " << p_friedlander << endl;
        count_pfried++;
        p_friedlander = pAmbient; // clamp
    }

    if (p_friedlander > pAmbient + SMALL)
    {
        // Post-shock density from adiabatic compression
        scalar ratio    = Foam::max(p_friedlander / Ps, 1e-12);
        scalar rho_adbt = rho_post * Foam::pow(ratio, 1.0/gamma);

        if (rho_adbt < SMALL && count_rhoadbt < 10)
        {
            Info << "Warning: rho_adbt very small at cell " << cellI << ": " << rho_adbt << endl;
            count_rhoadbt++;
            rho_adbt = SMALL; // clamp
        }

        // Safe radii
        scalar safe_rprime = Foam::max(Foam::mag(r_prime), 1e-6);
        scalar safe_rmag   = Foam::max(Foam::mag(r_mag), 1e-6);

        // Velocity magnitude
        scalar U_mag = U_post * (rho_amb / rho_adbt) * (r_shock / safe_rprime);

        // Smooth near focus
        if (r_mag < 0.5*r_shock)
            U_mag *= (r_mag / (0.5*r_shock));

        // Cap
        if (U_mag > U_CAP)
            U_mag = U_CAP;

        // Internal energy and temperature
        scalar e_adbt = (p_friedlander + pInf) / ((gamma-1.0)*rho_adbt);
        scalar T_adbt = e_adbt / C_heat;

        // Assign fields
        p[cellI]    = p_friedlander;
        rho1[cellI] = rho_adbt;

        scalar ux = x_shift / safe_rmag;
        scalar uy = y_shift / safe_rmag;
        U[cellI].x() = -U_mag * ux;
        U[cellI].y() = -U_mag * uy;

        T[cellI] = T_adbt;
    }
    else
    {
        p[cellI]    = pAmbient;
        rho1[cellI] = rho_amb;
        U[cellI]    = vector::zero;
        T[cellI]    = 298.0;
    }

    rho2[cellI] = 1.293; // air density
}

// ------------------------------------------------------------------------- //
// Write fields

U.write();
p.write();
alpha1.write();
rho1.write();
rho2.write();
T.write();

Info<< "Field initialization complete.\n" << endl;

/*---------------------------------------------------------------------------*\
                         End of createFields.H
\*---------------------------------------------------------------------------*/
