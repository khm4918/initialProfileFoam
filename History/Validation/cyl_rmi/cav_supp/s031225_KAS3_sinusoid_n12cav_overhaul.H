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

scalar r0      = 2.3e-03;
scalar a0      = 0.1e-03;
int n          = 12;

scalar Ps      = 380e6;
scalar t_star  = 0.9e-6;
scalar u_shock = 1900;

scalar cs      = 1500;
scalar rho_amb = 1030;			// according to J-boy

scalar P0      = 1.01e5;        // ambient pressure
scalar P2      = Ps;            // shock peak

scalar r_arr   = 15.5e-03;
scalar r_shock = r_arr - u_shock*t_star;

scalar x_off = 0;
scalar y_off = 0;

const scalar SMALL    = 1e-12;
const scalar U_CAP    = 1.0e5;
const scalar pAmbient = 1.01e5;

// Cavity parameters

int ncav = 12;
scalar r_cav=9.3e-03;

float x_p0={-8e-03};
float y_p0={0};
	
float x_p1={-9.1e-03};			//averaged from x,y positions of visible cavity vertices in shot
float y_p1={1.7e-03};
	
float x_p2={x_p1};
float y_p2={-1*y_p1};

// ------------------------------------------------------------------------- //

// Cavity defining equations and rotation info

std::vector<std::vector<float>> lines = {
					{(y_p1-y_p0)/(x_p1-x_p0),y_p1-(y_p1-y_p0)/(x_p1-x_p0)*x_p1}, //L0 - p0 to p1
					{(y_p0-y_p2)/(x_p0-x_p2),y_p0-(y_p0-y_p2)/(x_p0-x_p2)*x_p0}  //L1 - p0 to p2
					}; 

//rotation matrix params
scalar rot=M_PI*2/(ncav);

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

    // Interface perturbation and cavity zero
    scalar azim  = Foam::atan2(yy, xx);
    scalar pertR = r0 - a0*Foam::cos(n*azim);
    scalar r2    = xx*xx + yy*yy;
    alpha1[cellI] = ( (r2 <= pertR*pertR) || ((r2 <= r_cav*r_cav) and (yy<=lines[0][0]*xx+lines[0][1]) and (yy>=lines[1][0]*xx+lines[1][1])) ) ? 0.0 : 1.0;

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
        scalar U_mag = U_post; //* (rho_amb / rho_adbt) * (r_shock / safe_rprime);

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
		T[cellI] = T_adbt;


        scalar ux = -x_shift / safe_rmag;
        scalar uy = -y_shift / safe_rmag;
		U[cellI] = vector(U_mag*ux, U_mag*uy, 0.0);
	}
    else
    {
        p[cellI]    = pAmbient;
		rho1[cellI] = rho_amb;
		T[cellI]    = 298.0;
		U[cellI]    = vector::zero;
    }
    
    rho2[cellI] = 1.293; // air density
}

//create vector to store x, y coords of cavity zero
	List<vector> targetCoords;	

 //find cavity zero cell coords
          forAll(p,cellI){
    scalar xc=mesh.C()[cellI].x();
    scalar yc=mesh.C()[cellI].y();
	
					//go to radii outside of central air cylinder
                    if(
					
					((xc*xc + yc*yc) > (r0 + a0)*(r0 + a0))
					
					//find air cell
					and (alpha1[cellI]==0)
					
					)
					
						{
							//Info << "Found target cell " << endl;
							
							//calculate corresponding rotation cell coords 
							for (int n=1; n<int(ncav); n++) 
								
							{
										scalar xrot = xc*Foam::cos(rot*n)-yc*Foam::sin(rot*n);
										scalar yrot = xc*Foam::sin(rot*n)+yc*Foam::cos(rot*n);
										targetCoords.append(vector(xrot, yrot, 0.0));
	
							}
						
						}

                     }

//now target all cells nearest to identified coords and set alpha1 to 0
	forAll(targetCoords, i)
	{
		vector target = targetCoords[i];
		label nearest = -1;
		scalar minDist = 1e-6;

		forAll(mesh.C(), cellI)
		{
			scalar d = magSqr(mesh.C()[cellI] - target);
			if (d < minDist)
			{
				minDist = d;
				nearest = cellI;
			}
		}

		if (nearest != -1)
		{
			//Info << "Updating alpha at nearest cell " << nearest
				//<< " (coords " << mesh.C()[nearest] << ")" << endl;
			alpha1[nearest] = 0.0;
		}
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
