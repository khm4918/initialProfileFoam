 
 /*   dimensionedScalar R
    (
        transportProperties.lookup("R")
    );

    dimensionedScalar Cv
    (
        transportProperties.lookup("Cv")
    );*/

    Info<< "Reading field pressure\n" << endl;
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
    

    Info<< "Reading field rho\n" << endl;
    volScalarField rho1
    (
        IOobject
        (
            "rho1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
        volScalarField rho2
    (
        IOobject
        (
            "rho2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField alpha1    
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field T\n" << endl;
  volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
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

 Info<< "Finish initialization field" << endl;

// ---------- helper: point-in-polygon (ray casting) ----------
auto pointInPolygon = [](const std::vector<vector>& poly, const scalar x, const scalar y)
{
    bool inside = false;
    size_t nvert = poly.size();
    for (size_t i = 0, j = nvert - 1; i < nvert; j = i++)
    {
        scalar xi = poly[i].x();
        scalar yi = poly[i].y();
        scalar xj = poly[j].x();
        scalar yj = poly[j].y();

        bool intersect = ((yi > y) != (yj > y)) &&
                         (x < (xj - xi) * (y - yi) / (yj - yi + SMALL) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
};

// ---------- define geometry / params ----------
scalar r0 = 2.29e-03;
scalar a0 = 0.1e-03;
int n = 16;           // used for radial perturbation
int ncav = 16;        // number of cavities to create by rotation
scalar theta_lower = 0.68722;		//angular range of inserted defect
scalar theta_upper = 0.88357;

scalar r_in  = 15.25e-03;
scalar r_out = 15.75e-03;
scalar r_cav = 10.00e-03;

// rotation amount per cavity
scalar rot = M_PI*2 / scalar(ncav);

// ---------- build polygon for cavity zero (in CCW or CW order) ----------
// Use OpenFOAM vector for points; adjust vertices as needed.
// Example uses the points you gave for the base cavity (one polygon)
std::vector<vector> cavity0;
cavity0.push_back(vector(0.0, 9e-03, 0.0));
cavity0.push_back(vector(0.39e-03, 9.17e-03, 0.0));
cavity0.push_back(vector(0.84e-03, 9.35e-03, 0.0));
cavity0.push_back(vector(1.42e-03, 9.72e-03, 0.0));
cavity0.push_back(vector(-1.42e-03, 9.72e-03, 0.0));
cavity0.push_back(vector(-0.84e-03, 9.35e-03, 0.0));
cavity0.push_back(vector(-0.39e-03, 9.17e-03, 0.0));
// If your polygon is closed, you don't need to duplicate first point.

// ---------- INITIALIZE alpha1 to 1 (water) first ----------
forAll(alpha1, cellI)
{
    alpha1[cellI] = 1.0;
}

// ---------- main cell loop: set alpha1=0 for central air and rotated cavities ----------
forAll(p, cellI)
{
    scalar xx = mesh.C()[cellI].x();
    scalar yy = mesh.C()[cellI].y();

    // central air cylinder with radial perturbation & defect (use atan2)
    scalar ang = Foam::atan2(yy, xx);
    if (ang < 0) ang += M_PI*2; // normalize to [0, 2pi) if you prefer

    bool isCentralAir = false;
    {
        // radial perturbation
        scalar rCell = Foam::sqrt(xx*xx + yy*yy);
        scalar perturbedR = r0 - a0 * Foam::cos(n * Foam::atan2(yy, xx));
        if ((ang > theta_lower && ang < theta_upper) &&
            (rCell <= perturbedR * perturbedR ? rCell*rCell <= perturbedR*perturbedR : false))
        {
            isCentralAir = true;
        }
        // remaining circular region (outside defect sector)
        if (!isCentralAir && ((ang < theta_lower || ang > theta_upper) && (xx*xx + yy*yy) <= (r0*r0)))
        {
            isCentralAir = true;
        }
        // mirrored/in-place defect rule from your original code (left-bottom region) - keep or remove as intended
        if (!isCentralAir && (ang > theta_lower && ang < theta_upper) && (xx < 0) && (yy < 0) && ((xx*xx + yy*yy) <= (r0*r0)))
        {
            isCentralAir = true;
        }
    }

    if (isCentralAir)
    {
        alpha1[cellI] = 0.0;
    }
    else
    {
        // Now check the rotated cavities:
        // For each cavity index k, rotate the cell coordinate *back* by +rot*k
        // (equivalently rotate the cell by -rot*k and test in the base polygon).
		bool cellIsInAnyCavity = false;
		scalar rCell2 = xx*xx + yy*yy;

		for (int k = 0; k < ncav && !cellIsInAnyCavity; ++k)
		{
			scalar angleK = rot * scalar(k);

			scalar c = Foam::cos(angleK);
			scalar s = Foam::sin(angleK);
			scalar x_local = c*xx + s*yy;
			scalar y_local = -s*xx + c*yy;

			if (pointInPolygon(cavity0, x_local, y_local) && (rCell2 <= r_cav*r_cav))
			{
				cellIsInAnyCavity = true;
			}
		}

		
        if (cellIsInAnyCavity)
        {
            alpha1[cellI] = 0.0;
        }
        else
        {
            alpha1[cellI] = 1.0;
        }
    }

    // --- now the rest of your low-order initialization (shock annulus, p, rho1, U, etc.) ---
    if ((r_in*r_in <= (xx*xx + yy*yy)) && ((xx*xx + yy*yy) <= r_out*r_out))
    {
        p[cellI] = 8.15e08;
        rho1[cellI] = 1200;
        scalar rr = Foam::sqrt(xx*xx + yy*yy) + SMALL;
        U[cellI].x() = -(370.8/rr) * xx;
        U[cellI].y() = -(370.8/rr) * yy;
    }
    else
    {
        p[cellI] = 1.01e5;
        rho1[cellI] = 1000;
        U[cellI].x() = 0.0;
        U[cellI].y() = 0.0;
    }

    T[cellI]  = 298.0;
    rho2[cellI] = 1.293;
} // end forAll(p,cellI)


      // high order initialization    
      
        /* label np=40;
                    
         forAll(U,cellI)
        {
            scalar sumY0=0.0;
            scalar sumY1=0.0;
            scalar sumY2=0.0;
            //vector sumV=vector::zero;
            scalar sumNum=0.0;
            if(mixMeshType[cellI]==0)
            {
                label pt[4];
                pt[0]=mixVertLabels[cellI][0];
                pt[1]=mixVertLabels[cellI][1];
                pt[2]=mixVertLabels[cellI][2];
                pt[3]=mixVertLabels[cellI][3];

                scalar xM[4],yM[4];
                for(label m=0;m<4;m++)
                {
                    xM[m]=mesh.points()[pt[m]].x(); 
                    yM[m]=mesh.points()[pt[m]].y();
                }

                scalar N[4];
                for(int i=0; i<=np; i++)
                {
                    for(int j=0; j<=np; j++)
                    {
                        scalar wgI = -1.0+double(i)*2.0/double(np);
                        scalar wgJ = -1.0+double(j)*2.0/double(np);
                        N[0]=(1.0-wgI)*(1.0-wgJ)/4.0;
                        N[1]=(1.0+wgI)*(1.0-wgJ)/4.0;
                        N[2]=(1.0+wgI)*(1.0+wgJ)/4.0;
                        N[3]=(1.0-wgI)*(1.0+wgJ)/4.0;    
                        scalar xx = N[0]*xM[0]+N[1]*xM[1]+N[2]*xM[2]+N[3]*xM[3];
                        scalar yy = N[0]*yM[0]+N[1]*yM[1]+N[2]*yM[2]+N[3]*yM[3];
                        if(((xx-0.5)*(xx-0.5)+(yy-0.5)*(yy-0.5))<=0.0625)
                        {
                            sumY0+=1.0e-6;
                            sumY1+=1.0e-6;
                            sumY2+=0.999998;
                        }
                        else
                        {
                            sumY0+=0.999998;
                            sumY1+=1.0e-6;
                            sumY2+=1.0e-6;
                        }
                        sumNum+=1.0;
                    }
                }
            }
            else
            {
                label pt[3];
                pt[0]=mixVertLabels[cellI][0];
                pt[1]=mixVertLabels[cellI][1];
                pt[2]=mixVertLabels[cellI][2];

                scalar xM[3],yM[3];
                for(label m=0;m<3;m++)
                {
                   xM[m]=mesh.points()[pt[m]].x(); 
                   yM[m]=mesh.points()[pt[m]].y();
                } 

                scalar N[3];

           }

             Y[0][cellI] = sumY0/sumNum;
             Y[1][cellI] = sumY1/sumNum;
             Y[2][cellI] = sumY2/sumNum;
            p[cellI]=1.01e5;
            T[cellI]=298.0;  
           // rho[cellI] = sumR/sumNum;
        }
        */



    U.write();
    p.write();
    alpha1.write();
    rho1.write();
    rho2.write();
    T.write();
    U.write();
    

