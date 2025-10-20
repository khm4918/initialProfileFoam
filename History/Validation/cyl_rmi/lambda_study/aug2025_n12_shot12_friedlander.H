 
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

//fluid boundary geometry params
	scalar r0=2.29e-03;
	scalar a0=0.1e-03;
	int n = 16;

//friedlander shock params
	scalar Ps = 815e06;
	scalar t_star = 0.9e-06;

//measured shock params
	scalar u_shock = 2200;

//shock thermo params
	scalar cs = 1500; 			//sound speed in water
	scalar rho_amb = 1000;		//ambient density of water
	scalar gamma = 4.4;			//gamma from stiffened gas EoS
	scalar pInf = 6.0e8;   		//or appropriate for water stiffened gas	
	scalar C_heat = 4285;		//heat capacity	
	
	scalar P0 = 1.01e5 + pInf;
	scalar P2 = Ps + pInf;
	
//shock geometry params
	scalar r_arr=15.5e-03;						//initial position of wire array. currently used as 0 pressure point post shock
	scalar r_shock = r_arr - (cs * t_star);		//shock front position
	
	//offset from centre
	scalar x_off=-6*(32e-6);
	scalar y_off=-8*(32e-6);

 //low order initialization
          forAll(p,cellI){
    scalar xx=mesh.C()[cellI].x();
    scalar yy=mesh.C()[cellI].y();					
	
					//radial perturbation to fluid-gas boundary
                    if((xx*xx+yy*yy)<=(r0-a0*Foam::cos(n*Foam::atan(yy/xx)))*(r0-a0*Foam::cos(n*Foam::atan(yy/xx)))
					
					)
						
					{    //Foam::cos()

                       alpha1[cellI]=0;
					  
					}else{
						
						alpha1[cellI]=1.0; 
					
					//shock annulus in water
                    }    //Foam::cos()
						
					//shift to shock centre
					scalar x_shift = xx-x_off;
					scalar y_shift = yy-y_off;
					
					//shift out to higher radius
					scalar r_prime = Foam::sqrt(x_shift*x_shift + y_shift*y_shift) - r_shock;
					
					//friedlander equation for pressure
					scalar p_friedlander = Ps * (1 - (r_prime/(u_shock*t_star))) * Foam::exp(-1*(r_prime/(u_shock*t_star)));
					
					//immediate post shock quantities according to R-H conditions
					scalar rho_post = rho_amb * ((gamma+1)*(u_shock*u_shock))/(2*(cs*cs)+(gamma-1)*(u_shock*u_shock));
					scalar U_post = 2*(u_shock-(cs*cs)/(u_shock))/(gamma+1);
					
					//set initial fields for shock region assuming adiabatic shock interior
					if((p_friedlander <= Ps)and(p_friedlander >= 0)){
						
						//calculate derived post shock thermo quantities
						scalar rho_adbt = rho_post*Foam::pow(p_friedlander/Ps, 1.0/gamma);
						scalar U_mag = U_post*(rho_amb/rho_adbt)*(r_shock/r_prime);
						scalar e_edbt = p_friedlander/((gamma-1)*rho_adbt);
						
						//set thermo quantities in field
						p[cellI]=p_friedlander;						
						rho1[cellI]= rho_adbt;
						U[cellI].x()=-((U_mag)/Foam::sqrt(x_shift*x_shift+y_shift*y_shift))*x_shift;
						U[cellI].y()=-((U_mag)/Foam::sqrt(x_shift*x_shift+y_shift*y_shift))*y_shift;
						T[cellI]=e_edbt/C_heat;

                     }else{
                                           
						p[cellI]=1.01e5;
						rho1[cellI]=1000;
						U[cellI].x()=0.0;
						U[cellI].y()=0.0;
						T[cellI]=298.0; 					   
                       
					   }
        
                      rho2[cellI]=1.293;
       
                     }
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
