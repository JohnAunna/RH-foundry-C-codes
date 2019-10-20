//Drag& lift &Virtual mass force and k & epsilon sources  in VOF designed by hamid    2014/09/11-evening												

#include "udf.h"
																											
real Random();
real con(cell_t c, Thread *t);
real Cd(cell_t c,Thread*t);
real Re(cell_t c, Thread *t);
real d_bubble(cell_t c,Thread*t);
real Ux,Uy,Uz,Vx,Vy,Vz,dVxdx,dVxdy,dVxdz,dVydx,dVydy,dVydz,dVzdx,dVzdy,dVzdz,dUzdz,dUydy,dUxdx,Pb,g,l,xL,yL,zL,xD,x_dS,yD,y_dS,zD,z_dS,k,ep,xV,yV,zV;
real Ck=0.1;
real C1=1.44;
real Ce=0.08;
real Cm=0.09;
real r_lift=0.1;
real r_drag=0.1;
real r_k=0.01;
real r_ep=0.01;
real r_virt=0.1;
//real r_omeg=0.1;


DEFINE_EXECUTE_AT_END(EXE_end)
{
	Domain *d = Get_Domain(1);
	cell_t c;
	Thread *t,**pt;
	real xyz[ND_ND],z,x,y;			
																																																																																						//real omeg;
	 thread_loop_c(t,d)																												/*loop for all of cells*/ 
	{
		begin_c_loop(c,t)
		{
			pt = THREAD_SUB_THREADS(t);
			C_CENTROID(xyz,c,t);
			z=xyz[2];
			y=xyz[1];																															//coordniates of cell in order to limiting the sources on meshing domain
			x=xyz[0];
			if((C_VOF(c,pt[0])>0.1)&&(C_VOF(c,pt[1])>0.1)&&z>2.9&&z<4.6&&x>0.355&&x<1.005&&y>-0.325&&y<0.325&&sqrt((x-0.68)*(x-0.68)+y*y)<0.325)			//limiting the sources to special part of the meshing domain        
{																																				/*velocity and velocity derivatives equations*/ 
			Ux = C_U(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0);					
			Uy = C_V(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0);					//velocity of gas in every cell	 based on ratio of gas phase in every cell+random velocity
			Uz = C_W(c,t)*C_VOF(c,pt[1]) +Random()*sqrt(2.0*C_K(c,t)/3.0);
				Vx = C_U(c,t)*C_VOF(c,pt[0]);
				Vy = C_V(c,t)*C_VOF(c,pt[0]);																			//velocity of liquid in every cell	based on ratio of liquid phase in every cell
				Vz = C_W(c,t)*C_VOF(c,pt[0]);
					dVxdx=C_DUDX(c,t)*C_VOF(c,pt[0]); 
					dVxdy=C_DUDY(c,t)*C_VOF(c,pt[0]);
					dVxdz=C_DUDZ(c,t)*C_VOF(c,pt[0]);
					dVydx=C_DVDX(c,t)*C_VOF(c,pt[0]);
					dVydy=C_DVDY(c,t)*C_VOF(c,pt[0]);																//velocity derivative of liquid in every cell	
					dVydz=C_DVDZ(c,t)*C_VOF(c,pt[0]);
					dVzdx=C_DWDX(c,t)*C_VOF(c,pt[0]);
					dVzdy=C_DWDY(c,t)*C_VOF(c,pt[0]);
					dVzdz=C_DWDZ(c,t)*C_VOF(c,pt[0]);
					dUxdx=C_DUDX(c,t)*C_VOF(c,pt[1]); 
					dUydy=C_DVDX(c,t)*C_VOF(c,pt[1]);															//velocity derivative of gas in every cell	
					dUzdz=C_DWDX(c,t)*C_VOF(c,pt[1]); 

xL=0.1*C_R(c,t)*((Uy-Vy)*(dVydx-dVxdy)-(Uz-Vz)*(dVxdz-dVzdx));							//lift_x_equation
yL=0.1*C_R(c,t)*((Uz-Vz)*(dVzdy-dVydz)-(Ux-Vx)*(dVydx-dVxdy));							//lift_y_equation
zL=0.1*C_R(c,t)*((Ux-Vx)*(dVxdz-dVzdx)-(Uy-Vy)*(dVzdy-dVydz));							//lift_z_equation
xD= con(c,t)*fabs(Ux-Vx)*(Ux-Vx) ;																																			/*drag_x_equation*/
x_dS=2.*con(c,t)*fabs(Ux-Vx) ;																																				/*drag_x_ds equation*/
yD=con(c,t)*fabs(Uy-Vy)*(Uy-Vy) ;																																			/*drag_y_equation*/
y_dS=2.*con(c,t)*fabs(Uy-Vy) ;																																				/*drag_y_ds equation*/
zD= con(c,t)*fabs(Uz-Vz)*(Uz-Vz) ;																																			/*drag_z_equation*/
z_dS=2.*con(c,t)*fabs(Uz-Vz) ;																																				/*drag_z_ds equation*/
Pb=(3./4.)*C_MU_L(c,t)*Cd(c,t)*Re(c,t)*(pow((Ux-Vx),2.0)+pow((Uy-Vy),2.0)+pow((Uz-Vz),2.0))/pow(d_bubble(c,t),2.0)/C_VOLUME(c,t);			/*turbulence Pb equation*/
k=Ck*Pb;																																													/*k equatoin*/
ep=Pb*C1*Ce*C_D(c,t)/C_K(c,t);																																			/*ep equatoin*/ 
																																																																																																//omeg=ep/k/Cm;																									/*omega equatoin*/ 
xV=0.5*C_R(c,t)*((Vx-C_UDMI(c,t,14))/CURRENT_TIMESTEP-(Ux-C_UDMI(c,t,15))/CURRENT_TIMESTEP+Vx*dVxdx-Ux*dUxdx);		/* virtual mass_x_ equation*/
yV=0.5*C_R(c,t)*((Vy-C_UDMI(c,t,16))/CURRENT_TIMESTEP-(Uy-C_UDMI(c,t,17))/CURRENT_TIMESTEP+Vy*dVydy-Uy*dUydy); 		/* virtual mass_y_ equation*/
zV=0.5*C_R(c,t)*((Vz-C_UDMI(c,t,18))/CURRENT_TIMESTEP-(Uz-C_UDMI(c,t,19))/CURRENT_TIMESTEP+Vz*dVzdz-Uz*dUzdz);		/* virtual mass_z_ equation*/
C_UDMI(c,t,14)=Vx;																										//from last time step data for Virtual mass force equation
C_UDMI(c,t,15)=Ux;																										//from last time step data for Virtual mass force equation
C_UDMI(c,t,16)=Vy;																										//from last time step data for Virtual mass force equation
C_UDMI(c,t,17)=Uy;																										//from last time step data for Virtual mass force equation
C_UDMI(c,t,18)=Vz;																										//from last time step data for Virtual mass force equation
C_UDMI(c,t,19)=Uz;																										//from last time step data for Virtual mass force equation




//C_UDMI(c,t,20)+=r_omeg*(omeg - C_UDMI(c,t,20));
//Message("VM():%g\n" ,((Vx-C_UDMI(c,t,14))/CURRENT_TIMESTEP-(Ux-C_UDMI(c,t,15))/CURRENT_TIMESTEP+Vx*dVxdx-Ux*dUxdx));
//Message("xD:%g yD:%g zD:%g\n",xD,yD,zD);
//Message("dSxD:%gdS yD: g dSzD:%g\n",x_dS,y_dS,z_dS);
//Message("xL:%g yL:%g zL:%g\n",xL,yL,zL);
//Message("x:%g y:%g z:%g\n",x,y,z);
//Message("xV:%g yV:%g zV:%g\n",xV,yV,zV);
//Message("Pb:%g k:%g epsilon:%g omeg: %g\n",Pb,k,ep,omeg);
//Message("Re:%gCd:%g con:%g Rand_num:%g\n",Re(c,t),Cd(c,t),con(c,t),Random());
//Message("Pb:%g C_MU_water:%g Cd:%g Rey:%g pow((Ux-Vx),2.0): %g alpha gas:%g  alpha water:%g db:%g  cvolume:%g\n",Pb,C_MU_L(c,pt[0]),Cd(c,t),Re(c,t),pow((Ux-Vx),2.0),C_VOF(c,pt[1]),C_VOF(c,pt[0]),d_bubble(c,t),C_VOLUME(c,t));
}
else
{
	xL=yL=zL=xD=yD=zD=x_dS=y_dS=z_dS=xV=yV=zV=k=ep=0;

			}
	C_UDMI(c,t,20)=xL;	
	C_UDMI(c,t,21)=yL;
	C_UDMI(c,t,22)=zL;
	C_UDMI(c,t,23)=xD;
	C_UDMI(c,t,24)=x_dS;
	C_UDMI(c,t,25)=yD;
	C_UDMI(c,t,26)=y_dS;
	C_UDMI(c,t,27)=zD;
	C_UDMI(c,t,28)=z_dS;
	C_UDMI(c,t,29)=k;
	C_UDMI(c,t,30)=ep;
	C_UDMI(c,t,31)=xV;
	C_UDMI(c,t,32)=yV;
	C_UDMI(c,t,33)=zV;
		}
		end_c_loop(c,t)
	}
}

DEFINE_ADJUST(AD_just,d)													/*under relaxation factor equations for all sources for every iteration*/ 		
{
	cell_t c;
	Thread *t;
	thread_loop_c(t,d)																							/*loop for all of cells*/ 
	{
		begin_c_loop(c,t)
		{
			
C_UDMI(c,t,0)+=r_lift*(C_UDMI(c,t,20)- C_UDMI(c,t,0)); 
C_UDMI(c,t,1)+=r_lift*(C_UDMI(c,t,21) - C_UDMI(c,t,1));
C_UDMI(c,t,2)+=r_lift*(C_UDMI(c,t,22) - C_UDMI(c,t,2));
C_UDMI(c,t,3)+=r_drag*(C_UDMI(c,t,23)- C_UDMI(c,t,3));
C_UDMI(c,t,4)+=r_drag*(C_UDMI(c,t,24) - C_UDMI(c,t,4));
C_UDMI(c,t,5)+=r_drag*(C_UDMI(c,t,25) - C_UDMI(c,t,5));
C_UDMI(c,t,6)+=r_drag*(C_UDMI(c,t,26) - C_UDMI(c,t,6));
C_UDMI(c,t,7)+=r_drag*(C_UDMI(c,t,27) - C_UDMI(c,t,7));
C_UDMI(c,t,8)+=r_drag*(C_UDMI(c,t,28) - C_UDMI(c,t,8));
C_UDMI(c,t,9)+=r_k*(C_UDMI(c,t,29)- C_UDMI(c,t,9));
C_UDMI(c,t,10)+=r_ep*(C_UDMI(c,t,30) - C_UDMI(c,t,10));
C_UDMI(c,t,11)+=r_virt*(C_UDMI(c,t,31)- C_UDMI(c,t,11));
C_UDMI(c,t,12)+=r_virt*(C_UDMI(c,t,32) - C_UDMI(c,t,12));
C_UDMI(c,t,13)+=r_virt*(C_UDMI(c,t,33) - C_UDMI(c,t,13));
		}end_c_loop(c,t)
	}
}




DEFINE_SOURCE(Lift_xmom,c,t,dS,eqn)   													   /*lift Xmom source*/       
{
		 return (C_UDMI(c,t,0));
}

DEFINE_SOURCE(Lift_ymom,c,t,dS,eqn)													   /*lift Ymom source*/  
{
	
	return (C_UDMI(c,t,1));
}

DEFINE_SOURCE(Lift_zmom,c,t,dS,eqn)													   /*lift Zmom source*/  
{
	return (C_UDMI(c,t,2));
}

DEFINE_SOURCE(Drag_xmom,c,t,dS,eqn)												   /*Drag XZmom source*/  
{
		dS[eqn]=C_UDMI(c,t,4);
	return C_UDMI(c,t,3);
}

DEFINE_SOURCE(Drag_ymom,c,t,dS,eqn)												    /*Drag Ymom source*/  
{
		dS[eqn]=C_UDMI(c,t,6);
	return C_UDMI(c,t,5);
}
DEFINE_SOURCE(Drag_zmom,c,t,dS,eqn)													  /*Drag Zmom source*/  
{
		dS[eqn]=C_UDMI(c,t,8);
		return C_UDMI(c,t,7);
}

DEFINE_SOURCE(Virtual_mass_xmom,c,t,dS,eqn)   										/*Virtual mass Xmom source*/       
{
		 return C_UDMI(c,t,11);
}

DEFINE_SOURCE(Virtual_mass_ymom,c,t,dS,eqn)											  /*Virtual mass Ymom source*/  
{
	
		return C_UDMI(c,t,12);
}

DEFINE_SOURCE(Virtual_mass_zmom,c,t,dS,eqn)											/*Virtual mass Zmom source*/  
{
		return C_UDMI(c,t,13);
}

DEFINE_SOURCE(kinetic_turb_energy,c,t,dS,eqn)											  /* k source*/ 
{
	return C_UDMI(c,t,9);
}
	
DEFINE_SOURCE(epsilon,c,t,dS,eqn)																	    /* epsilon source*/ 
{
	return C_UDMI(c,t,10);
}

//DEFINE_SOURCE(omega,c,t,dS,eqn)															    /* omega  source*/ 
//{
//	return C_UDMI(c,t,20);
//}

																							//Functions:

real Re(cell_t c,Thread *t)																									/*Reynolds number function for every cell*/ 
{
	Thread **pt = THREAD_SUB_THREADS(t);
			real Ugx,Ugy,Ugz,Vlx,Vly,Vlz;
			Ugx = C_U(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0);						
			Ugy = C_V(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0);
			Ugz = C_W(c,t)*C_VOF(c,pt[1]) +Random()*sqrt(2.0*C_K(c,t)/3.0);
				Vlx = C_U(c,t)*C_VOF(c,pt[0]);
				Vly = C_V(c,t)*C_VOF(c,pt[0]);
				Vlz = C_W(c,t)*C_VOF(c,pt[0]);
	return (sqrt(pow((Ugx-Vlx),2.) +pow((Ugy-Vly),2.)+pow((Ugz-Vlz),2.))*d_bubble(c,t)*C_R(c,pt[0])/C_MU_T(c,pt[0]));		
}

real Cd(cell_t c,Thread*t)																										/*Cd function for Drag force fo different particle shape*/ 
 {
real Cdcap,Cdellipse,Cdsphere,Cddist,Cdrag, Rey_num;
Thread **pt=THREAD_SUB_THREADS(t);
Rey_num=Re(c,t);
Cdcap=8./3.;																																	/* cap particle*/ 
Cdellipse=(2./3.)*sqrt(9.81*(C_R(c,pt[0])-C_R(c,pt[1]))*d_bubble(c,t)*d_bubble(c,t)/1.2);			/* ellipse particle*/ 
Cdsphere=24.0*(1.0+0.15*pow(Rey_num,0.687))/Rey_num;														/*sphere particle*/ 
if (Cdellipse<Cdcap)
Cddist=Cdellipse;
else Cddist=Cdcap;
if(Cdsphere>Cddist)
Cdrag=Cdsphere;
else Cdrag=Cddist;	
return Cdrag;
}

real con(cell_t c,Thread*t)																								/*Constant number function for Drag force*/ 
 {
real Alg;
Thread **pt=THREAD_SUB_THREADS(t);
Alg=C_VOF(c,pt[0])*C_VOF(c,pt[1])*C_VOLUME(c,t)/d_bubble(c,t);									//CFX theory guide interphase drag 
return (Cd(c,t)*C_R(c,t)*Alg);																			
 }
	
real d_bubble(cell_t c,Thread*t)																						/*equivalent diameter of bubble in every cell */ 
 {
Thread **pt=THREAD_SUB_THREADS(t);
return pow(6.0*C_VOF(c,pt[1])*C_VOLUME(c,t)/3.141592653589793238463,(1./3.));																			
 }
	



real Random()												/*random numbr function for Gaussian distribution*/ 
{
	real random_number,u,v;
	srand((unsigned)time(NULL));	
	random_number = rand();
	u = random_number*(1./(float)RAND_MAX);
	random_number = rand();
	v = random_number*(1./(float)RAND_MAX);
	return ((pow((-2.*log(u)),0.5)*cos(2*M_PI*v))/5.);
}

DEFINE_ON_DEMAND(monitor_x_sources)		//for monitoring the x-sources and k & epsilon, this macro would produce a text file which has the value of  sources
{
Domain *d = Get_Domain(1);
	cell_t c;
	Thread *t,**pt;
	real xyz[ND_ND],z,x,y;			
	FILE *tpL;
	FILE *tpL_detail;
	FILE *tpD;
	FILE *tpD_detail;
	FILE *tpdS;
	FILE *tpV;
	FILE *tpV_detail;
	FILE *tpk;
	FILE *tpk_detail;
	FILE *tpep;
	FILE *tpep_detail;
	FILE *tptime;
	tpL = fopen("F:\\output-L.txt","a+");
	tpL_detail = fopen("F:\\output-L_detail.txt","a+");
	tpD = fopen("F:\\output-D.txt","a+");
	tpD_detail = fopen("F:\\output-D_detail.txt","a+");
	tpdS = fopen("F:\\output-dS.txt","a+");
	tpV = fopen("F:\\output-VM.txt","a+");
	tpV_detail = fopen("F:\\output-VM_detail.txt","a+");
	tpk = fopen("F:\\output-k.txt","a+");
	tpk_detail = fopen("F:\\output-k_detail.txt","a+");
	tpep = fopen("F:\\output-ep.txt","a+");
	tpep_detail = fopen("F:\\output-ep_detail.txt","a+");
	tptime=fopen("F:\\output-time.txt","a+");
	 thread_loop_c(t,d)																												//loop for all of cells
	{
		begin_c_loop(c,t)
		{
			pt = THREAD_SUB_THREADS(t);
			C_CENTROID(xyz,c,t);
			z=xyz[2];
			y=xyz[1];
			x=xyz[0];
			if((C_VOF(c,pt[0])>0.1)&&(C_VOF(c,pt[1])>0.1)&&z>2.9&&z<4.6&&x>0.355&&x<1.005&&y>-0.325&&y<0.325&&sqrt((x-0.68)*(x-0.68)+y*y)<0.325)			//limiting the sources to special part of the meshing        
			{
			fprintf(tpL,"L:%g \n",C_UDMI(c,t,0));
			fprintf(tpL_detail,"X:%g  Y:%g  Z:%g  Density%g  u%g  v%g  %g  w%g\n",x,y,z,C_R(c,t),C_U(c,t),C_V(c,t),C_W(c,t),con(c,t));							//((C_V(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0))-(C_V(c,t)*C_VOF(c,pt[0]))),((C_DVDX(c,t)*C_VOF(c,pt[0]))-(C_DUDY(c,t)*C_VOF(c,pt[0]))),((C_W(c,t)*C_VOF(c,pt[1]) +Random()*sqrt(2.0*C_K(c,t)/3.0))-(C_W(c,t)*C_VOF(c,pt[0]))),((C_DUDZ(c,t)*C_VOF(c,pt[0]))-(C_DWDX(c,t)*C_VOF(c,pt[0]))));			//PRINTING ALL OF THE VALUE IN LIFT FORCE
			fprintf(tpD,"D:%g \n",C_UDMI(c,t,3));
			fprintf(tpD_detail,"X:%g  Y:%g  Z:%g  Cd%g  u%g  v%g  w%g\n",x,y,z,Cd(c,t),C_U(c,t),C_V(c,t),C_W(c,t),con(c,t));							 //fabs((C_U(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0))-(C_U(c,t)*C_VOF(c,pt[0]))),((C_U(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0))-(C_U(c,t)*C_VOF(c,pt[0]))));				//PRINTING ALL OF THE VALUE IN DRAG FORCE
			fprintf(tpdS,"dS:%g \n",C_UDMI(c,t,4));
			fprintf(tpV,"VM:%g \n",C_UDMI(c,t,11));
			fprintf(tpV_detail,"X:%g  Y:%g  Z:%g  Den:%g  u:%g  v:%g  w:%g  con:%g\n",x,y,z,C_R(c,t),C_U(c,t),C_V(c,t),C_W(c,t),con(c,t));									//(((C_U(c,t)*C_VOF(c,pt[0]))-C_UDMI(c,t,14))/CURRENT_TIMESTEP-((C_U(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0))-C_UDMI(c,t,15))/CURRENT_TIMESTEP+(C_U(c,t)*C_VOF(c,pt[0]))*C_DUDX(c,t)*C_VOF(c,pt[0])-C_U(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0)*C_DUDX(c,t)*C_VOF(c,pt[1])));		//PRINTING ALL OF THE VALUE IN VIRTUAL MASS FORCE
			fprintf(tpk,"%g \n",C_UDMI(c,t,9));
			fprintf(tpk_detail,"X:%g  Y:%g  Z:%g  Visc%g  Cd:%g  Re:%g  k:%g  dbub:%g  cvolum:%g   con:%g\n",x,y,z,C_MU_L(c,t),Cd(c,t),Re(c,t),C_K(c,t),d_bubble(c,t),C_VOLUME(c,t),con(c,t));		//(pow(((C_U(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0))-(C_U(c,t)*C_VOF(c,pt[0]))),2.0)+pow(((C_V(c,t)*C_VOF(c,pt[1]) + Random()*sqrt(2.0*C_K(c,t)/3.0))-(C_V(c,t)*C_VOF(c,pt[0]))),2.0)+pow(((C_W(c,t)*C_VOF(c,pt[1]) +Random()*sqrt(2.0*C_K(c,t)/3.0))-(C_W(c,t)*C_VOF(c,pt[0]))),2.0)),d_bubble(c,t),C_VOLUME(c,t));		//PRINTING ALL OF THE VALUE IN k SOURCE
			fprintf(tpep,"%g \n",C_UDMI(c,t,10));
			fprintf(tpep_detail,"X:%g  Y:%g  Z:%g  den:%g  k:%g  con:%g \n",x,y,z,C_D(c,t),C_K(c,t),con(c,t));	//PRINTING ALL OF THE VALUE IN EPSILON SOURCE
			fprintf(tptime,"timestep%g \n",CURRENT_TIMESTEP);
			}
		}
		end_c_loop(c,t) 
	}
	
}
