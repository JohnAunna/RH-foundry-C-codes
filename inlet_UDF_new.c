/***************************************
		UDF for inlet velocity
***************************************/
#include "udf.h"
cell_t nozzle_c[12];
Thread *nozzle_t[12];
real nozzle_coor[12][3] = {-0.741,-0.32194,3.588,-0.4651897,-0.15577,3.588,-0.4711897,0.1661684,3.588,-0.753,0.3219446,3.588,-1.0284,0.1557761,3.588,-1.02281,-0.166168,3.588,-0.58083,-0.27581,3.713,-0.425055,0.006,3.713,-0.5912238,0.28181,3.713,-0.913168,0.27581,3.713,-1.06894,-0.006,3.713,-0.9027761,-0.28181,3.713};//  input coordinates of inlets
real nozzle_vel[12][3];
real Q[12]={0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,0.0022336,}; //gas flow rate (m3/s) for every inlet/ 
DEFINE_ON_DEMAND(Find_coor)
{
	Domain *dd = Get_Domain(1);
	Thread *t,*tf;
	cell_t c;
	face_t f;
	Node *bfn[4],*pn;
	real d[12] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
	real nx[ND_ND],dx[ND_ND],bfx[ND_ND],A_derec[ND_ND],A_v[ND_ND],distance,A;
	int j,n,m,k,point,point_vec;
	thread_loop_c(t,dd) //find the cell and thread of inlet cell for the 12 nozzles
	{
		begin_c_loop(c,t)
		{
			C_CENTROID(nx,c,t);
			for (j=0;j<12;j++)
			{
				NV_VV(dx, =, nozzle_coor[j], -, nx);
				distance = NV_MAG(dx);
				if (distance < d[j])
				{
					d[j] = distance;
					nozzle_c[j] = c;
					nozzle_t[j] = t;
				}
			}
		}
		end_c_loop(c,t)
	}
	thread_loop_c(t,dd)//calcualte the velocity for the inlet cell
	{
		begin_c_loop(c,t)
		{
			C_CENTROID(nx,c,t);
			for (j=0;j<12;j++)
			{
				if ((nozzle_c[j] == c)&&(nozzle_t[j] == t))
				{
					c_face_loop(c,t,n)
					{
						if (BOUNDARY_FACE_THREAD_P(t))
						{
							f = C_FACE(c, t, n);
							tf = C_FACE_THREAD(c, t, n);
							f_node_loop(f,tf,m)
								bfn[m] = F_NODE(f,tf,m);
						}
					}
					c_face_loop(c,t,n)
					{
						point = 1;
						f = C_FACE(c,t,n);
						tf = C_FACE_THREAD(c,t,n);
						f_node_loop(f,tf,m)
						{
							pn = F_NODE(f,tf,m);
							for (k = 0; k < 4; k++)
							{
								if (bfn[k] == pn)
								{
									point = 0;
									break;
								}
							}
							if (point == 1)
							{
								F_AREA(A_v,f,tf);
								F_CENTROID(bfx,f,tf);
								NV_VV(A_derec, =, bfx, -,nx);
								point_vec = 1;
								for (k = 0; k < 3; k++)
								{
									if ((A_v[k]*A_derec[k]) < 0)
										point_vec = 0;
								}
								A = NV_MAG2(A_v);
								if (point_vec)
								{
									for (k = 0; k < 3; k++)
										nozzle_vel[j][k] = Q[j]*A_v[k]/A;
								}
								else
								{
									for (k = 0; k < 3; k++)
										nozzle_vel[j][k] = Q[j]*A_v[k]*(-1.)/A;
								}
							}
						}
					}
				}
			}
		}
		end_c_loop(c,t)
	}
}
DEFINE_EXECUTE_AT_END(Tapping)  //defining velocity and Vod of inlet cell
{
	Domain *d = Get_Domain(1);
	Thread *t,**pt;
	cell_t c;
	real xc[ND_ND];
	int j;
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
			for (j = 0; j<12; j++)
			{
				if ((nozzle_c[j] == c)&&(nozzle_t[j] == t))
				{
					pt = THREAD_SUB_THREADS(t);
					C_VOF(c,pt[0]) = 0.;//
					C_VOF(c,pt[1]) = 1.;//
					C_U(c,t) = nozzle_vel[j][0];
					C_V(c,t) = nozzle_vel[j][1];
					C_W(c,t) = nozzle_vel[j][2];
				}
			}
		}
		end_c_loop(c,t)
	}
}