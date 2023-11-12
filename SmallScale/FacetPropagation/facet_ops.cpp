
void gen_subfacet(facet &afacet, long double nxdivs, long double nydivs, int xth, int yth, facet &asubfacet)
	{
	long double kx1, kx2, ky1,ky2;
	int i;
	long double temp, dxy, coef, modfactor;
	
	kx1=xth/nxdivs;
	kx2=(xth+1)/nxdivs;
	ky1=yth/nydivs;
	ky2=(yth+1)/nydivs;
	asubfacet.x0=afacet.x0*(1.-kx1)+kx1*afacet.x1;
	asubfacet.x1=afacet.x0*(1.-kx2)+kx2*afacet.x1;
	asubfacet.y0=afacet.y0*(1.-ky1)+ky1*afacet.y1;
	asubfacet.y1=afacet.y0*(1.-ky2)+ky2*afacet.y1;
	asubfacet.centre[0]=(asubfacet.x0+asubfacet.x1)/2.;
	asubfacet.centre[1]=(asubfacet.y0+asubfacet.y1)/2.;
	
	// recall here z01 is x0y1, but z[1] is x1y0
	asubfacet.z00=afacet.z00*(1-kx1)*(1.-ky1) + afacet.z01*(1-kx1)*ky1 + afacet.z10*kx1*(1.-ky1) + afacet.z11*kx1*ky1;
	asubfacet.z10=afacet.z00*(1-kx2)*(1.-ky1) + afacet.z01*(1-kx2)*ky1 + afacet.z10*kx2*(1.-ky1) + afacet.z11*kx2*ky1;
	asubfacet.z01=afacet.z00*(1-kx1)*(1.-ky2) + afacet.z01*(1-kx1)*ky2 + afacet.z10*kx1*(1.-ky2) + afacet.z11*kx1*ky2;
	asubfacet.z11=afacet.z00*(1-kx2)*(1.-ky2) + afacet.z01*(1-kx2)*ky2 + afacet.z10*kx2*(1.-ky2) + afacet.z11*kx2*ky2;
	asubfacet.centre[2]=(asubfacet.z00+asubfacet.z10+asubfacet.z01+asubfacet.z11)/4.;
	fill_facet_points(asubfacet);
	get_facet_normal(asubfacet);
	
	asubfacet.xproj=pow(1.+ asubfacet.norm[0] * asubfacet.norm[0]/( asubfacet.norm[2] * asubfacet.norm[2]), 0.5l);
	asubfacet.yproj=pow(1.+ asubfacet.norm[1] * asubfacet.norm[1]/( asubfacet.norm[2] * asubfacet.norm[2]), 0.5l);
	modfactor=pow(asubfacet.norm[2]*(asubfacet.xproj*asubfacet.yproj),0.5l);
	asubfacet.xproj/=modfactor;
	asubfacet.yproj/=modfactor;
	
	asubfacet.dx=(asubfacet.x1-asubfacet.x0)*asubfacet.xproj;
	asubfacet.dy=(asubfacet.y1-asubfacet.y0)*asubfacet.yproj;
			
	asubfacet.area=asubfacet.dx*asubfacet.dy;
	asubfacet.sqrtA=pow(asubfacet.area,0.5l);
	
	if (USE_WHICH==1)
		{
		gen_vec3_norm_to_vec2_from_vec1(XBASIS, asubfacet.norm, asubfacet.xbasis);
		cross(asubfacet.norm, asubfacet.xbasis, asubfacet.ybasis);
		}
	else if (USE_WHICH == 2)
		{
		gen_vec3_norm_to_vec2_from_vec1(YBASIS, asubfacet.norm, asubfacet.ybasis);
		cross(asubfacet.ybasis, asubfacet.norm, asubfacet.xbasis);
		}
	else
		{
		for (i=0; i<3; i++)
			{
			asubfacet.xbasis[i]=XBASIS[i]-asubfacet.norm[i]*asubfacet.norm[0];
			asubfacet.ybasis[i]=YBASIS[i]-asubfacet.norm[i]*asubfacet.norm[1];
			}
		normalise(asubfacet.xbasis);
		normalise(asubfacet.ybasis);
		dxy=dot(asubfacet.xbasis,asubfacet.ybasis);
		if (dxy == 0.)
			{
			coef = 0.;
			}
		else
			{
			coef = (-1.+pow(1.-dxy,0.5l))/dxy; // all div by 2
			}
		for (i=0; i<3; i++)
			{
			temp = asubfacet.xbasis[i]+coef*asubfacet.ybasis[i];
			asubfacet.ybasis[i] = asubfacet.ybasis[i]+coef*asubfacet.xbasis[i];
			asubfacet.xbasis[i] = temp;
			if (isnan(asubfacet.xbasis[i]))
				{
				cout<<"fuck it!";
				}
			}
		normalise(asubfacet.xbasis);
		normalise(asubfacet.ybasis);
		}
	}

void fill_facet_points(facet &afacet)
	{
	afacet.p[0][0]=afacet.x0;
	afacet.p[0][1]=afacet.y0;
	afacet.p[0][2]=afacet.z00;
	afacet.p[1][0]=afacet.x1;
	afacet.p[1][1]=afacet.y0;
	afacet.p[1][2]=afacet.z10;
	afacet.p[2][0]=afacet.x0;
	afacet.p[2][1]=afacet.y1;
	afacet.p[2][2]=afacet.z01;
	afacet.p[3][0]=afacet.x1;
	afacet.p[3][1]=afacet.y1;
	afacet.p[3][2]=afacet.z11;
	}



// x is slow, y is fast

// x is slow, y is fast
void init_facet(float *surface, int nx, int ny, int i, int j)
	{
	int k;
	long double coef, dxy, temp, modfactor;
	NFX=1;
	NFY=1;
	NFACETS=1;
	if (i==0 && j == 0) {FACETS = new facet*; FACETS[0] = new facet;}
	
	FACETS[0][0].x0=i*DD;
	FACETS[0][0].x1=(i+1)*DD;
	FACETS[0][0].y0=j*DD;
	FACETS[0][0].y1=(j+1)*DD;
			
	FACETS[0][0].z00=surface[i*ny+j];
	FACETS[0][0].z10=surface[(i+1)*ny+j];
	FACETS[0][0].z01=surface[i*ny+j+1];
	FACETS[0][0].z11=surface[(i+1)*ny+j+1];
			
	FACETS[0][0].p[0][0]=i*DD;
	FACETS[0][0].p[0][1]=j*DD;
	FACETS[0][0].p[0][2]=surface[i*ny+j];
	FACETS[0][0].p[1][0]=(i+1)*DD;
	FACETS[0][0].p[1][1]=j*DD;
	FACETS[0][0].p[1][2]=surface[(i+1)*ny+j];
	FACETS[0][0].p[2][0]=i*DD;
	FACETS[0][0].p[2][1]=(j+1)*DD;
	FACETS[0][0].p[2][2]=surface[i*ny+j+1];
	FACETS[0][0].p[3][0]=(i+1)*DD;
	FACETS[0][0].p[3][1]=(j+1)*DD;
	FACETS[0][0].p[3][2]=surface[(i+1)*ny+j+1];
			
	FACETS[0][0].centre[0]=(i+0.5)*DD;
	FACETS[0][0].centre[1]=(j+0.5)*DD;
	FACETS[0][0].centre[2]=0.25*(surface[(i+1)*ny+j+1] + surface[(i)*ny+j+1] + surface[(i+1)*ny+j] + surface[i*ny+j]);
			
	get_facet_normal(FACETS[0][0]);
	FACETS[0][0].area=DD*DD/FACETS[0][0].norm[2];
			
			
			// the facet dx, dy, and area is the real area, not the x-y plane projected area. It is always equal to or greater than [DD/DD^2 as appropriate].
	FACETS[0][0].xproj=pow(1.+ FACETS[0][0].norm[0] * FACETS[0][0].norm[0]/( FACETS[0][0].norm[2] * FACETS[0][0].norm[2]), 0.5l);
	FACETS[0][0].yproj=pow(1.+ FACETS[0][0].norm[1] * FACETS[0][0].norm[1]/( FACETS[0][0].norm[2] * FACETS[0][0].norm[2]), 0.5l);
	modfactor=pow(FACETS[0][0].norm[2]*(FACETS[0][0].xproj*FACETS[0][0].yproj),0.5l);
	FACETS[0][0].xproj/=modfactor;
	FACETS[0][0].yproj/=modfactor;
	
	FACETS[0][0].dx=DD*FACETS[0][0].xproj;
	FACETS[0][0].dy=DD*FACETS[0][0].yproj;
	
	FACETS[0][0].sqrtA=pow(FACETS[0][0].area,0.5l);
			
	if (USE_WHICH==1)
		{
		gen_vec3_norm_to_vec2_from_vec1(XBASIS, FACETS[0][0].norm, FACETS[0][0].xbasis);
		cross(FACETS[0][0].norm, FACETS[0][0].xbasis, FACETS[0][0].ybasis);
		}
	else if (USE_WHICH == 2)
		{
		gen_vec3_norm_to_vec2_from_vec1(YBASIS, FACETS[0][0].norm, FACETS[0][0].ybasis);
		cross(FACETS[0][0].ybasis, FACETS[0][0].norm, FACETS[0][0].xbasis);
		}
	else
		{
		// creates projections of x-y vectors onto the plane defined by .norm
		for (k=0; k<3; k++)
			{
			FACETS[0][0].xbasis[k]=XBASIS[k]-FACETS[0][0].norm[k]*FACETS[0][0].norm[0];
			FACETS[0][0].ybasis[k]=YBASIS[k]-FACETS[0][0].norm[k]*FACETS[0][0].norm[1];
			}
		normalise(FACETS[0][0].xbasis);
		normalise(FACETS[0][0].ybasis);
		dxy=dot(FACETS[0][0].xbasis,FACETS[0][0].ybasis);
		if (dxy == 0.)
			{
			coef = 0.;
			}
		else
			{
			coef = (-1.+pow(1.-dxy,0.5l))/dxy; // all div by 2
			}
		for (k=0; k<3; k++)
			{
			temp = FACETS[0][0].xbasis[k]+coef*FACETS[0][0].ybasis[k];
			FACETS[0][0].ybasis[k] = FACETS[0][0].ybasis[k]+coef*FACETS[0][0].xbasis[k];
			FACETS[0][0].xbasis[k] = temp;
			if (isnan(FACETS[0][0].xbasis[k]))
				{
				cout<<"fuck it!";
				}
			}
		normalise(FACETS[0][0].xbasis);
		normalise(FACETS[0][0].ybasis);
		}
	}

int init_facets(float *surface, int nx, int ny)
	{
	int i,j,k, OK=1;
	long double coef, dxy, temp, modfactor;
	NFX=nx-1;
	NFY=ny-1;
	NFACETS=NFX*NFY;
	cout<<"The memory required for facet initialisation is "<<(sizeof(facet)*NFX*NFY/1024/1024)<<" MB"<<endl;
	FACETS = new facet*[NFX];
	
	
	
	if (FACETS == 0) {OK=0; return OK;}
	for (i=0; i<NFX; i++)
		{
		FACETS[i] = new facet[NFY]; if (FACETS[i]==0) {OK=0; return OK;}
		for (j=0; j<NFY; j++)
			{
			FACETS[i][j].x0=i*DD;
			FACETS[i][j].x1=(i+1)*DD;
			FACETS[i][j].y0=j*DD;
			FACETS[i][j].y1=(j+1)*DD;
			
			FACETS[i][j].z00=surface[i*ny+j];
			FACETS[i][j].z10=surface[(i+1)*ny+j];
			FACETS[i][j].z01=surface[i*ny+j+1];
			FACETS[i][j].z11=surface[(i+1)*ny+j+1];
			
			FACETS[i][j].p[0][0]=i*DD;
			FACETS[i][j].p[0][1]=j*DD;
			FACETS[i][j].p[0][2]=surface[i*ny+j];
			FACETS[i][j].p[1][0]=(i+1)*DD;
			FACETS[i][j].p[1][1]=j*DD;
			FACETS[i][j].p[1][2]=surface[(i+1)*ny+j];
			FACETS[i][j].p[2][0]=i*DD;
			FACETS[i][j].p[2][1]=(j+1)*DD;
			FACETS[i][j].p[2][2]=surface[i*ny+j+1];
			FACETS[i][j].p[3][0]=(i+1)*DD;
			FACETS[i][j].p[3][1]=(j+1)*DD;
			FACETS[i][j].p[3][2]=surface[(i+1)*ny+j+1];
			
			FACETS[i][j].centre[0]=(i+0.5)*DD;
			FACETS[i][j].centre[1]=(j+0.5)*DD;
			FACETS[i][j].centre[2]=0.25*(surface[(i+1)*ny+j+1] + surface[(i)*ny+j+1] + surface[(i+1)*ny+j] + surface[i*ny+j]);
			
			get_facet_normal(FACETS[i][j]);
			FACETS[i][j].area=DD*DD/FACETS[i][j].norm[2];
			
			
			// the facet dx, dy, and area is the real area, not the x-y plane projected area. It is always equal to or greater than [DD/DD^2 as appropriate].
			FACETS[i][j].xproj=pow(1.+ FACETS[i][j].norm[0] * FACETS[i][j].norm[0]/( FACETS[i][j].norm[2] * FACETS[i][j].norm[2]), 0.5l);
			FACETS[i][j].yproj=pow(1.+ FACETS[i][j].norm[1] * FACETS[i][j].norm[1]/( FACETS[i][j].norm[2] * FACETS[i][j].norm[2]), 0.5l);
			modfactor=pow(FACETS[i][j].norm[2]*(FACETS[i][j].xproj*FACETS[i][j].yproj),0.5l);
			FACETS[i][j].xproj/=modfactor;
			FACETS[i][j].yproj/=modfactor;
			
			FACETS[i][j].dx=DD*FACETS[i][j].xproj;
			FACETS[i][j].dy=DD*FACETS[i][j].yproj;
			
			FACETS[i][j].sqrtA=pow(FACETS[i][j].area,0.5l);
			
			if (USE_WHICH==1)
				{
				gen_vec3_norm_to_vec2_from_vec1(XBASIS, FACETS[i][j].norm, FACETS[i][j].xbasis);
				cross(FACETS[i][j].norm, FACETS[i][j].xbasis, FACETS[i][j].ybasis);
				}
			else if (USE_WHICH == 2)
				{
				gen_vec3_norm_to_vec2_from_vec1(YBASIS, FACETS[i][j].norm, FACETS[i][j].ybasis);
				cross(FACETS[i][j].ybasis, FACETS[i][j].norm, FACETS[i][j].xbasis);
				}
			else
				{
				// creates projections of x-y vectors onto the plane defined by .norm
				for (k=0; k<3; k++)
					{
					FACETS[i][j].xbasis[k]=XBASIS[k]-FACETS[i][j].norm[k]*FACETS[i][j].norm[0];
					FACETS[i][j].ybasis[k]=YBASIS[k]-FACETS[i][j].norm[k]*FACETS[i][j].norm[1];
					}
				normalise(FACETS[i][j].xbasis);
				normalise(FACETS[i][j].ybasis);
				dxy=dot(FACETS[i][j].xbasis,FACETS[i][j].ybasis);
				
				if (dxy == 0.)
					{
					coef = 0.;
					}
				else
					{
					coef = (-1.+pow(1.-dxy,0.5l))/dxy; // all div by 2
					}
				for (k=0; k<3; k++)
					{
					temp = FACETS[i][j].xbasis[k]+coef*FACETS[i][j].ybasis[k];
					FACETS[i][j].ybasis[k] = FACETS[i][j].ybasis[k]+coef*FACETS[i][j].xbasis[k];
					FACETS[i][j].xbasis[k] = temp;
					if (isnan(FACETS[i][j].xbasis[k]))
						{
						cout<<"fuck it!";
						}
					}
				normalise(FACETS[i][j].xbasis);
				normalise(FACETS[i][j].ybasis);
				}
			}
		}
	return OK;
	}

void get_facet_normal(facet &afacet)
	{
	int i;
	for (i=0; i<3; i++)
		{
		afacet.xvec[i] = (afacet.p[1][i]+afacet.p[3][i] - afacet.p[2][i]-afacet.p[0][i])/2.;
		afacet.yvec[i] = (afacet.p[2][i]+afacet.p[3][i] - afacet.p[0][i]-afacet.p[1][i])/2.;
		}
	
	cross(afacet.xvec, afacet.yvec, afacet.norm);
	normalise(afacet.norm);
	}

