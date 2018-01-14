#include "QualityCal.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <algorithm>


QualityCal::QualityCal()
{
}

QualityCal::QualityCal(float * RedData, float * NirData, float * SolarZenith, float * ViewZenith, float * RaAngle,
	unsigned char landcover, int sizeMODIS, int size, unsigned char * cloudMASK)
{
	
	n1 = n2 = n3 = n4 = 0;
	
	for (int i = 0; i < size; i++)
	{
		if (RedData[i] <= 0.001 || RedData[i] >= 1 || NirData[i] <= 0.001 || NirData[i] >= 1 || SolarZenith[i] < 0 ||
			ViewZenith[i] < 0 || ViewZenith[i] > 65 || (RaAngle[i] + 999) < EP || landcover == 0 || cloudMASK[i] == 1)
		{
		//	qcref_new[i] = 5;
		}
		else
		{
			refb1_c.push_back(RedData[i]);
			refb2_c.push_back(NirData[i]);
			ndvi_c.push_back((NirData[i] - RedData[i]) / (NirData[i] + RedData[i]));
			sangle_c.push_back(SolarZenith[i]);
			vangle_c.push_back(ViewZenith[i]);
			rangle_c.push_back(RaAngle[i]);
			tagT1.push_back(i);
			qcref_new.push_back(5);
			
			n1++;
			n2++;
			if (i < sizeMODIS)
			{
				countMODIS++;
			}
		}

		if (!(RedData[i] <= 0.001 || RedData[i] >= 1 || NirData[i] <= 0.001 || NirData[i] >= 1 || SolarZenith[i] < 0 ||
			ViewZenith[i] < 0 || ViewZenith[i] > 45 || (RaAngle[i] + 999) < EP || landcover == 0 || cloudMASK[i] == 1))
		{
			refb1_a1.push_back(RedData[i]);
			refb2_a1.push_back(NirData[i]);
			ndvi_a1.push_back((NirData[i] - RedData[i]) / (NirData[i] + RedData[i]));
			n3++;
		}
		if (!(RedData[i] <= 0.001 || RedData[i] >= 1 || NirData[i] <= 0.001 || NirData[i] >= 1 || SolarZenith[i] < 0 ||
			ViewZenith[i] < 0 || ViewZenith[i] > 45 || (RaAngle[i] + 999) < EP || landcover == 0 ))
		{
			refb1_a2.push_back(RedData[i]);
			refb2_a2.push_back(NirData[i]);
			ndvi_a2.push_back((NirData[i] - RedData[i]) / (NirData[i] + RedData[i]));
			n4++;
		}
	}//end for

	if (n3 >= 1)
	{
		ndvi_mvc_mask = max_k_select(ndvi_a1, 1);
	}
	if (n3 >= 2)
	{
		//cv_mvc??
	}
	else if (n3 == 1)
	{
		ndvi_cv = ndvi_a1[0];
	}
	if (n4 >= 1)
	{
		ndvi_mvc_all = max_k_select(ndvi_a2, 1);
	}

}


QualityCal::~QualityCal()
{
}

float QualityCal::mean(vector<float> data)
{
	float sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += data[i];
	}
	return sum / data.size();
}

void QualityCal::Process(float& newvi,float& newvi_qc)
{
	int countMODIST2 = countMODIS;
	if (n1 >= 5)
	{
		vector<float> RedDataT2, NirDataT2, SZ_T2, VZ_T2, RA_T2, NDVI_T2;//第一次剔除无效值
		vector<int> tagT2;
		float ndvi_m = max_k_select(ndvi_c, 2);   //2nd largest
		for (int i = 0; i < n1; i++)
		{
			if (ndvi_c[i] < (ndvi_m - 0.3) || ndvi_c[i] < 0.01)
			{
				qcref_new[i] = 3;
				n2--;
			}
			else
			{
				RedDataT2.push_back(refb1_c[i]);
				NirDataT2.push_back(refb2_c[i]);
				SZ_T2.push_back(sangle_c[i]);
				VZ_T2.push_back(vangle_c[i]);
				RA_T2.push_back(rangle_c[i]);
				NDVI_T2.push_back(ndvi_c[i]);
				tagT2.push_back(tagT1[i]);
			}
		}
		float ndvi_m_c = mean(NDVI_T2);  //mean ndvi

		if (n2 >= 5)
		{
			BRDF_Total(RedDataT2, NirDataT2, SZ_T2, VZ_T2, RA_T2, NDVI_T2);
			if (refnew_b1 > 0 && refnew_b1 < 1 && refnew_b2>0 && refnew_b2 < 1 && newndvi > 0 && newndvi < 1)
			{
				newvi = newndvi;
				newvi_qc = 0;
				
			}
			else
			{

			}
		}
		else if (n2 > 4)
		{

		}

	}
}


float QualityCal::max_k_select(vector<float> a, int k)  //select the kth max number
{
	sort(a.begin(), a.end()); //asceding sequence;
	return (a[a.size()-k]);

}

void QualityCal::BRDF_Total(fvec RedDataT2, fvec NirDataT2, fvec SZ_T2, fvec VZ_T2, fvec RA_T2, fvec NDVI_T2)
{
	
	int m = SZ_T2.size();   //
	int n = 3;
	float ndvi_m_c = sum(NDVI_T2) / m;
	fvec kvol(m);
	fvec kgeo(m);
	fvec kiso(m, fill::ones);
	brdf_forward_RLM(SZ_T2, VZ_T2, RA_T2, kvol, kgeo);
	fmat B_mds = {join_rows(join_rows(kiso,kvol),kgeo)};   //refer to A m*3 
	fmat lb1, lb2;
	lb1 = RedDataT2;
	lb2 = NirDataT2;

	fmat Pllb1(m, m, fill::zeros), Pllb2(m, m, fill::zeros);
	for (int m1 = 0; m1 < m; m1++)
	{
		Pllb1(m1, m1) = NDVI_T2[m1] / ndvi_m_c;
		Pllb2(m1, m1) = NDVI_T2[m1] / ndvi_m_c;
	}
	fvec xils_b1(n),xils_b2(n);
	frowvec Pfnew_b1(m), Pfnew_b2(m);
	fmat refnewb1, refnewb2;
	fvec  ndvi_mds_cal;
	fmat Pllb1_tmp, Pllb2_tmp;
	for (int iter = 0; iter < 10; iter++)
	{
		iterls(B_mds, Pllb1, lb1, xils_b1, Pfnew_b1);
		iterls(B_mds, Pllb2, lb2, xils_b2, Pfnew_b2);
		/*cout << xils_b1;
		cout << xils_b2;*/
		refnewb1 = B_mds*xils_b1;
		refnewb2 = B_mds*xils_b2;
	
		ndvi_mds_cal = (refnewb2 - refnewb1)/(refnewb2 + refnewb1);
		cout << ndvi_mds_cal;
		fvec Pn = NDVI_T2/ndvi_mds_cal;
		cout << NDVI_T2;
		for (int mm = 0; mm < m; mm++)
		{
			if (Pn(mm) <= 0) 
			{
				Pn(mm) = 0.1;
			}
		}
		cout << Pn;
		Pllb1_tmp = Pllb1;
		Pllb2_tmp = Pllb2;
		for (int m1 = 0; m1 < m; m1++)
		{
			Pllb1(m1, m1) = Pfnew_b1(m1)*Pn(m1);
			Pllb2(m1, m1) = Pfnew_b2(m1)*Pn(m1);
		}
		if (iter >= 1 && (abs((Pllb1 - Pllb1_tmp).max()) < 0.001)&&(abs((Pllb2 - Pllb2_tmp).max())<0.001))
		{
			break;
		}
	}//end iter

	
	fvec newkvol(1), newkgeo(1);
	//angle normalization
	fvec newsz(1),newvz(1),newra(1);
	newsz(0) = 45.;
	newvz(0) = 0.;
	newra(0) = 0.;
	brdf_forward_RLM(newsz,newvz,newra, newkvol, newkgeo);

	refnewb1 = 1 * (xils_b1[0]) + newkvol * (xils_b1[1]) + newkgeo * (xils_b1[2]);
	refnewb2 = 1 * (xils_b2[0]) + newkvol * (xils_b2[1]) + newkgeo * (xils_b2[2]);
	refnew_b1 = refnewb1(0, 0);
	refnew_b2 = refnewb2(0, 0);
//	cout << refnewb1 << "****" << refnewb2 << endl;
	newndvi = (refnew_b2 - refnew_b1) / (refnew_b1 + refnew_b2);
}

void QualityCal::brdf_forward_RLM(fvec SZ, fvec VZ, fvec RA, fvec& kvol, fvec& kgeo)
{
	int n = SZ.size();
	for (int i = 0; i < n; i++)
	{
		float sz, vz, ra;
		sz = SZ[i] * PI / 180.0;
		vz = VZ[i] * PI / 180.0;
		ra = RA[i] * PI / 180.0;
		kvol[i] = Ross_thick(sz, vz, ra);
		kgeo[i] = Li_Transit(sz, vz, ra);
	}
	
}


//小窦
float QualityCal::Ross_thick(float sz, float vz, float relaz)
{
	float K_vol = 1, cosxi, xi, cosvz, cossz;
	cosvz = cos(vz);
	cossz = cos(sz);
	cosxi = cossz*cosvz + sin(sz)*sin(vz)*cos(relaz);
	xi = acos(cosxi);
	K_vol = ((PI / 2 - xi)*cosxi + sin(xi))*(1 + 1 / (1 + fabs(xi) / (PI*1.5 / 180.0))) / (cossz + cosvz) - PI / 4;
	return K_vol;
}
//小窦
float QualityCal::Li_Transit(float sZenith, float vZenith, float rAzimuth)
{
	float brratio = 1;  //	b/r=1:
	float hbratio = 2;  //	h/r=2:
	float t1, theta_ip, t2, theta_vp, temp1, temp2, cosxip, D1, D, cost1, cost2, temp3, cost, t, O, B, k;
	t1 = brratio*tan(sZenith);
	theta_ip = atan(t1);
	t2 = brratio*tan(vZenith);
	theta_vp = atan(t2);
	temp1 = cos(theta_ip);
	temp2 = cos(theta_vp);
	cosxip = temp1*temp2 + sin(theta_ip)*sin(theta_vp)*cos(rAzimuth);
	D1 = tan(theta_ip)*tan(theta_ip) + tan(theta_vp)*tan(theta_vp) - 2 * tan(theta_ip)*tan(theta_vp)*cos(rAzimuth);
	D = sqrt(D1);
	cost1 = tan(theta_ip)*tan(theta_vp)*sin(rAzimuth);
	cost2 = D1 + cost1*cost1;
	temp3 = 1 / temp1 + 1 / temp2;
	cost = hbratio*sqrt(cost2) / temp3;

	if (cost > 1)
		cost = 1;
	t = acos(cost);
	O = (t - sin(t)*cost)*temp3 / PI;
	B = temp3 - O;
	//    if (B > 2)
	//        k = (1+cosxip)/(temp2*temp1*B)-2;
	//    else
	k = -B + (1 + cosxip) / (2 * temp2*temp1);
	return k;
}

void QualityCal::iterls(fmat B, fmat Pll, fmat l, fvec& xils, frowvec& Pfnew)
{
	cout << B;
	cout << Pll;
	cout << l;

	int m = Pll.n_rows;
	int n = Pll.n_cols;
	fmat Qll(n,n,fill::zeros);
	for (int i = 0; i < n; i++)
	{
		Qll(i, i) = 1 * 1.0 / Pll(i, i);
	}
	cout << "Qll：" << Qll;
	fmat ibtpb;
	int fl = arma::inv(ibtpb, B.t()*Pll*B);
	if (fl == 0)
	{
		ibtpb.set_size(3, 3);
		ibtpb.fill(1);
		cout << ibtpb;
	}
	fmat Qvv = Qll - B*ibtpb*B.t();
	fmat QP = Qvv*Pll;
	frowvec rs(n, fill::zeros);
	for (int i = 0; i < n; i++)
	{
		rs(i) = QP(i, i);
	}
	float r = trace(QP);
	
	fmat x = ibtpb* (B.t()*Pll*l);
	xils = x;
	fmat V = B*x - l;
	
	fmat sig02 = (V.t()*Pll*V) / r;
	fmat sig2_sum = (V.t()*V) / r;
	frowvec sig2(n, fill::zeros);
	for (int i = 0; i < n; i++)
	{
		sig2(i) = pow(V(i, 0), 2) / rs(i);
	}
	Pfnew.fill(fill::zeros);
	float T;
	for (int i = 0; i < n; i++)
	{
		T = sig2(i) / sig02(0,0);
		int tagID = ((int)(r + 0.5)) - 1;
		if (tagID >= 100 || tagID <0)
		{
			Pfnew[i] = 1;
			continue;
		}
		if (T < finvData[((int)(r + 0.5)) - 1])
		{
			Pfnew[i] = 1;
		}
		else
		{
			if (V[i] > -0.00001 && V[i] < 0.00001)
			{
				Pfnew[i] = 1;
			}
			else {
				Pfnew[i] = sig02(0,0)*rs[i] / (V[i] * V[i]);
			}
		}
	}


}





