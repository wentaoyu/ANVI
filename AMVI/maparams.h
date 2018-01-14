#pragma once
#define ANGLE_FILLVALUE -999
#define G1  2.5
#define L1  1
#define C1 6.0
#define C2 7.5
#define EP 0.00001
#define PI 3.1415926
class maparams
{
public:
    //��ͼ��Ĳ�����ϢRed/NIR/Blue
    float *Red;
    float *NIR;
    float *Blue;
    float *NDVI;
    float *EVI;
    //��ͼ��ĽǶ���Ϣphi_view/phi_zenith/phi_sun
    float *viewZenith;
    float *solarZenith;
    float *RaAngle;
    unsigned char* cloudMask;
};
