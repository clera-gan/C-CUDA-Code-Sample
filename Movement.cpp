#include "dempacking.h"
//#include "Movement.h"

int IsTimeToTranslateHost(uint it, double real_dt, uint NumOfTranslation)
{
 if(NumOfTranslation>0)
{
 uint i;
 for(i=0;i<NumOfTranslation;i++)
 {
  if(it> (int)(m_hTranslation[i].starttime/real_dt) && it<= (int)(m_hTranslation[i].endtime/real_dt) )
  return i;
 }
 return NumOfTranslation;
 }
 else 
 return NumOfTranslation;
}

 int IsTimeToRotateHost(uint it, double real_dt, uint NumOfRotation)
{
 if(NumOfRotation>0)
{
 uint i;
 for(i=0;i<NumOfRotation;i++)
 {
  if(it> (int)(m_hRotation[i].starttime/real_dt) && it<= (int)(m_hRotation[i].endtime/real_dt) )
 return i;
 }
 return NumOfRotation;

 }
 else 
 return NumOfRotation;

}


 int IsTimeToVibrateHost(uint it, double real_dt, uint NumOfVibration)
{
 if(NumOfVibration>0)
{
 uint i;
 for(i=0;i<NumOfVibration;i++)
 {
  if(it> (int)(m_hVibration[i].starttime/real_dt) && it<= (int)(m_hVibration[i].endtime/real_dt) )
 return i;
 }
 return NumOfVibration;
 }
 else 
 return NumOfVibration;
}

void RotationPositionConvertHost(double posix,double posiy,double posiz,
                                 double nWx,double nWy,double nWz, double theta,
                                 double *newposx,double *newposy,double *newposz)
{
    double c = cos(theta);
    double s = sin(theta);
    *newposx = (nWx*nWx*(1 - c) + c) * posix + (nWx*nWy*(1 - c) - nWz*s) * posiy + (nWx*nWz*(1 - c) + nWy*s) * posiz;
    *newposy = (nWy*nWx*(1 - c) + nWz*s) * posix + (nWy*nWy*(1 - c) + c) * posiy + (nWy*nWz*(1 - c) - nWx*s) * posiz;
    *newposz = (nWx*nWz*(1 - c) - nWy*s) * posix + (nWy*nWz*(1 - c) + nWx*s) * posiy + (nWz*nWz*(1 - c) + c) * posiz;
}
