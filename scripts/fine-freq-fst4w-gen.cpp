// fine-freq-fst4w-gen: synthetic FST4W WAV at an exact frequency, independent double-precision GFSK.
// usage: fine-freq-fst4w-gen PERIOD F_LOWEST_TONE_HZ SNR_DB SEED OUT.wav [fspread_hz]
// Prints the message used to stdout.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#include "fst4w_ldpc.h"
#include "fst4w_unpack.h"
using namespace fst4w;
static int nsps_of(int p){switch(p){case 120:return 8200;case 300:return 21504;case 900:return 66560;case 1800:return 134400;}return 0;}
int main(int argc,char**argv){
  if(argc<6){fprintf(stderr,"usage\n");return 2;}
  int period=atoi(argv[1]); double f=atof(argv[2]); double snr=atof(argv[3]); unsigned seed=strtoul(argv[4],0,10); const char*out=argv[5];
  int nsps=nsps_of(period); if(!nsps){fprintf(stderr,"bad period\n");return 2;}
  std::mt19937_64 rng(seed);
  // message: random 50-bit payload that unpacks as a plain type-1 message
  int8_t m74[74]; std::string msg;
  for(int tries=0;;++tries){
    uint8_t pl[50]; for(int i=0;i<50;i++) pl[i]=rng()&1;
    Unpacker u; auto r=u.unpack(pl);
    if(!r.ok) continue; std::string s=r.msg; while(!s.empty()&&s.back()==' ')s.pop_back();
    if(s.find('<')!=std::string::npos||s.find('/')!=std::string::npos) continue;
    int sp=0; for(char c:s) sp+=c==' '; if(sp!=2) continue;
    for(int i=0;i<50;i++) m74[i]=pl[i]; for(int i=50;i<74;i++) m74[i]=0;
    int crc=get_crc24(m74,74); for(int i=0;i<24;i++) m74[50+i]=(crc>>(23-i))&1;
    if(get_crc24(m74,74)!=0){fprintf(stderr,"crc fail\n");return 3;}
    msg=s; break;
  }
  int8_t cw[240]; encode240_74(m74,cw);
  const int NN=160, ND=120; int itone[NN]; int itmp[ND];
  for(int i=1;i<=ND;i++){int is=cw[2*i-1]+2*cw[2*i-2]; itmp[i-1]= is<=1?is:(is==2?3:2);}
  static const int s1[8]={0,1,3,2,1,0,2,3}, s2[8]={2,3,1,0,3,2,0,1};
  for(int i=0;i<8;i++){itone[i]=s1[i];itone[38+i]=s2[i];itone[76+i]=s1[i];itone[114+i]=s2[i];itone[152+i]=s1[i];}
  for(int i=0;i<30;i++){itone[8+i]=itmp[i];itone[46+i]=itmp[30+i];itone[84+i]=itmp[60+i];itone[122+i]=itmp[90+i];}
  const double fs=12000, baud=fs/nsps, twopi=2*M_PI;
  const double c=M_PI*std::sqrt(2.0/std::log(2.0)), bt=2.0;
  std::vector<double> pulse(3*nsps);
  for(int i=1;i<=3*nsps;i++){double t=(i-1.5*nsps)/nsps; pulse[i-1]=0.5*(std::erf(c*bt*(t+0.5))-std::erf(c*bt*(t-0.5)));}
  const long nmax=(long)period*12000; const long i0=(long)(1.0*fs); // nominal start 1 s
  std::vector<double> x(nmax,0.0);
  const double A=std::sqrt(2.0*2500.0/6000.0*std::pow(10.0,snr/10.0));
  double phi=std::uniform_real_distribution<double>(0,twopi)(rng);
  const long nsig=(long)NN*nsps;
  for(int ib=1;ib<=NN;ib++){
    for(int m=0;m<nsps;m++){
      double d=0; for(int j=std::max(1,ib-1);j<=std::min(NN,ib+1);j++) d+=pulse[(ib-j+1)*nsps+m]*itone[j-1];
      double inst=f + d*baud; // lowest tone f, tone k at f+k*baud
      long k=(long)(ib-1)*nsps+m; double env=1.0;
      if(k<nsps/4) env=(1-std::cos(twopi*k/(nsps/2.0)))/2;
      long kd=nsig-1-k; if(kd<nsps/4) env=(1-std::cos(twopi*kd/(nsps/2.0)))/2;
      if(i0+k<nmax) x[i0+k]=A*env*std::cos(phi);
      phi+=twopi*inst/fs; if(phi>twopi) phi-=twopi;
    }
  }
  std::normal_distribution<double> g(0.0,1.0);
  std::vector<int16_t> s(nmax);
  for(long i=0;i<nmax;i++){double v=100.0*(x[i]+g(rng)); if(v>32767)v=32767; if(v<-32768)v=-32768; s[i]=(int16_t)std::lrint(v);}
  FILE*o=fopen(out,"wb"); uint32_t db=nmax*2, rs=36+db, sr=12000, br=24000; uint16_t fmt=1,ch=1,ba=2,bps=16; uint32_t fl=16;
  fwrite("RIFF",1,4,o);fwrite(&rs,4,1,o);fwrite("WAVEfmt ",1,8,o);fwrite(&fl,4,1,o);fwrite(&fmt,2,1,o);fwrite(&ch,2,1,o);
  fwrite(&sr,4,1,o);fwrite(&br,4,1,o);fwrite(&ba,2,1,o);fwrite(&bps,2,1,o);fwrite("data",1,4,o);fwrite(&db,4,1,o);fwrite(s.data(),2,nmax,o);fclose(o);
  printf("%s\n",msg.c_str());
  return 0;
}
