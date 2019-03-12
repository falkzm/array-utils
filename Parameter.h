#include <iostream>
using std::cout;
//using std::printf;


typedef long double Real;
  Real absReal(const Real real)
  {
      return real<0 ? -real : real;
  }
typedef int         Count;
//componentenId's & zugehoerige Schleifen: griechisch^=4d RaumZeit, lateinisch^= 3d EichfeldRaum
#define for_spaceComponentId(mu)      for(int mu=0;mu<dimensionSCount;mu++)
#define for_gaugeFieldComponentId(i)  for(int i=0;i<gaugeFieldSCount;i++)

#define for_i         for(int i=0;i<gaugeFieldSCount;i++) /*for_gaugeFieldComponentId(i)*/
#define for_j         for(int j=0;j<gaugeFieldSCount;j++) /*for_gaugeFieldComponentId(j)*/
#define for_i_j       for_i for_j
#define for_mu        for_spaceComponentId(mu) /*for_spaceTimeComponentId(mu)*/
#define for_nu        for_spaceComponentId(nu) /*for(int nu=0;nu<dimensionSCount;nu++)*/
#define for_mu_nu     for_mu for_nu
#define for_mu_i      for_mu for_i
#define for_mu_nu_i   for_mu_nu for_i
#define for_mu_nu_i_j for_mu_nu for_i_j

#define processorCount 8


class CMonteCarlo //used by dataAnalysis
{
public:
  inline static Real random(Real maximum=1)
  {return static_cast<Real>(rand())/static_cast<Real>(RAND_MAX)*maximum;}

  inline static void randomSSeed(unsigned int seed)
  {srand(seed);}
};


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
