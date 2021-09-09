// some basic explicit fit models...

// So-far, there are three fundamental models implemented which show up in statistics regulary:
// - a sum of exponentials
// - a gaussian distribution
// - lorentz distribution
//
// The sum of exponentials and lorentz distribution allow superpositions, and the model for lorentz distributions is the most flexible implementation. It may serve as a protoype for simimalar model

// To costumize the class
//   public template<class Real=double,class Count=unsigned int>class StandardFit
// from standardFit.h. These basic models are directly included into the class definition.
//
// Each model must be derived from the abstract interface class Model defined in standardFit.h:
//
//   class Model //pure abstract class serving as interface specification
//   {
//   public:
//     Count n; //number of parameters
//     Data& data;
//     double chiSqPerDOF;
//
//     Model(Data &dataD, Count nD)
//      : data(dataD), n(nD)
//     {  }
//
//     virtual void        gsl(const gsl_vector* gsl_vector) = 0; // interface routines to exchange parameters via gsl_vector
//     virtual gsl_vector* gsl()                             = 0;
//     virtual Real        y(Real x)                         = 0; // the fit function
//     virtual Real        dyDivDj(Real x,int j)             = 0; // the derivative of the fit function ( to calc the Jacobian)
//
//     void checkImplementationI(Count i,Array<Real> x, Real di=1) // Provided method: checks for the i-th parameter whether the function's differential quotient converges against the implemented derivative for a series of decreasing argument spacing.
//   };
//



// 1fst Model: sum of expoentials.

// specifies the fit function of the following form:
// sum [ A_i*( exp(b_i*(x+c_1_i)) + c_2_i*exp(-b_i*(x)) ) ] via A_i,b_i with i from 0 to N
// TODO: c_{1,2}_i == {+,-}inf |=> drop corresponding term

// with:

// N......nAddends
// A_i....amplitude[i]
// b_i....exponents[i]
// c_1_i..c1s[i]
// c_2_i..c2s[i]

class ModelExponentials:public Model
{
protected:

  void resize(Count nAddends)
  {
    n=nAddends*2;
    if(nAddends!=this->nAddends)
    {
      if(amplitudes!=NULL)
      {
        delete[] amplitudes;
        delete[] exponents;
        delete[] c1s; delete[] c2s;
      }
      this->nAddends=nAddends;

      if(nAddends>0)
      {
        amplitudes=new Real[nAddends]; exponents=new Real[nAddends]; c1s=new Real[nAddends]; c2s=new Real[nAddends];
      }
      else
      {
        amplitudes=NULL; exponents=NULL; c1s=NULL; c2s=NULL;
      }
    }
  }

public:
  using Model::n;// using Model::data; // Why C++, WHY?
  Count nAddends;
  Real* amplitudes;
  Real* exponents;
  Real* c1s,* c2s;
  bool coshs;


  ModelExponentials(Count nAddends,Data& dataD,bool cosh=true)
    : Model(dataD,0), nAddends(0), amplitudes(NULL), exponents(NULL), c1s(NULL), c2s(NULL), coshs(coshs)
  {
    resize(nAddends);
    for(Count i=0;i<nAddends;i++)
    {
      amplitudes[i] =0.1;   // simple initial values
      exponents[i]  =i*0.1;
      c1s[i]=0; c2s[i]=0;
    }
  }

  ~ModelExponentials()
  {  resize(0);  }


  // interface routines to exchange parameters via gsl_vector

  virtual void gsl(const gsl_vector* gsl_vector)
  {
    resize(Count(gsl_vector->size/2));
    for(Count i=0;i<nAddends;i++)
    {
	    amplitudes[i] = gsl_vector->data[2*i];
	    exponents[i]  = gsl_vector->data[2*i+1];
      //cout<<"gsl<<"<<amplitudes[i]<<" "<<exponents[i]<<endl;
    }
  }

  virtual gsl_vector* gsl()
  {
    gsl_vector* result=gsl_vector_alloc(2*nAddends);

    for(Count i=0;i<nAddends;i++)
    {
      result->data[2*i]   = amplitudes[i];
      result->data[2*i+1] = exponents[i];
    }

    return result;
  }


  // the fit function
  virtual inline Real y(Real x)
  {
    Real y=0;
    for(Count i=0;i<nAddends;i++)
      y+=amplitudes[i]*( exp(exponents[i]*(x+c1s[i])) + ( coshs ? exp(-exponents[i]*(x+c2s[i])) : 0 ) );
    //cout<<amplitudes[0]<<" "<<exponents[0]<<" | "<<c1s[0]<<" "<<c2s[0]<<"  "<<y<<" "<<this<<endl;
    return y;
  }

  // the derivative of the fit function
  // to calc the Jacobian
  virtual inline Real dyDivDj(Real x,int j)
  {
    int i=j/2;

    if(j%2==0)
      return exp(exponents[i]*(x+c1s[i])) + ( coshs ? exp(-exponents[i]*(x+c2s[i])) : 0 );
    return amplitudes[i]*( (x+c1s[i])*exp(exponents[i]*(x+c1s[i])) - ( coshs ? (x+c2s[i])*exp(-exponents[i]*(x+c2s[i])) : 0 ) );
  }
};


class ModelSingleGaussian:public Model
{
public:
  using Model::n;// using Model::data; // Why C++, WHY?
  Real amplitude, exponent, x0;//x0 is the position of the center
  bool x0Fixed;


  ModelSingleGaussian(Data &data,Real amplitude=1,Real exponent=-1, Real x0=0, bool x0Fixed=true) //x0 is the position of the center
    : Model(data,2+Count(!x0Fixed)), amplitude(amplitude), exponent(exponent), x0(x0), x0Fixed(x0Fixed)
  {  }


  // interface routines to exchange parameters via gsl_vector

  virtual void gsl(const gsl_vector* gsl_vector)
  {
                 amplitude = gsl_vector->data[0];
                 exponent  = gsl_vector->data[1];
    if(!x0Fixed) x0        = gsl_vector->data[2];
  }

  virtual gsl_vector* gsl()
  {
    gsl_vector* result=gsl_vector_alloc(n);
                   result->data[0] = amplitude;
                   result->data[1] = exponent;
      if(!x0Fixed) result->data[2] = x0;
    return result;
  }


  // the fit function
  virtual inline Real y(Real x)
  {
    return amplitude*exp(exponent*pow(x-x0,2.));
  }

  // the derivative of the fit function
  // to calc the Jacobian
  virtual inline Real dyDivDj(Real x,int j)
  {
    if(j==0)
      return exp(exponent*pow(x-x0,2.));
    if(j==1)
      return amplitude*pow(x-x0,2.)*exp(exponent*pow(x-x0,2.));
    return amplitude*exponent*(-2.)*(x-x0)*exp(exponent*pow(x-x0,2.));
  }
};



template<Count nParameter> class AdvancedModel: public Model // an generic tool for advanced models: parameters are organised into a 2dim. Array with cyclic access: [i%size] "parametersPeaks"
                                                             //                                      rows represent different categories of fit parameters
                                                             //                                      columns cycle through all representants of the corresponding category. We assume each column to belong to an individual peak in the model and,
                                                             //                                         therefore, refer to representants as peaks. Derived models may assign other meanings.
{
public:
  using Model::n;
  enum ParameterId{aId,scaleId,x0Id};//must be initialised before the other members

  template<class T>class ModuloArray:public Array<T> //index parameter is the integer modulo group {unsigned int}/Array<.>::n()
  {
  public:
    using Array<T>::Array; using Array<T>::operator=;

    T& operator[](const Count i) noexcept  {
      return Array<T>::operator[](i%Array<T>::n());  } //static_cast<Array<T>&>(*this)[i%Array<T>::n()];  } // ATTENTION: casting to non-& type from *this results in a copy for the this parameter of Array<T>::[]

    const T& operator[](const Count i) const noexcept  {
      return Array<T>::operator[](i%Array<T>::n());  }
  };

protected:
  ModuloArray<ModuloArray<Real>>    parametersPeaks; // see class description.
  const std::array<bool,nParameter> fixedParameters; // simple const carray may help for optimisations
  const Count nPeaks;

  inline Count getParameterIdAndParameterSPeakId(Count &j) // j is in and output: input : id for internal parameter handling for the  gsl interface and y(x,j) and dyDivDj(x,j) methods
                                                           //                     output: ParameterSPeakId refering to the j-th peak
  {
    Count parameterId;
      for(parameterId=0; fixedParameters[parameterId]||j>=parametersPeaks[parameterId].n(); ++parameterId)
        if(!fixedParameters[parameterId])
          j-=parametersPeaks[parameterId].n();
    return parameterId;
  }

public:

  AdvancedModel (Data &data,ModuloArray<ModuloArray<Real>> parametersPeaks, Count nPeaks,std::array<bool,nParameter> fixedParameters)
    : Model(data,[&parametersPeaks,&fixedParameters](){ Count n=0; for(Count i=0;i<parametersPeaks.n();++i) n+=Count(!fixedParameters[i])*parametersPeaks[i].n();   return n; }()),
      parametersPeaks(parametersPeaks),
      fixedParameters(fixedParameters),
      nPeaks(nPeaks)
  {
    //std::cout<<"These are the model's parameter: "; parametersPeaks.print(); cout<<endl;
  }


  // interface routines to exchange parameters via gsl_vector

  virtual void gsl(const gsl_vector *gsl_vector)
  {
    Count parameterSI0=0, parameterId=0;
      while(fixedParameters[parameterId])
        ++parameterId;

    for(auto i=0;i<n;++i)
    {
      //cout<<"(i,pid,i-i0) "<<i<<", "<<parameterId<<", "<<i-parameterSI0<<endl;
      parametersPeaks[parameterId][i-parameterSI0] = gsl_vector->data[i];

      if(i-parameterSI0==parametersPeaks[parameterId].n()-1)
      {
        do  ++parameterId;
          while(parameterId<parametersPeaks.n()&&fixedParameters[parameterId]);
        parameterSI0=i+1;
      }
    }
  }

  virtual gsl_vector* gsl()
  {
    gsl_vector *result=gsl_vector_alloc(n);
    Count       parameterSI0=0, parameterId=0;
      while(fixedParameters[parameterId])
        ++parameterId;

    for(auto i=0;i<n;++i)
    {
      result->data[i] = parametersPeaks[parameterId][i-parameterSI0];

      if(i-parameterSI0==parametersPeaks[parameterId].n()-1)
      {
        do ++parameterId;
          while(parameterId<parametersPeaks.n()&&fixedParameters[parameterId]);
        parameterSI0=i+1;
      }
    }

    return result;
  }
};






class ModelGaussian:public AdvancedModel<3>
{
protected:
  template<class T> using ModuloArray = typename AdvancedModel<3>::template ModuloArray<T>;

  Count getNPeaks(Count aSN, Count scaleSN, Count x0SN) // the func is called before the member references as &a are initialised => need of the functions parameters
  {
    Count n=( aSN > scaleSN ? aSN : scaleSN );
    return ( n > x0SN ? n : x0SN );
  }

public:
  enum ParameterId{aId,scaleId,x0Id};//must be initialised before the other members
  ModuloArray<Real> &a, &scale, &x0, // a is the amplitude, scale is sigma or std deviation, x0 is the position of the center
                    &amplitude;      // reference to a
  bool symmetrize;

  ModelGaussian(Data &data, Array<Real> a=1, Array<Real> scale=1, Array<Real> x0=0, bool aFixed=false,bool scaleFixed=false,bool x0Fixed=true,bool symmetrize=false) //x0 is the position of the center
    : AdvancedModel<3>( data,
                        {a,scale,x0}, getNPeaks(a.n(), scale.n(), x0.n()), {aFixed,scaleFixed,x0Fixed} ),
       a(this->parametersPeaks[aId]), scale(this->parametersPeaks[scaleId]), x0(this->parametersPeaks[x0Id]),
          amplitude(this->a),
      symmetrize(symmetrize)
 {  }


  // the fit function
  inline Real gaussian(Real x,Count i)
  {
    return a[i]*exp(-0.5*pow((x-x0[i])/scale[i],2));
  }
  virtual inline Real y(Real x)
  {
    Real result=0;
      for(auto i=0;i<this->nPeaks;++i)
        result += gaussian(x,i) + ( symmetrize ? gaussian(-x,i) : 0 );
    return result;
  }

  // the derivative of the fit function
  // to calc the Jacobian
  virtual inline Real dyDivDj(Real x,int j)
  {
    Count parameterId=this->getParameterIdAndParameterSPeakId(j);

    //cout<<"PID "<<parameterId<<" / j: "<<j<<endl;
    Real res=0;
      for(auto s=0; s<=symmetrize; ++s,x*=-1.)
        for(Count i=j;i<this->nPeaks;i+=this->parametersPeaks[parameterId].n())
          res +=
          ( parameterId==aId ?
              gaussian(x,i)/a[i]
            :
              ( parameterId==scaleId ?
                  gaussian(x,i)*pow(x-x0[i],2)/pow(scale[i],3)
                :
                  gaussian(x,i)*(x-x0[i])/pow(scale[i],2)
              )
          );
    return res;
  }
};




// In the case of the lorentzian model, multiple peaks are supported: pass an Array<Real> for the DOF which should be fitted individually for each peak and a plain Real for a shared DOF
class ModelLorentzian:public AdvancedModel<3>
{
protected:
  template<class T> using ModuloArray = typename AdvancedModel<3>::template ModuloArray<T>;

  Count getNPeaks(Count amplitudeSN, Count scaleSN, Count x0SN) // the func is called before the member references as &amplitude are initialised => need of the functions parameters
  {
    Count n=( amplitudeSN > scaleSN ? amplitudeSN : scaleSN );
    return ( n > x0SN ? n : x0SN );
  }

public:
  enum ParameterId{amplitudeId,scaleId,x0Id};//must be initialised before the other members
  ModuloArray<Real> &amplitude, &scale, &x0; //x0 is the position of the center of the peak
  bool symmetrize;

  ModelLorentzian(Data &data, Array<Real> amplitude=1, Array<Real> scale=1, Array<Real> x0=0, bool amplitudeFixed=false,bool scaleFixed=false,bool x0Fixed=true, bool symmetrize=false) //x0 is the position of the center
    : AdvancedModel<3>(  data,
                         {amplitude,scale,x0}, getNPeaks(amplitude.n(), scale.n(), x0.n()), {amplitudeFixed,scaleFixed,x0Fixed}  ),
      amplitude(this->parametersPeaks[amplitudeId]), scale(this->parametersPeaks[scaleId]), x0(this->parametersPeaks[x0Id]),
      symmetrize(symmetrize)
  {  }


  // the fit function
  inline Real lorentzian(Real x,Count i)
  {
    return abs(amplitude[i])/( M_PI*abs(scale[i])*(1.+pow((x-x0[i])/scale[i],2.)) );
  }
  virtual inline Real y(Real x)
  {
    Real result=0;
      for(auto i=0;i<this->nPeaks;++i)
        result += lorentzian(x,i) + ( symmetrize ? lorentzian(-x,i) : 0 );
    return result;
  }

  // the derivative of the fit function
  // to calc the Jacobian
  virtual inline Real dyDivDj(Real x,int j)
  {
    Count parameterId=this->getParameterIdAndParameterSPeakId(j);

    //cout<<"PID "<<parameterId<<" / j: "<<j<<endl;
    Real res=0;
      for(auto s=0; s<=symmetrize; ++s,x*=-1.)
        for(Count i=j;i<this->nPeaks;i+=this->parametersPeaks[parameterId].n())
          res +=
          ( parameterId==amplitudeId ?
              ( amplitude[i]<0 ? -1. : 1. )/( M_PI*abs(scale[i])*(1.+pow((x-x0[i])/scale[i],2.)) )
            :
              ( parameterId==scaleId ?
                  lorentzian(x,i)*(-1./scale[i]+lorentzian(x,i)*2.*M_PI/abs(amplitude[i])*pow((x-x0[i])/scale[i],2.)*( scale[i]<0 ? -1. : 1. ))
                :
                  pow(lorentzian(x,i),2.)*2.*M_PI/(abs(amplitude[i])*abs(scale[i]))*(x-x0[i])
              )
          );
    return res;
  }
};
