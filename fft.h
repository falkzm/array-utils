//semi-efficient FFT implementation for use with Array<> system
//memory usage ~2*Memory(input data) assuming input is of same type as output
//Input Type T -> input Array is of Type Array T; output Type Tf => output Array is of Type Array< complex<Tf> >

//use FFTW or other library std::Types

//pScale defines the momentum scale: in std DFT the momentum resolution does not scale with the resolution of the input data but the maximal momentum increases proportionally
//  this behaviour invoked via default parameter pScale=1
//  in order to rescale momentum and input data proportionally, set pScale to the inverse scaling of the input data
//  the factor is absorbed into the exponentials of the "twiddle factors"/wave function space

#include "linearAlgebra.h"
#include <complex>
#include <cmath>


template<class T=double,class Tf=T >class Dft
{
public:

  class Fft
  {
  protected:
    Array<T> a;//access this input data field only as const object in order to prevent data copying due to multiple, manipulative usage.
    double pScale;
    double pScaleInverse;

    Array< complex<Tf> > fft(Count s, Count o)
      //s and o are the spacing and offset for data access of the input field a
    {
      if(s==a.n())
        return Array< complex<Tf> >(ceil(pScaleInverse),static_cast< complex<Tf> >(static_cast< const Array<T> >(a)[o]),NULL);
        //return static_cast< complex<Tf> >(static_cast< const Array<T> >(a)[o%a.n()]);
      else
      {
         Count           n =ceil(double(a.n())*pScaleInverse/double(2*s)); // n/2
         complex<double> ww=complex<double>(0,-M_PI/n);                    //"twiddle" factor exponential.
         Array< complex<Tf> > even=fft(s*2,o),
                              odd =fft(s*2,o+s),
                              b(2*n,NULL);
         //if(o==0)cout<<"GG "<<n<<endl;

         for(Count p=0;p<n;++p)
         {
           complex<double> w=exp(ww*double(p));
           //b[  p]=even[p]+odd[p]*w;
             b[  p]=even[p]+complex<Tf>(odd[p].real()*w.real()-odd[p].imag()*w.imag(),odd[p].imag()*w.real()+odd[p].real()*w.imag());
           //b[n+p]=even[p]-odd[p]*w;
             b[n+p]=even[p]-complex<Tf>(odd[p].real()*w.real()-odd[p].imag()*w.imag(),odd[p].imag()*w.real()+odd[p].real()*w.imag());
         }

         return b;
      }
    }

  public:

    Array< complex<Tf> > calc_2PowN(const Array<T> a_,Tf pScale_=1)
    {
      if(!is2PowN(a_.n()))
      {
        std::cerr<<"ERROR: Fft::calc_2PowN called for field with length != 2^n."<<endl;
        return Array< complex<Tf> >();// a blank field
      }
      a=a_; pScale=pScale_; pScaleInverse=1./double(pScale);

      return fft(1,0);
    }


    static bool is2PowN(Count x)
    {
      Count y=1;
        for(;y<x;y*=2);
      return y==x;
    }


    //there is no naive stretching or truncating mechanism to generalise the FFT to abitrary lengths n - doing so would add artificial frequencies to generate the artifical constant intervals - these will occur in the whole spectrum.

    static complex<Tf> interpolateToContineous(const Array< complex<Tf> > b, double x,double pScale=1)
      //interpolate between discrete data points at x using the fourier coefficients b (as provided from calc*)
      //this is the plain inverse DFT.
      //there should be a O(n log n) approach to double the resolution via an inverse fft...
    {
      complex<Tf> a =0,
                  ww=complex<Tf>(0,M_PI*Tf(2)/Tf(b.n())*pScale);            //"twiddle" factor exponential.
      for(Count p=0;p<b.n();++p)
        a+=exp(ww*Tf(p))*b[p];
      return a/Tf(b.n());
    }
  };


  static complex<Tf> calc(const Array<T> a, double p, double pScale)
  {
    complex<Tf>     b(a[0]*0);
    complex<double> ww(0,-M_PI*2./a.n()*pScale);
    for(Count x=0;x<a.n();++x)
      b+=complex<Tf>(a[x]*exp(ww*double(x)*p).real(),a[x]*exp(ww*double(x)*p).imag());
    return b;
  }

  static Array< complex<Tf> > calc(const Array<T> a,double pScale=1)
  {
    if(Fft::is2PowN(a.n()))
    {  Fft fft; return fft.calc_2PowN(a,pScale);  }
    else
    {
      Array< complex<Tf> > b(a.n(),NULL);
      for(Count p=0;p<b.n();++p)
        b[p]=calc(a,p,pScale);
      return b;
    }
  }


  static complex<Tf> interpolateToContineous(const Array< complex<Tf> > b, double x,Tf pScale=1)
  {  return Fft::interpolateToContineous(b,x,pScale);  }


  static void test()
  {
    std::cout<<"FFT-Test: a(length=1024); a[i]=sin(2Pi i/10) -> T=10,f=1/10; arg_max_j( calc(a)[j] ) = "<<std::flush;
    Array<T> a(1024,NULL); for_vec(a,i) a[i]=std::sin(M_PI*2./10.*i);
    auto b=calc(a);
    Count j=0; for_vec(b,i) if(i<b.n()/2 && abs(b[i])>abs(b[j])) j=i;
    std::cout<<j<<"\n"<<std::flush; // 102 // => =L/T
  }

};
