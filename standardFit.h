#pragma once

#define HAVE_INLINE
#define GSL_RANGE_CHECK_OFF

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit_nlin.h>


// performs exponential fits:
// sum [ A_i*( exp(b_i*(x+c_1_i)) + c_2_i*exp(-b_i*(x)) ) ] via A_i,b_i with i from 0 to N
// see class ModelExponentials for details
template<class Real=double,class Count=unsigned int>class StandardFit
{
public:
  #define GSL_HAPPY GSL_SUCCESS


  class Data
  {
  public:
    Count n;
    Real *xs, *ys, *sigmas;

    Data(Count n,bool useXs=false)
    {
      this->n=n;
      xs=( useXs ? new Real[n] : NULL); ys=new Real[n]; sigmas=new Real[n];

      for(Count i=0;i<n;i++)
      {
        ys[i]     =0;
        sigmas[i] =0;
      }
    }
    template<class T> Data(const Array<T> xs, const Array<T> ys, const Array<T> sigmas) //the positions x can be set to zero, or any Array shorter than ys: in this case the ys are assumed to be located equidistant at 0,1, ... , n-1
    {
      n=ys.n();
      this->xs=( xs.n()>=n ? new Real[n] : NULL ); this->ys=new Real[n]; this->sigmas=new Real[n];

      for(Count i=0;i<n;i++)
      {
        if(this->xs!=NULL)  this->xs[i]     = xs[i];
                            this->ys[i]     = ys[i];
                            this->sigmas[i] = sigmas[i];
      }
    }

    ~Data()
    {  delete[] xs; delete[] ys; delete[] sigmas;  }
  };


  class Model //pure abstract class serving as interface specification
  {
  public:
    Count n; //number of parameters 
    Data& data;
    double chiSqPerDOF;

    Model(Data &dataD, Count nD)
     : data(dataD), n(nD)
    {  }

    virtual void        gsl(const gsl_vector* gsl_vector) = 0; // interface routines to exchange parameters via gsl_vector
    virtual gsl_vector* gsl()                             = 0;
    virtual Real y(Real x)                         = 0; // the fit function
    virtual Real dyDivDj(Real x,int j)             = 0; // the derivative of the fit function ( to calc the Jacobian)

    void checkImplementation(Count i,Array<Real> x, Real di=1)  // checks for the i-th parameter whether the function's differential quotient converges against the implemented derivative for a series of decreasing argument spacing.
    {
      auto *gsl_vectorCopy=gsl(),
           *gsl_vector    =gsl();
      Real y1=0, dy=0;
          for(auto j=0;j<x.n(); ++j)
          {  y1+=y(x[j]); dy+=dyDivDj(x[j],i);  }
 
      std::cout<<"testing convergence (y(x+dx)-y(x))/(y'(x)*dx)->1 for i=="<<i<<": ";
      for(Count e=20;e>=-20;e-=2)
      {
        gsl_vector->data[i]=gsl_vectorCopy->data[i]+di*pow(10.,e);
          gsl(gsl_vector);
        Real y2=0;
          for(auto j=0;j<x.n(); ++j)
            y2+=y(x[j]);

         std::cout<<" "<<(y2-y1)/(dy*di*pow(10.,e));//<<"("<<y1<<","<<y2<<","<<dy*di*pow(10.,e)<<")";
      }
      gsl_vector->data[i]=gsl_vectorCopy->data[i];
        gsl(gsl_vector);
      std::cout<<"\n";

      delete gsl_vector, gsl_vectorCopy;
    }
  };


protected:
  // a lot of gsl stuff ...

  static int gsl_model(const gsl_vector* x, void* modelP,  gsl_vector* f)
  {
    Model& model=*static_cast<Model*>(modelP);
      model.gsl(x);
    Data& data=model.data;

    for(Count i=0;i<data.n;++i)
    {
      gsl_vector_set
      (
        f,i,
        (model.y( data.xs!=NULL ? data.xs[i] : i )-data.ys[i])/data.sigmas[i]
      );
    }

    return GSL_HAPPY;
  }


  static int gsl_modelJacobian(const gsl_vector* x, void* modelP,  gsl_matrix* J)
  {
    Model& model=*static_cast<Model*>(modelP);
      model.gsl(x);
    Data& data=model.data;

    for(Count d=0;d<data.n;++d)
      for(Count j=0;j<model.n;++j) //2 for amplitude and exponent
        gsl_matrix_set( J,d,j, model.dyDivDj(( data.xs!=NULL ? data.xs[d] : d ),j)/data.sigmas[d] );

    return GSL_HAPPY;
  }


  static int gsl_modelAndJacobian(const gsl_vector* x, void* modelP, gsl_vector* f, gsl_matrix* J)
  {
    gsl_model(x,modelP,f);
    gsl_modelJacobian(x,modelP,J);

    return GSL_HAPPY;
  }


  static void printState(int nIteration, gsl_multifit_fdfsolver* s, int status)
  {
    cout<<"iter: "<<nIteration<<" x: ";
    for(Count i=0;i<s->x->size;i++)
    {
      cout<<gsl_vector_get(s->x,i)<<"  ";
    }
    printf ("%s\n", gsl_strerror (status));
  }

public:

  //some basic explicit models...
  #include"standardFitSModels.h"


  // the lovely main routine ...

  static void calc(Model& model,int verbosity=1)
  {
    int nIterMax=10000;

    Data& data=model.data;
    if(verbosity>0)  printf("Fitting ...\n");
    if(verbosity>2)  for(auto i=0;i<data.n;++i)
                       std::cout<<( data.xs==NULL ? Real(i) : data.xs[i] )<<"   "<<data.ys[i]<<"   "<<data.sigmas[i]<<"\n";
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver            *s = gsl_multifit_fdfsolver_alloc (T, data.n, model.n);
    gsl_multifit_function_fdf f;

    gsl_vector &x = *(model.gsl());
    //double parametersInitial[model.n];
    //  gsl_vector *parametersInitial_gsl=model.gsl();
    //  for(Count p=0;p<model.n;p++)
    //    parametersInitial[p] = parametersInitial_gsl[p]; 
    //gsl_vector_view x = gsl_vector_view_array (parametersInitial, model.n);

    f.f      = &gsl_model;
    f.df     = &gsl_modelJacobian;
    f.fdf    = &gsl_modelAndJacobian;
    f.n      = data.n;
    f.p      = model.n;
    f.params = &model;
 
    gsl_multifit_fdfsolver_set (s, &f, &x);

    {
      int nIter=0,status=0;
      if(verbosity>1)  printState(nIter,s,status);
      do
      {
        nIter++;
        status = gsl_multifit_fdfsolver_iterate (s);
     
        if(verbosity>1)  printState (nIter,s,status);
  
        if (status)
          break;
     
        status = gsl_multifit_test_delta (s->dx, s->x, 1e-8, 1e-8);
      }
      while (status == GSL_CONTINUE && nIter < nIterMax);
      if(verbosity>=0&& nIter >= nIterMax ) printf("WARNING: fit did not converge within limit of iteration nr.\n");
    }

    {
      Real chi = gsl_blas_dnrm2(s->f);
      Real dof = data.n-model.n;
        model.chiSqPerDOF = chi*chi / dof;
        if(verbosity>0)  printf("chisq/dof = %g\n", model.chiSqPerDOF);
    }
    gsl_multifit_fdfsolver_free (s);
    if(verbosity>0)  printf("Fitting done.\n");
  }
};



/* An example:

void fit_removeExcitedStates(const gsl_matrix* bootstrap_corr,unsigned int t0,unsigned int t1,Real initialExponent,gsl_matrix* bootstrap_corr_relaxed)
{
  //gsl_matrix_alloc(bootstrapping, 2*T+4)
  cout<<"defining exp fit ["<<t0<<","<<t1<<"]...\n";
  if(t0>=bootstrap_corr->size2) t0=bootstrap_corr->size2-1; if(t1>=bootstrap_corr->size2) t1=bootstrap_corr->size2-1; if(t0>t1) t1=t0;
  NonLinearFit::Data data(t1-t0+1);

  printf("preparing corr data...\n");
  for (unsigned int t=0; t<data.n; t++)
  {
    cout<<t<<endl;
    data.ys[t]    =gsl_matrix_get(bootstrap_corr,0,t+t0);
    data.sigmas[t]=gsl_stats_sd(gsl_matrix_const_column(bootstrap_corr,t+t0).vector.data,1,bootstrap_corr->size2); // gsl_stats_sd (const double data[], size_t stride, size_t n)
  }

  NonLinearFit::ModelExponentials modelExponentials(1,data);
  for(Count a=0;a<modelExponentials.nAddends;a++)
  {
    modelExponentials.c1s[a]=Real(t0); // A*exp(b*(t+t0))
    modelExponentials.c2s[a]=-Real(T)+Real(t0); // A*exp(-b*(-T+t+t0))

    modelExponentials.exponents[a]=Real(a+1)*initialExponent;
    modelExponentials.amplitudes[a]=1.-2.*int( data.ys[0]<0 );//initial amplitude should have same sign as data
  }


  NonLinearFit::calc(modelExponentials);
  for (unsigned int t=0; t<data.n; t++)
    cout<<t<<" "<<data.ys[t]<<" "<<modelExponentials.y(t)<<endl;

  gsl_matrix_memcpy (bootstrap_corr_relaxed, bootstrap_corr);
  for (unsigned int t=0; t<bootstrap_corr->size2; t++)
    gsl_vector_add_constant (&gsl_matrix_column (bootstrap_corr_relaxed, t).vector, modelExponentials.y(Real(t)-Real(t0))-gsl_matrix_get(bootstrap_corr,0,t));
}

*/
