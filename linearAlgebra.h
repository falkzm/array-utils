#pragma once
#include "array.h"



//////////////////////////////////////////Vektor-Ableitung
#define                                       Vector                 Array
#define                                       Vec                    Vector
#define                                       Matrix(ComponentSType) Array< Array< ComponentSType > >
#define                                       Mat                    Matrix
template<Count rank,class FundamentalT> using Tensor                 = typename nested< rank, Array, FundamentalT >::type;

#define for_vec(vec,i)   for(int i=0;i<vec.n();++i)
#define for_mat(mat,i,j) for_vec(mat,i) for(int j=0;j<mat.atFast(0).n();++j)
//////////////////////////////////////////Matrix-Ableitug




namespace LinearAlgebra // Matrix opr-func for Array are OK, but there should be a derived Type with corresponding opr to have the desired behaviour on default.
                    // There should be an main matrix type with explizit typecasts and a derived compatibility type with implicit cast behaviour (eg. for type-generic algorithm)
                    // The Array operations should still apply to LinearAlgebra types (for general tensor procedures, circumventing definition of unchanged opr as + ), but not the other way round.
                    // Types could be introduced as
                    //  - Mat & M
                    //  - Vec & V
                    //  - non-generic matrices (such as hermitian, eg.) should provide insert their own methods into the LinearAlgebra scope
                    //    & should be convertable into the generic versions to remove the need of reimplementing all methods.
                    // circumventing compatibility w Matrix & Vector issues. These are unlikly when Matrix and Vector do not provide implicit conversion to Array: unwanted * opr are then fetched by type safety.
                    // -some opr can be continued recursively for sub-vectors/matrices: an example is sd.
                    //  such opr can be either a) limited to lowest decomposition of tensor rank possible (->default?), b) maximal decomposition, c) setting the depth optionaly via template argument: either type or dimension e) as c) but dynamical
                    //  it would be optimal to provide dynamical an static control -> need 2 interfaces: template argument & default parameter.
{
  template<class SubSettable, class Scalar=fundamental_t<SubSettable> > inline static constexpr Count rank  =          Array<bool>::DimensionCount<SubSettable,Scalar>::value;
  template<Count rank,template<class ...> class Container, class... FundamentalT> using nested_t            = typename Array<bool>::nested< rank, Container, FundamentalT... >::type;
  template<Count rank,                                     class    FundamentalT> using tensor_t            = typename Array<bool>::nested< rank, Array    , FundamentalT    >::type;


  namespace Internal
  {
    template<bool recursive,class... K, std::enable_if_t< recursive>* =nullptr> inline decltype(auto) mult(K&&... k); // needed in order to control whether Array-valued fields K should be traeted as defined by operator* (->elementwise) or using the local linear algebra multiplication (=>set template argument recursive=true)
    template<bool recursive,class... K, std::enable_if_t<!recursive>* =nullptr> inline decltype(auto) mult(K&&... k); // LinearAlgebra::multiplication(.) can not be used directly bc in a namespace functions can not be invoked before their declaration. Using multiplication directly by other linearAlgebra::functions would require each varient of the overload set to be forward declared. Internal::mult - wrapper circumvent that need.

    template<bool recursive> inline decltype(auto) matrixAddScalarExpansion(auto&& matrixA,auto&& matrixB)
    {
      if constexpr(!recursive||!is_subsettable_v<decltype(matrixA)>)
        return matrixA+=matrixB;
      else
      {
        if(matrixA.n()==matrixB.n())
          return matrixA+=matrixB;
        if(matrixB.n()==1) // we assume compatible sizes or size 1 for scalars
        {
          const auto s=matrixB.atFast(0).atFast(0);
          for(Count i=0;i<matrixA.n();++i)
            matrixA[i][i]+=s;
        }
        else
        {
          const auto s=matrixA.atFast(0).atFast(0);
          matrixA=matrixB;
          for(Count i=0;i<matrixA.n();++i)
            matrixA[i][i]+=s;
        }
        return matrixA;
      }
    }

    template<bool recursive,class KA,class KB> inline static KA scalarProductScalarExpansion(const Array<KA> &vector, const KB &scalar) // used for rare blockwise-optimisation of the scalar product
    {
      __ARRAY_DEBUG__printV__("LinearAlgebra::scalarProductScalarExpansion",3);
      KA result=mult<recursive>(vector[0],scalar);
      for(Count i=1;i<vector.n();++i)
        result += mult<recursive>( vector[i], scalar );

      return result;
    }
  };


  template<class K> static Array< Array<remove_cr_t<K>> > inline createMatrix(const Count rowSCount,const Count columnSCount)
  {
    return Array<Array<remove_cr_t<K>>>( Array<Count>{rowSCount,columnSCount}, NULL );
  }


  template<class K> static void id(Array< Array<K> >& quadraticMatrix)
  {
    quadraticMatrix=0;
    for(Count i=0;i<quadraticMatrix.n();++i)
      quadraticMatrix[i][i]=1;
  }

  template<class K> static Array< Array<remove_cr_t<K>> > id(const Count rowSCount)
  {
    auto quadraticMatrix=createMatrix<remove_cr_t<K>>(rowSCount,rowSCount);
    id(quadraticMatrix);
    return quadraticMatrix;
  }


  template<class K> static Array<K> diag(const Array< Array <K> > matrix) // extracts diagonal
  {
    Array<K> vector(min(matrix.size(),matrix.atFast(0).n()),NULL);
    for(Count i=0;i<vector.n();++i)
      vector[i]=matrix[i][i];
    return vector;
  }

  template<class K> static Array< Array<K> > diag(const Array<K> vector, Count rowSCount=0, Count columnSCount=0) // creates diagonal matrix from vector
  {
    if( !rowSCount    )     rowSCount = vector.size();
    if( !columnSCount )  columnSCount = vector.size();

    Array< Array <K> > diagonalMatrix = ( createMatrix<K>( rowSCount, columnSCount )=0 );

    for(Count i=0;i<min(min(rowSCount,columnSCount),vector.size());++i)
      diagonalMatrix.atFast(i).atFast(i)=vector[i];
    return diagonalMatrix;
  }


  template<class K> static Array<K> conjugate(Array<K> a)
  {
    if constexpr(!is_complex_v<std::decay_t<fundamental_t<K>>>)
      for(Count i=0;i<a.n();++i)
        a[i]=conj(std::move(a[i]));
    return a;
  }
  template<class K, class=enable_if_notInheritsArray_t<K>> static complex<K> conjugate(complex<K> scalar)
  {  return std::conj(scalar);  }
  template<class K, class=enable_if_notInheritsArray_t<K>> static K          conjugate(K          scalar)
  {  return scalar; }


  template<class K> static Array< Array<K> > transpose(Array< Array<K> > matrix)
  {
    const Count&  rowSCount   =matrix.atFast(0).n(),
                  columnSCount=matrix.n();

    if( rowSCount!=columnSCount || matrix.n()>1 )
    {
      Array< Array<K> > resultMatrix=createMatrix<K>(rowSCount,columnSCount);

      for(Count i=0;i<rowSCount;++i)
        for(Count j=0;j<columnSCount;++j)
          resultMatrix.atFast(i).atFast(j)=std::as_const(matrix)[j][i];
      return resultMatrix;
    }

    for(Count i=0;i<rowSCount;++i) // in-place transposition
    {
      Array<K> &rowISPComponents=matrix[i];
      for(Count j=0;j<i;++j)
        std::swap(rowISPComponents[j],matrix[j][i]);
    }
    return matrix;
  }
  template<class K> static Array< Array<K> > transpose(const Array<K> vector)
  {
    Array< Array<K> > resultMatrix=createMatrix<K>(1,vector.n());

    for(Count i=0;i<vector.n();++i)
        resultMatrix.atFast(0).atFast(i)=vector[i];

    return resultMatrix;
  }
  template<class K,std::enable_if_t<!is_subsettable_v<K>>* = nullptr > K transpose(K k) //trivial transpose for scalar types (used by (blocK) matrix multipliction)
  {
    return k;
  }

  template<class P> static auto adjoint(P&& p)
  {
    return conjugate(transpose(std::forward<P>(p)));
  }


  template<bool recursive=false,class... K> decltype(auto) inline multiplication(K&&... k) // this might be used (recursively) by more the other (more specialised) multiplication algorithm and serves as trivial recursion anchor.
  {  return (std::forward<K>(k) * ...);  }


  template<bool recursive=false,bool scalarExpansion=false,class KA,class KB> static KA scalarProduct(const Array<KA> vectorA, const Array< KB > vectorB) // recursive==true enables elementwise multipliction be treated as linearAlgebra multiplication, scalarReduction==false (default) turns of a feature for blockwise multiplications: expansion of single-element-vectors to full-sized vectors since such a optimisation might be unnecessary for a rank-1-operation
  {
    __ARRAY_DEBUG__printV__("LinearAlgebra::scalarProduct",3);
    if constexpr (recursive&&scalarExpansion)
    {
      if(!vectorA.n()||!vectorB.n())
        return KA();
      if(vectorA.n()==1&&vectorB.n()!=1)
        return Internal::scalarProductScalarExpansion<recursive>(std::move(vectorB),vectorA[0]);
      if(vectorB.n()==1&&vectorA.n()!=1)
        return Internal::scalarProductScalarExpansion<recursive>(std::move(vectorA),vectorB[0]);
    }
    else
      __ARRAY_DEBUG__assert__(vectorA.n()<vectorB.n(),"ERROR: LinearAlgebra::scalarProduct: vectorA.n=="<<vectorA.n()<<" but vectorB.n=="<<vectorB.n());
    KA result=Internal::mult<recursive>(vectorB[0],vectorA[0]);

    for(Count i=1;i<vectorA.n();++i)
      result += Internal::mult<recursive>( vectorB[i], vectorA[i] );

    return result;
  }
  template<bool recursive=false,class KA,class KB> static KA multiplication(const Array<KA> vectorA, const Array< KB > vectorB)
  {  return scalarProduct<recursive>(std::move(vectorA),std::move(vectorB));  }


  template<class RUser=void, bool conjA=true, bool recursive=false, class KA,class KB, keyClasses(Key,R,RUser,KA)> static R generalisedScalarProduct(const KA a,const KB b,Key key=0)
  {
    if constexpr (is_same_v<remove_cr_t<KA>,remove_cr_t<R>>)
      if constexpr (conjA)
        return conjugate(a)*b;
      else
        return a*b;
    else
    {
      Count i=min(a.size(),b.size())-1;
      if constexpr (is_same_v<remove_cr_t<decltype(a[0])>,remove_cr_t<R>>)
        if constexpr ( !conjA )
        {
          R r{ Internal::mult<recursive>(a.template at<R>(a[i],key),a.template at<R>(b[i],key)) };
          for( ; i-->0; )
            r += Internal::mult<recursive>(a.at(a[i],key),b.at(b[i],key));
          return r;
        }
        else
        {
          R r{ Internal::mult<recursive>(conjugate(a.template at<R>(a[i],key)),a.template at<R>(b[i],key)) };
          for( ; i-->0; )
            r += Internal::mult<recursive>(conjugate(a.template at<R>(a[i],key)),b.template at<R>(b[i],key));
          return r;
        }
      else
      {
        R r{ generalisedScalarProduct<R,conjA,recursive>(a[i],b[i],key) };
        for( ; i-->0; )
          r += generalisedScalarProduct<R,conjA,recursive>(a[i],b[i],key);
        return r;
      }
    }
  }


  template<bool recursive=false,class KA,class KB> static Array< Array<KA> > matrixScalarMultiplication(Array< Array<KA> > matrix, const KB scalar) //needed by matrixMultiplication for blockwise optimisation for trivial/scalar blocks
  {
    __ARRAY_DEBUG__printV__("LinearAlgebra::matrixScalarMultiplication",3);
    const Count&  rowSCount=matrix.n(), columnSCount=matrix[0].n();
    for(Count i=0;i<rowSCount;++i)
    {
      Array<KA> &matrixSRow=matrix.atFast(i); // atFast bco matrix[0].n() is a non-const access => Array::derefer()
      for(Count j=0;j<columnSCount;++j)
        matrixSRow[j] = Internal::mult<recursive>(std::move(matrixSRow[j]),scalar);
    }
    return matrix;
  }


  template<bool recursive=false,class KA,class KB> static Array< Array<KA> > matrixMultiplication(Array< Array<KA> > matrixA, const Array< Array<KB> > matrixB) // setting recursive true makes element-multiplications being evaluated with linearAlgebra's multiplication routines. Hence, Array<Array<.>> invokes again matraix multiplication for the elements, eg.
                                                                                                                                                                // it also enables checks for one factor beeing zero-sized or scalar-like in order to enable blockwise multiplications with certain elements being simplified-
  {
    __ARRAY_DEBUG__printV__("LinearAlgebra::matrixMultiplication",3);
    if constexpr (recursive)
    {
    //cout<<"m1 "<<MathematicaArray(matrixA.getSize())<<" "<<MathematicaArray(matrixB.getSize())<<" "<<MathematicaArray(matrixA)<<endl<<flush;
      //cout<<"m1\n"<<flush;
      if(!matrixA.n()||!matrixB.n())
        return decltype(matrixA)();
      //cout<<"m2\n"<<flush;
      if(matrixA.atFast(0).n()==1)
        return matrixScalarMultiplication<recursive>(std::move(matrixB),matrixA.atFast(0).atFast(0));
      //cout<<"m3\n"<<flush;
      if(matrixB.n()==1)
        return matrixScalarMultiplication<recursive>(std::move(matrixA),matrixB.atFast(0).atFast(0));
    }
      //cout<<"m4\n"<<flush;
    const Count&  rowSCount   =matrixA.n(), // writing result into matrix A
                  columnSCount=matrixB.atFast(0).n();

    if(matrixA.n()>1||columnSCount!=matrixA.atFast(0).n()) // A referenced multiple times or B is not quadratic => generic multiplication
    {
      //cout<<"m5\n"<<flush;
      auto resultMatrix=createMatrix<KA>(rowSCount,columnSCount);

      for(Count i=0;i<rowSCount;++i)
      {
        /**/ if constexpr(rank<KA> > 0) cout<<"m5.1 "<<i<<endl;
        const auto& matrixASRow=static_cast<const decltype(matrixA)>(matrixA)[i];
        for(Count j=0;j<columnSCount;++j)
        {
          auto &resultMatrixSComponentIJ=resultMatrix.atFast(i).atFast(j);

          resultMatrixSComponentIJ = Internal::mult<recursive>( matrixASRow[0], matrixB[0][j] );
          for(Count k=1;k<matrixB.n();++k)
            Internal::matrixAddScalarExpansion<recursive>( resultMatrixSComponentIJ, Internal::mult<recursive>( matrixASRow[k], matrixB[k][j] ) );
      //cout<<"m6\n"<<flush;
        }
      //cout<<"m7 "<<i<<endl<<flush;
      }

      //cout<<"m8\n"<<flush;
      return resultMatrix;
    }
    else // this branch potentially reuses memory of matrixA if A is referenced only once (here).
    {
      //cout<<"m9\n"<<flush;
      KA v[columnSCount];
      for(Count i=0;i<rowSCount;++i)
      {
        auto &matrixASRow=matrixA.atFast(i);
        for(Count j=0;j<columnSCount;++j)
        {
          KA &v_j=v[j];

          v_j=0;
          for(Count k=0;k<columnSCount;++k)
            Internal::matrixAddScalarExpansion<recursive>( v_j, Internal::mult<recursive>( static_cast<const decltype(matrixASRow)>(matrixASRow)[k], matrixB[k][j] ) );
        }
        for(Count j=0;j<columnSCount;++j) // writing result into matrix A
          matrixASRow[j]=v[j];
      }

      //cout<<"m10\n"<<flush;
      return matrixA;
    }
  }
  template<bool recursive=false,class KA,class KB> static Array< Array<KA> > multiplication(const Array< Array<KA> > matrixA, const Array< Array<KB> > matrixB)
  {  return matrixMultiplication<recursive>(std::move(matrixA),std::move(matrixB));  }


  template<bool recursive=false,class KA,class KB> static Array<KA> matrixVectorMultiplication(const Array< Array<KA> > matrix, const Array<KB> vector)
  {
    __ARRAY_DEBUG__printV__("LinearAlgebra::matrixVectorMultiplication",3);
    if constexpr (recursive)
    {
      if(!matrix.n()||vector.n())
        return Array<KA>();
      if(matrix.n()==1&&matrix.atFast(0).n()==1&&vector.n()!=1)
      {
        Array<KA> resultVector(vector.n(),NULL);
        const KA &m=matrix[0][0];
        for(Count i=0;i<vector.n();++i)
          resultVector.atFast(i)=Internal::mult                        <recursive>(vector[i],m);
        return resultVector;
      }
      if(vector.n()==1&&(matrix.n()!=1||matrix.atFast(0).n()!=1))
      {
        Array<KA> resultVector(matrix.n(),NULL);
        const KA &v=vector[0];
        for(Count i=0;i<matrix.n();++i)
          resultVector.atFast(i)=Internal::scalarProductScalarExpansion<recursive>(matrix[i],v);
        return resultVector;
      }
    }
    const Count& rowSCount=matrix.n();

    Array<KA> resultVector(rowSCount,NULL);

    for(Count i=0;i<rowSCount;++i)
      resultVector.atFast(i)=scalarProduct<recursive>(matrix[i],vector);

    return resultVector;
  }
  template<bool recursive=false,class KA,class KB> static Array<KA> multiplication(const Array< Array<KA> > matrix,const Array<KB> vector)
  {  return matrixVectorMultiplication<recursive>(std::move(matrix),std::move(vector));  }


  template<class K> static K trace(const Array< Array<K> > quadraticMatrix)
  {
    K resultK(quadraticMatrix[0][0]);

    for(Count i=1;i<quadraticMatrix.n();++i)
      resultK+=quadraticMatrix[i][i];

    return resultK;
  }


  template<bool recursive=false,class M> static M matrixPow(const M matrix,Count n) // Only for n in Z>=0.
                                                               // O(k^3*log(n)) for k x k matrix. NOTE: if the matrix has full rank & n>>0 another method should be uses: diagonalise A=Q^T D Q => A^n=Q^T D^n Q which costs only ~O(k^3) for the diagonalisation
                                                               // Here we use the recursive approach pow(A,n)=pow(A,n/2)^2 <=>n even and pow(A,n)=A*pow(A,n-1) <=>n odd .
  {
    //cout<<"mp"<<endl;
    if(n<=1)
      if(n==1)
        return matrix;
      else
        return id<decltype(matrix[0][0])>(matrix.n());
    else
      if(n%2)
        return matrixMultiplication<recursive>(matrixPow<recursive>(matrix,n-1),matrix);
      else
      {
        auto b=matrixPow<recursive>(matrix,n/2);
        //cout<<"\n\n\n < "<<b<<" > "<<endl;
        return matrixMultiplication<recursive>(b,b);
      }
  }
  template<class K,class N> static auto pow(K&& k,N&& n)
  {
    if constexpr( is_subsettable_v<K> )
      return matrixPow(std::forward<K>(k),std::forward<N>(n));
    else
      return std::pow (std::forward<K>(k),std::forward<N>(n));
  }


  template<bool highPrecision=true,class K> static K determinant(Array<Array<K>> m,K eps=1E-8) // generic O(n^3) det computation ( using det(A)=det(A+b*e_j*A_i^T) (invariance under adding a multiple of arow(col) to other rows(cols)) to achieve upper triangle form. )
  {
    const Count n=m.n();
    K det=1;
    //det=ArrayTools::max<K>(m);m/=det; det=pow(det,n);
    //if constexpr(highPrecision) // we do this even if no highPrecision is needed: bco derefer() will be called for each column and element therein, we could use atFast everywhere below.
      for(Count repeat=1;repeat>0;repeat--) // balancing //n=1 w final transpose lead to best results in a testcase. might have been random.
      {
        for(Count i=0; i<n; ++i)
        {
          auto& m_i = m[i];
          K t=0;
            for(Count j=0; j<n; ++j)
               t+=abs(std::as_const(m_i)[j]);
            t/=n;
            det*=t;
          //cout<<"     "<<t<<" -- "<<m_i<<endl;
          m_i/=t;
        }
        if(repeat>1) // seem to help even wo repeatition
          m = transpose(std::move(m));
      }
    for(Count i=0; i<m.n(); ++i)
    {
      const auto& m_i = std::as_const(m)[i];
      if constexpr (highPrecision) //search absolut element in ther lower right (n-i)x(n-i) matrix closest to 1 & swap rows and columns st. it becomes m_ii. Below, the i-th/m_ii will be used to zero colum i's elements below the i-th one (->upper triangle form). 1/_mii~1 will not decrease the spread of magnitudes within the matrix. This helps for matrices with det close to zero: the rows or columns are then almost parralel and the little deviations determine the determinant.
      {                            //Hence those deviations should not be swallowed by adding/subtr. large and tiny absolut numeric values.
        Count jBest=i,kBest=i; K vBest=abs(1-abs(m_i[i]));
        for(Count j=i; j<n; ++j)
          for(Count k=i; k<n; ++k)
        {
          K v=abs(1-abs(m.atF(j).atF(k)));
          if( v>vBest )
            { vBest=v; jBest=j; kBest=k; }
        }
        std::swap(m.atF(i),m.atF(jBest));
        if(kBest!=i)
          for(Count j=i; j<n; ++j)
            std::swap(m.atF(j).atF(kBest),m.atF(j).atF(i));

        if((jBest+kBest)%2) det=-det;
      }
      else
        if(abs(m_i.atF(i))<eps) [[unlikely]] // too small diagonal element -> search for one closer to 1 and swap rows.
        {
          Count jBest=i; K vBest=abs(1-abs(m_i.atF(i)));
          for(Count j=i+1; j<n; ++j)
          {
            K v=abs(abs(1-m.atF(j).atF(i)));
            if( v<vBest )
              { vBest=v; jBest=j; }
          }
          std::swap(m.atF(i),m.atF(jBest));
          if((jBest-i)%2) det=-det;
        }

      // This is the actual elemmination to achieve the upper-triangle form
      const K& a = m_i.atF(i);
      det*=a;
      if(det==0) // catches a==0
        return 0;
      //cout<<"  m_(i=="<<i<<") "<<m[i]<<endl<<"  d="<<det<<endl;
      for(Count j=i+1; j<n; ++j)
      {
        auto&    m_j = m.atF(j);
        const K  b   = std::as_const(m_j)[i]/a;
        K t=0;
        for(Count k=i+1; k<n; ++k)
          t+=abs( m_j.atF(k) -= b*m_i.atF(k) );
        if constexpr (highPrecision) // renormalise remaining row to 1
          if((t/=n)!=0)
        {
          for(Count k=i+1; k<n; ++k)
            m_j.atF(k)/=t;
          det*=t;
        }
      }
    }
    return det;
  }


  class Norm
    // todo: norms
    // move basic-tensor operations from DataAnalysis & LinearAlgebra like conjugate, abs, re, imag, etc. into a comman base - like ArrayTools. => prevent double implementatins
    // - canonical from SP
    // - std/froben ||.||_F
    // - eucledian/least squares/std ||.||_2
    // - maximum element ||.||_me
    // - maximum orientation ||A||_mo:=max{ A v | \fa v in VectorSpace \w ||v||_2==1 }
  {
  public:

    template<class K> static inline K norm(K k) // generalisation of std::norm.  ATTENTION: it lacks the square root in comparision to the mathematical standard norm
    {
      return mult(conjugate(k),k);
    }


    template<class K> static auto frobenius(const Matrix(K) m)
    {
      decltype(norm(m[0][0])) r=0;

      for(Count i=0;i<m.n();++i)
        for(Count j=0;j<m[i].n();++j)
          r+=norm(m[i][j]);
      return sqrt(r);
    }


    template<class K> static auto canonical(const Vector<K> v)//scalar product should be redefinable -> automatically achieved, when * becomes specialised for Vector types
    {
      return scalarProduct(conj(v),v);
    }


    template<class K> static auto eucledian(const Vector<K> v)
    {
      auto r=norm(v[0]);
      for(Count i=1;i<v.n();++i)
        r+=norm(v[i]);
      return sqrt(r);
    }
  };


  // definition of Internal::mult making use of the overload-set of multiplication functions
  template<bool recursive,class... K, std::enable_if_t<!recursive>* =nullptr> inline decltype(auto) Internal::mult(K&&... k)   {  return (std::forward<K>(k) * ...);  }
  template<bool recursive,class... K, std::enable_if_t< recursive>* =nullptr> inline decltype(auto) Internal::mult(K&&... k)   {  return multiplication<true>(std::forward<K>(k)...);  }

  // short hands ////////////////////////////////////////////////////////////////////////////////////////////////////

  template <typename    P> inline static auto conj  (P&&    p)  {  return conjugate      (std::forward<P>(p));  }
  template <typename... P> inline static auto sp    (P&&... p)  {  return scalarProduct  (std::forward<P>(p)...);  }
  template <typename... P> inline static auto mult  (P&&... p)  {  return multiplication (std::forward<P>(p)...);  }
  template <typename    P> inline static auto tr    (P&&    p)  {  return trace          (std::forward<P>(p));  }
  template <typename    P> inline static auto transp(P&&    p)  {  return transpose      (std::forward<P>(p));  }
  template <typename    P> inline static auto trp   (P&&    p)  {  return transpose      (std::forward<P>(p));  }
  template <typename    P> inline static auto adj   (P&&    p)  {  return adjoint        (std::forward<P>(p));  }
  template <typename    P> inline static auto det   (P&&    p)  {  return determinant    (std::forward<P>(p));  }


  static void check()
  {
    using namespace std;
    typedef double K;

    cout<<"LinearAlgebra::check(): sanity check using a minimal example...\n\n"<<flush;

    Matrix(K) m={Array<K>{1,2},Array<K>{3,4}};

    cout<<"Matrix m: "<<m<<endl;

    cout<<"id(2)*m = "<<mult(id<K>(2),m)<<endl;
    cout<<"m*id(2) = "<<mult(m,id<K>(2))<<endl;
    cout<<"transp(m) = "<<transp(m)<<endl;

    Vector<K> v{10,20};

    cout<<"Vector v: "<<v;
    cout<<"\ntransp(v)*v = "<<mult(transp(v),v);
    cout<<"\n<v,v> = "<<sp(v,v);
    cout<<"\nid(2)*v = "<<mult(id<K>(2),v);
    cout<<"\ntransp(v)*id(2)*v ="<<mult(transp(v),mult(id<K>(2),v));
    cout<<"\nm*v ="<<mult(m,v);
    cout<<"\ntransp(v)*m*v ="<<mult(transp(v),mult(m,v))<<endl;

    cout<<"Repeatation with temporary m ...\n";
      cout<<"  id(2)*m = "<<mult(id<K>(2),Array<Array<K>>{Array<K>{1,2},Array<K>{3,4}})<<endl;
      cout<<"  m*id(2) = "<<mult(Array<Array<K>>{Array<K>{1,2},Array<K>{3,4}},id<K>(2))<<endl;
      cout<<"  transp(m) = "<<transp(Array<Array<K>>{Array<K>{1,2},Array<K>{3,4}})<<endl;
      cout<<"\n  m*v ="<<mult(Array<Array<K>>{Array<K>{1,2},Array<K>{3,4}},v);
      cout<<"\n  transp(v)*m*v ="<<mult(transp(v),mult(Array<Array<K>>{Array<K>{1,2},Array<K>{3,4}},v))<<endl;

    cout<<"Check of non-quadratic matrices.\n";
    Matrix(K) n={Array<K>{1,2},Array<K>{3,4},Array<K>{5,6}};

    cout<<"  Matrix n: "<<m<<endl;
    Vector<K> w{10,20,30};

    cout<<"  Vector w: "<<w<<endl;

    cout<<"  id(3)*n = "<<mult(id<K>(3),n)<<endl;
    cout<<"  n*id(2) = "<<mult(n,id<K>(2))<<endl;
    cout<<"  transp(n) = "<<transp(n)<<endl;
    cout<<"\n  n*v ="<<mult(n,v);
    cout<<"\n  transp(w)*n*v ="<<mult(transp(w),mult(n,v))<<endl;
  }
};
