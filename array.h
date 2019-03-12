#pragma once
#include <complex>
#include <iostream>



template<class ComponentSType> class Array; //forward declaration needed by arrayTraits.h
#include "arrayTraits.h"
template<class T> inline static constexpr bool is_complex_v = boost::is_complex<T>::value;

template<class T> using remove_cr_t = std::remove_const_t<std::remove_reference_t<T>>;

// NOTES
//
// constness is crucial for speed and memory usage: any non-const reference to the components lead to a mem copy (Array::derefer() call) if multiple containers (Array) exists
// + this is especially true for Array::opr[]
// copy is usually avoided for parameter passing and temporary Arrays within formulas.
//
// the data types should behave as algebras opr(ie. *,+) must preserve the structure
// using std opr for non-std types enables usage of template algorithm with these algebras
// own algebra enable special features as precision and error tracking ...
//
// The internal memory handling can be debugged by defining
//     #define __ARRAY_DEBUG__
//
//   and setting the desired verbosity via the global variable
//     unsigned lineareAlgebraSBloedeKommentare=0;  //verbosity
//
//   to any number >0. Use this variable to enable debugging locally ( in order to avoid unnecessary & potentially excessive debugging ) and set the level of verbosity.
//   The definitions can be found below this initial documentation.
//
//
// - Not initialized arrays:
// -- operators are not defined for not initizialised array, however resize(), size(), ==, != are legal and the member pData will point to zero.
// -- initialisation of multidimensional structure like Arra< Array < .. > > can be done by calling resize(...) with an index array that specifies the size in the desired number of dimensions.
// -- an array x<T> can be initialised to 0 (if the underlying type can hold 0): this results in an array of size 1 and value x[0] == static_cast<T>(x) == 0
// -- an array x<T> can be initialised to size 0. Comparison operations are defined but any access x[0], static_cast<T>(x) is undefined
//
// TODO
//
// simple modernisation?
//  - have a look into boost, eigen (& other advanced Tensor/LinAlg API's) to grasp potentially useful further concepts.
//  - use && per default
//  - dynamic memory model as 2nd, seldomly-to-use scenario, preprocessor optimisations
//  - swap & move
//  - provide external function templates for getFundamental, size, n, [], ... : template<class T> f(T t) { if constexpr (is_array<T>::value) return t.f(t); else return [TrivialScalar case]; }
//    This enables templates action on Scalars as well as on Array types
//  - Avoidance of inhertance for manipulated ('derived') Array types <=> avoidance of virtual operators
//    => use #include modules /usage of namespaces instead?
//    [] etc might become external opr (non-member) st. overload resolution takes place
//    if [] must be a member it might be altered by traits: it invokes a non-member func
//      the downside is access to private data from non-member: 'friend' statements help
//      also: some opr must be members, however a generic
//      #define GenericOpr(Container,Element) Element& Container::[](int i)
//              { return index(i); }
//      helps, when the template non-member templ<class C<>,class E,enable_if_inherits<Array,C> > E& index(C<E>* p,int i) { ret ... }
//      contains the specific implemenetation of Array::[]. Usage of templates ensures type forwarding for internal operator calls (in contrast to index(Array<E>*,int) w upcasts)
//      Usage of concepts/enable_if_inherits qualifies the 'more specialised' overload soluation.
//
//  -- basic opr may use swap to select useCount=1-Parameter for return-val reusage (currently eg. opr+(Array,const Array): there is a bet that the 1fst param might be reused)
//  - if all elementwise opr would return a temporary obj and overloads for these obj, than there could
//  -- an evaluation when an access occurs
//  -- skip of evaluation at ==, <, > level via symbolical evaluation
//  --- here 0 & id matrices might become handy
//  +++ Support for algebras, group representations?
//  +++ Note: is becoming to complicated. What about existing systems? Sage?
//  - array should be selfexpanding for simple support of template algrorithem \w such things as T x=0; x+a; where T is an array an a an T-parameter \w a.n()>1
//  -- in selfassignment opr (*=, +=, ...) there should be a derefer_resize_swap call.
//  ++ this might indeed be one func to keep the opr themself simple.
//  ++ The +, - , / external opr should mostly inherit the behaviour from calling selfassigning opr.
//  - opr?(T1,T2) using selfassigning opr may fail if complexity(T1.getFundamentalType)!=complexity(T2.getFundamentalType)
//    eg., T1=Array<double>, T2=complex(double) => must be promoted to Array<complex<double>>.
//     exitence of conversion opr must select a preferable solution or, if no rank exists, take the first type.
//  + only if T1==T2 swap is possible
//  - adressing of type resolutions
//  -- a lot of std algorithm (mean, sum) can be made to resolve the arrays up to template type paramter T
//     by using AT::map (or its )
//  --- AT::map may act on rvalues for abitrary resolution depths:
//  ----> a) forced forward of referenceType when going [] b) wrap all calls of map such that paramters become lvalues (prefered)
//  -- Array< ... <Array<T> ... > representation via typecast as Array<... Scalar<Array<T> ...> ?
//  --- allows operators ( calls without option to specify template arg ) to act up to a specified resolution depths
//  ----- static_cast< Array<Array<Scalar< decltype(a[0][0]) >>>(a)*b;
//  ----- a.scalarView<2>*b
//  ----- a.scalarView<T>*b \w T:=delctype(a[0][0])
//  -- support for getSize, setSize(indeVector), at, [IndexVec]
//  -- swap at<Depth> to at<T> and provide adapter depth<Array,T>::value
//
//  - simple remapping view?
//  -- eg. shift index, matrix.view -> specialisations
//  -- should not copy for easy cases
//  -- shift, intervall, n-dim <-> m.dim
//  --- give typecast-like representation for 'basic' transformations: Vec<T> v; static_cast< Matrix<T,R=10> >(v) and dynamically: Matrix<T>(v,R=10)
//
//
// in that sense: array should be a very basic type (field), thus
// opr should, as far as possible, be defined elementwise
// -> * should be moved from sp to elemtwise behaviour of *=
// -.- mat/vec mult should be done via explicit type convertisation to mats/vecs
// simplification to Scalar should be resolved via static_cast to the first element

// define an abstract interface class ArrayInterface
// + derive Array and ArrayViews
// + ArrayView: it takes &Data of a (nested) Array type and a view method (type) defining the behaviour of [].
// + external entities should use the abstract interface type in order to support all derived classes

// Vectors and Matrices should be additionaly implemented by deriving Array and replacing the * Opr
// -> a huge class of Algebras could be than realised by Matrix representation
// - initialisation, = with Scalars s should be transformed into id*s -> 1, 0 would behave as neutral elements of +,- in template algorithm
// vec and mat class should have a selection std norms - especially for conversion to Scalar
// + and a default func norm() that uses the norm specified in a public member var -> for each algorithm the type of norm can be selected in advance
// + external opr should be derived by including "template" headers w type-to-be-specified by #define #undef
//   some opr could be grouped: opr(Scalar,Array) & opr(Array,Scalar) into one file

// In order to allow for easly derived Representations, acting like references, but with modified operator[]
// - Introduce generic Reference Array: it should reference-bind to another array
// - the modified operator[] should used operator[](.) of the referenced array => constness of data via derefer calls is preserved automatically
// - external and member operators should be converted into full-template with is_array-trait use: template<class ComponentSType> f(Array<ComponentSType>&) ~> template<class ArrayT, enable_if_is_array<ArrayT>::type=0> f(ArrayT&)
//   in order to support 'derived' Array types
// - how to forward non-operator[] operators into the 'derived' Array type while making them use the modified operator[]?
// - operators should use to operator[] with const-casts for read access in order for automatic derefer-calls
// <- Tests with multiple-Container layers shows g++'s and clang's ability to optimise the such an interface away.
//
  // [Why not using virtual operators]
  //
  // virtual member functions often lead to simpler expressions for modifable container (via inheritance) in contrast to a template operator/concept interface making use of an is_Array treat (or is_myDerivedType.)
  // The reason for this is the need of a conditional jumb in order to select the appropiate virtual member operator from the v-table.
  // V-tables might be realised very similar to function pointers.
  // The argument is, that in general the conditional jump can not be optimized away, since:
  // - the adress of the operator could be taken and stored in dynamically allocated objects.
  // - the code could contain a comparison of the stored function pointer with the operators address.
  // -> there can not be multiple implementations of the operator
  // - even if the container consist of a pure abstract interface class and one class with the implementation of virtual operators, such that the v-table becomes trivial, there could be always non-trivial container, which must be processable by the same function,
  //   dealing with the abstract interface class.
  // - there could be only a 2nd implementation for inlineable situations where the argument containers are known to be of the trivial kind.
  // - this is not the case for any dynamically allocated container.




// there might be abstract algebras for speedup and more theoretical problems, eg. a List and RSV representation of numbers for quantisation

#ifdef DEBUG
  #define __ARRAY_DEBUG__
#endif
unsigned lineareAlgebraSBloedeKommentare=0; //verbosity


#ifdef __ARRAY_DEBUG__
  #define __ARRAY_DEBUG__print__(TextStream)                     {  if(lineareAlgebraSBloedeKommentare)                     std::cout<<TextStream<<std::endl<<std::flush;  }
  #define __ARRAY_DEBUG__printV__(TextStream,requiredVerbosity)  {  if(lineareAlgebraSBloedeKommentare>=requiredVerbosity)  std::cout<<TextStream<<std::endl<<std::flush;  }
  #define __ARRAY_DEBUG__assert__(condition,TextStream)          {  if(condition)                                           std::cout<<TextStream<<std::endl<<std::flush;  }
#else
  #define __ARRAY_DEBUG__print__(TextStream)
  #define __ARRAY_DEBUG__printV__(TextStream,requiredVerbosity)
  #define __ARRAY_DEBUG__assert__(condition,TextStream)
#endif



// macros.
#define for_Array(array,index) for(Count index=0; index<array.size(); ++index)
//suggestion: #define for_A for_Array

// reservation for an for_each version for arrays: should provide numa-local parallelisation
// - note: lambda [&]... can be used to capture the complete scope and forward into sub routines
// + direct function calls allow for limited overload resolution, eg for_each(Range,processorN) where Range could be a type  with construction Range(Count i0,Count i1,d1=1)
//   Preferable over for_each_range bco possible template scenarios can employ both implementations.
// - should use STL parallelism
// - sugeestion: provide a version which allows to specify nr of processes. #define for_each_array_n
//   The experimental/theoretical Array implementation with DataModels might house a processing plan, too.
#define for_each_array for



template<class ComponentSType> class Array
{
public:
  typedef ComponentSType ArraySType;// ComponentSType could not be used directly; only for Array< T >::Type to recover SubComponentSType if Array< T > is given as single expression (eg. template type)a
                                    // BUT decltype(<array[0]>) would do the same (?) - no: also defined T!=Array< > but an other Type that has an operator []

  class Data
  {
  public:
    ComponentSType* __restrict__ pComponents;
    Count           componentSCount;
    mutable Count   useCount;


    Data() noexcept
      : pComponents(NULL),useCount(0)
    {}

    Data(const Count componentSCount) noexcept
      : pComponents( new ComponentSType[componentSCount] ),componentSCount( componentSCount ),useCount(0)
    {
      #ifdef __ARRAY_DEBUG
       if(componentSCount==0)
         std::cout<<"WARNING: requesting Array of size 0"<<std::endl<<std::flush;
       if(pComponents==NULL)
         std::cout<<"ERROR: Array storage could not be allocated."<<std::endl<<std::flush;
       #endif
    }
    // The following constructor would allow for direct & individual initialization using std::initializer_list. But: I do not know a good method to pass the list to new ComponentSType[n] {list} since the syntax requires(?) the {}-brackets. But {list} has type initializer_list<initializer_list>(list).
    // The circumvention below uses a second, non-allocation call of new, however, this unusal form might be inefficient directly, or prevent standard optimisation
    //
    //template<class Scalar> Data(std::initializer_list<Scalar> list) noexcept
    //  : pComponents( static_cast<ComponentSType*>(::operator new[](list.size()*sizeof(ComponentSType))) ), //pComponents( new ComponentSType[list.size()]{}{list} ),
    //     componentSCount( list.size() ),useCount(0)
    //{
    //  auto l=list.begin();
    //  for(auto i=0; i<componentSCount; l++, ++i)
    //    ::new (pComponents+i) ComponentSType(*l);
    //}
    //template<class... Scalars> Data(Salars &&... list)
    //  : pComponents( new ComponentSType[sizeof...(list)]{std::forward<Scalars>(list)...} ),componentSCount( sizeof...(list) ),useCount(0)
    //{}

    ~Data() noexcept
    {
      delete[] pComponents;
    }
  };


  Data* pData;

protected:

  inline void useData(Data* _pData)
  {
    if(_pData!=NULL)
    {
      __ARRAY_DEBUG__print__("Array::useData with useCount="<<_pData->useCount);

      ++(_pData->useCount);
    }
    pData=_pData;
  }

  inline void finishDataSUse()
  {
    if(pData!=NULL)
    {
       __ARRAY_DEBUG__print__("Array::finishDataSUse with useCount="<<pData->useCount);

      if(pData->useCount<=1)
        delete pData;
      else
        --(pData->useCount);
    }
  }


  template< class ScalarUser=void, Count i=0, class ArrayT, class Scalar=first_nonVoid_t<ScalarUser,fundamental_t<ArrayT>> >
    static void resize(ArrayT& array,const Array<Count>& sizeArray,bool keepComponents)
  {
    __ARRAY_DEBUG__print__("res"<<i);

    if constexpr(!std::is_same_v<ArrayT,remove_cr_t<Scalar>>)
    {
      if(i<sizeArray.pData->componentSCount)
      {
       __ARRAY_DEBUG__print__("res"<<i<<" -> "<<sizeArray[i]);
        array.resize(sizeArray[i],keepComponents);
        for(Count j=0;j<array.pData->componentSCount;++j)
          resize<Scalar,i+1,remove_cr_t<decltype(array[0])>,Scalar>(array[j],sizeArray,keepComponents);
      }
    }
  }

  template< class ScalarUser=void, Count i=0, class ScalarIndex, class ArrayT, class Scalar=first_nonVoid_t<ScalarUser,fundamental_t<ArrayT>> >
    static Array<Count>  getSize(const ArrayT& array,const Array<ScalarIndex>& indexArray)
  {
    if constexpr(std::is_same_v<ArrayT,remove_cr_t<Scalar>>)
      return Array<Count>(i,NULL);
    else
    {
      __ARRAY_DEBUG__print__("getSize("<<i);

      if(array.pData==NULL||array.pData->componentSCount==0)
      {
        __ARRAY_DEBUG__print__("NULL");
        Array<Count> sizeArray(i+1,NULL);
          sizeArray[i]=0;
        return sizeArray;
      }
      //cout<<( indexArray.pData==NULL||i>=indexArray.pData->componentSCount ? 0 : indexArray[i] )<<flush;
      Array<Count> sizeArray = getSize<void,i+1,ScalarIndex,remove_cr_t<decltype(array[0])>,Scalar> //defining all template param should reduce compile time.
                                 (array[ ( indexArray.pData==NULL||i>=indexArray.pData->componentSCount ? 0 : indexArray[i] ) ],indexArray);
        sizeArray[i]=array.pData->componentSCount;
      //cout<<array.pData->componentSCount<<endl<<flush;
      return sizeArray;
    }
  }


  template<bool noDerefer,class SubSettable>inline static decltype(auto) basicAt(SubSettable &&array, auto&& i)
  {
    if constexpr(inheritsArray_v<SubSettable>&&(std::is_const_v<SubSettable>||noDerefer))
    {
      __ARRAY_DEBUG__printV__("Array::[](Count="<<i<<")",2);
      __ARRAY_DEBUG__assert__(i<0||i>=array.size(),"ERROR: Array::[]("<<i<<") (basicAt) but size=="<<array.size());

      //array.pData->pComponents[i]; // forwardMember makes not much sense for Array<.> types due to the smartpointer-feature. However, we propagate the rvalue attribute as expected.
      return forwardMember<SubSettable>(array.pData->pComponents[i]); // forwardMember makes not much sense for Array<.> types due to the smartpointer-feature. However, we propagate the rvalue attribute as expected.
    }
    else
      //return array[i];
      return forwardMember<SubSettable>(array[i]);
  }

public:

  #include "arrayTraits.h"


  template<class T, class=typename enable_if_notInheritsArray<T>::type>
                    static T&& getFundamentalType(T&&);  // use these functions only to deduce the fundamental type of an Array<..> via decltype( getFundamentalType(array) )
                                                         // new version: use decltype(opr[]) to deduce cv and reference qualifiers -> more centralised

  template<class T> static auto getFundamentalType(      Array<T>&  array) -> decltype( getFundamentalType(array[0]) );
  template<class T> static auto getFundamentalType(const Array<T>&  array) -> decltype( getFundamentalType(array[0]) );

  // and some && functions for rvalue references: defining a wrapper function using static_cast<Array<T> &>(array) will result in an infinite self-loop (probably due to return auto and maybe bco T& operator[]() resulting in a prvalue/xvalue)
  template<class T> static auto getFundamentalType(      Array<T>&& array) -> decltype( getFundamentalType(array[0]) );
  template<class T> static auto getFundamentalType(const Array<T>&& array) -> decltype( getFundamentalType(array[0]) );

                    static auto getFundamentalType(ComponentSType c=std::declval<ComponentSType>()) -> decltype( getFundamentalType(c) );



protected:

  template<int lineBreakDim=0,int d=0,class Scalar>void print(const Scalar scalar) const
  {
    cout<<scalar;
  }
  template<int lineBreakDim=0,int d=0,class Scalar>void print(const std::complex<Scalar> complexScalar) const
  {
    cout<<complexScalar.real()<<"+i*"<<complexScalar.imag();
  }
  template<int lineBreakDim=0,int d=0,class SubComponentSType>void print(const Array<SubComponentSType> array) const
  {
    if(array.pData!=NULL)
      for(Count i=0;i<array.pData->componentSCount;i++)
    {
      print<lineBreakDim,d+1>(array[i]);
      if(d==lineBreakDim)
        cout<<"\n";
      else
        cout<<" ";
    }
  }

public:

  Array() noexcept
    :pData(NULL)
  {//zur spaeteren Verwendung in den Operatoren & fuer Kompatibiltaet als ComponentSType: denn die Komponenten sind sowieso Eigentum der UeberKomponenten
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1a");
  }

  Array(const Array&  array) noexcept
  {//zur formgerechten Rückgabewerterzeugung in Operatoren: - alle noetigen Merkmale sind im paramArray, sodass abgeleitete Klassen keinen anderen OperatorDefinitionen brauchen
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1c");
    useData(array.pData);
  }
  Array(      Array&& array) noexcept
    : pData(array.pData)
  {//move constructor: will be only slightly faster.
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1c2");
    array.pData=NULL;
  }

  Array(const Count componentSCount,int *nullP) noexcept
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1b "<<componentSCount);
    useData(new Data(componentSCount));
  }

  Array(const Array<Count>& sizeArray,int *nullP) noexcept
    : pData(NULL)
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1b2\n  Array(&sizeArray) "<<sizeArray);
    resize(sizeArray,false);
    //cout<<"N_Arr(sArr,NULL)\n"<<flush;
  }
  // ambigeous wrt contructor Array(const Count componentSCount,const ComponentSType v,int *nullP)
  // Array(const Array<Count>& sizeArray,auto v,int *nullP) noexcept
  //  : pData(NULL)
  // {
  //   __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1b2\n  Array(&sizeArray) "<<sizeArray);
  //   resize(sizeArray,false);
  //   *this=v;
  //   //cout<<"N_Arr(sArr,NULL)\n"<<flush;
  // }

  Array(const Count componentSCount,const ComponentSType v,int *nullP) noexcept
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1b2 "<<componentSCount);

    useData(new Data(componentSCount));
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]=v;
  }

  template<class Value> Array(const Value v) noexcept // allows to mimic built-in scalar types
                                                      // also initialisation with Array of different depths or fundamentalType
    : pData(NULL) // needed by opr= -- assumed to be optimised away when building template instances with if(pData==NULL) statements in opr=
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1b3b "<<v);

    *this=v;
    // useData(new Data(1)); // Old definition for scalar types only: less assigments/no conditional jumps but these are potentially optimised away.
    // pData->pComponents[0]=scalar;
  }

  template<Count componentSCount> Array(const ComponentSType (&components)[componentSCount]) noexcept// allows construction from constant Arrays
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1b4 "<<componentSCount);

    useData(new Data(componentSCount));
    for(Count i=0;i<pData->componentSCount;++i)
       pData->pComponents[i]=ComponentSType(components[i]);
  }

  //These constructors would allow for direct & individual element initialisation, but there is a syntax problem of passing a variable of type initializer_list to new instead of using the {param...} form.
  //template<class Scalar> Array(std::initializer_list<Scalar> list) noexcept  //without this constructor, it is still possible to construct an Array<.> from an initializer_list: it is automatically converted into an c-array and the Array<.> construction appears subsequently.
  //{
  //  __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1d_initializer_list"<<list.size());
  //  useData(new Data(list));
  //}
  //template<class... Scalars> Array(Scalars &&... list) noexcept  //without this constructor, it is still possible to construct an Array<.> from an initializer_list: it is automatically converted into an c-array and the Array<.> construction appears subsequently.
  //{
  //  __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1d_initializer_list"<<list.size());
  //
  //  useData(new Data(std::forward<Scalars>(list)...));
  //}

  template<class Scalar> Array(std::initializer_list<Scalar> list) noexcept  //without this constructor, it is still possible to construct an Array<.> from an initializer_list: it is automatically converted into an c-array and the Array<.> construction appears subsequently.
    : pData(NULL)
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1d_initializer_list"<<list.size());

    *this=list;
  }

  template<class Invocable, class... Parameter, std::enable_if_t<std::is_invocable_r<ComponentSType,Invocable,Count,Parameter...>::value>* =nullptr >
    Array(Count componentSCount,Invocable &&invocable,Parameter... parameter)
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_1d_lambda");

    useData(new Data(componentSCount));
    for(auto i=0;i<componentSCount;++i)
      pData->pComponents[i]=invocable(i,parameter...);
  }


  ~Array() noexcept
  {//zur formgerechten Rueckgabewerterzeugung in Operatoren: - alle noetigen Merkmale sind im paramArray, sodass abgeleitete Klassen keinen anderen OperatorDefinitionen brauchen
    //cout<<"destructor; ";
    finishDataSUse();
  }

/*
  static Array create(const Count componentSCount) noexcept
  {
    //cout<<"<create1\n";
    Array array(componentSCount,NULL);
    return array;
  }
  static Array create(const Array<Count>& sizeArray) noexcept
  {
    //cout<<"<create2\n";
    Array array(sizeArray,NULL);
    return array;
  }
  static Array                  create(const Count componentSCount,const ComponentSType v) noexcept
  {  return (create(componentSCount)=v);  }
  template<class T>static Array create(const Array<Count>& sizeArray,const T v) noexcept
  {  return (create(sizeArray)=v);  }
  template<class T, Count componentSCount>static Array create(const T (&components)[componentSCount]) noexcept
  {
    //cout<<"<create5\n";
    Array array(componentSCount,NULL);
    for(Count i=0;i<componentSCount;++i)
      array[i]=ComponentSType(components[i]);
    return array;
  }
  template<class T>static Array create(const T *pComponents,const Count componentSCount) noexcept
  {
    //cout<<"<create6\n";
    Array array(componentSCount,NULL);
    for(Count i=0;i<componentSCount;++i)
      array[i]=ComponentSType(pComponents[i]);
    return array;
  }
*/

  void resize(const Count& componentSCount,bool keepComponents=false) noexcept
  {//bei Aenderung derefer() ueberpruefen
    __ARRAY_DEBUG__print__("Array::resize("<<componentSCount<<","<<int(keepComponents)<<") with useCount="<<( pData==NULL ? -1 : pData->useCount));

    if(pData==NULL)
      useData(new Data(componentSCount));
    else
    {
      //cout<<"resize,useCount="<<pData->useCount<<"; ";
      ComponentSType* pNewComponents=new ComponentSType[componentSCount];

      if(keepComponents)
      {
        //cout<<"resize,copy; ";
        for(int i=( componentSCount<pData->componentSCount ? componentSCount : pData->componentSCount )-1;i>=0;--i)
          pNewComponents[i]=pData->pComponents[i];
      }

      if(pData->useCount<=1)
        delete[] pData->pComponents;
      else
      {
        finishDataSUse();
        useData(new Data);
      }

      pData->pComponents     =pNewComponents;
      pData->componentSCount =componentSCount;
    }

    #ifdef __ARRAY_DEBUG
      if(componentSCount==0)
        std::cout<<"WARNING: requesting Array of size 0"<<std::endl<<std::flush;
      if(pComponents==NULL)
        std::cout<<"ERROR: Array storage could not be allocated."<<std::endl<<std::flush;
    #endif
  }
  template< class ScalarUser=void> inline void resize(const Array<Count> &sizeArray,bool keepComponents=false) noexcept
  {
    resize<ScalarUser>(*this,sizeArray,keepComponents);
  }


  Count size() const noexcept //returns size of this Array, for compatibility with std::array
  {  return ( pData==NULL ? 0 : pData->componentSCount );  }


  template<class ScalarUser=void,class ScalarIndex=Count>Array<Count> getSize(Array<ScalarIndex> indexArray=0) const noexcept
                                                                                               //accesses the array at index specified by the first entry in the index array, than repeats this for subarrays but increases the index of the value taken from the index array by one each time until either there is no subarray anymore or all entries of the index array were used.
  {                                                                                            //before each access the size of (sub)array is stored in an array, that will be returned eventually, at same index as the index used to access the index array.
    #ifdef __ARRAY_DEBUG__
    /**/if(lineareAlgebraSBloedeKommentare){
      cout<<"lineareAlgebraSBloedeKommentare_getSize()\n"<<flush;
      Array<Count> tmp=getSize<ScalarUser>(*this,indexArray);
      cout<<tmp<<endl<<flush;
      cout<<"lineareAlgebraSBloedeKommentare_n_getSize()\n"<<flush;
      return tmp;
    }
    #endif
    return getSize<ScalarUser>(*this,indexArray);
  }


  //an unifified short-handed interface for the resize/size methods
  template<class ScalarUser=void> std::conditional_t<is_same_v<ScalarUser,void>,Count,Array<Count>>   n()  const noexcept // Non-auto return type prevents deduction errors ('used before deduced' when n() is called/argument of decltype before it's definition)
  {                                                                                                                       // returns size of this Array, for compatibility with std::array
    if constexpr(is_same_v<ScalarUser,void>)
      return size();
    else
      return getSize<ScalarUser>();
  }
  template<class ScalarUser=void,class N> Array<ComponentSType>& n(N n,bool keepComponents=false) noexcept
  {
    if constexpr(!is_subsettable_v<N>)
      resize            (n,keepComponents);
    else
      resize<ScalarUser>(n,keepComponents);

    return *this;
  }



  void derefer() noexcept
  {
    //cout<<"derefer; ";
    __ARRAY_DEBUG__print__("Array::derefer with useCount="<<( pData==NULL ? -1 : pData->useCount ));

    if(pData->useCount>1)
      resize(pData->componentSCount,true);// legt eigenen Speicher an und kopiert Inhalt
  }

  void finishUse() noexcept
  {//dadurch gibt der Container die Referenz auf den Inhalt, wodurch die Variable effektiv frueher dekonstruiert wird, was zB. bei Parameter-Rueckgabe-Beziehungen sinnvoll sein kann
    finishDataSUse();
    pData=NULL;
  }


  // the index operator versions:
  // We do not want a rvalue version since we always can return lvalue references and passing them costs the same as passing rvalues. Is it possible to use subset a rvalued array wo an operator[]() && ?
  // Indexing proxies as "views"

  ComponentSType& operator[](const Count i) & noexcept
  {
    __ARRAY_DEBUG__printV__("Array::[](Count="<<i<<")",2);
    __ARRAY_DEBUG__assert__(i<0||i>=n(),"ERROR: Array::[]("<<i<<") but size=="<<size());

    derefer();//if a non-const element is needed, a modification must be asumed. Transports copy-data-call from base array with multiple uses to all sub...subarrays along the index chain
    return(pData->pComponents[i]); // Modifications must be also made at basicAt().
  }
  ComponentSType&& operator[](const Count i) && noexcept
  {
    __ARRAY_DEBUG__printV__("Array::[](Count="<<i<<") &&",2);
    __ARRAY_DEBUG__assert__(i<0||i>=n(),"ERROR: Array::[]("<<i<<") but size=="<<size());

    derefer();//if a non-const element is needed, a modification must be asumed. Transports copy-data-call from base array with multiple uses to all sub...subarrays along the index chain
    return(static_cast<ComponentSType&&>(pData->pComponents[i])); // Modifications must be also made at basicAt().
  }

  const ComponentSType& operator[](const Count i) const & noexcept
  {
    __ARRAY_DEBUG__printV__("Array::[](Count="<<i<<") const",2);
    __ARRAY_DEBUG__assert__(i<0||i>=n(),"ERROR: Array::[]("<<i<<") but size=="<<size());

    return(pData->pComponents[i]); // Modifications must be also made at basicAt().
  }
  const ComponentSType&& operator[](const Count i) const && noexcept
  {
    __ARRAY_DEBUG__printV__("Array::[](Count="<<i<<") const &&",2);
    __ARRAY_DEBUG__assert__(i<0||i>=n(),"ERROR: Array::[]("<<i<<") but size=="<<size());

    return(static_cast<const ComponentSType&&>(pData->pComponents[i])); // Modifications must be also made at basicAt().
  }
  // The at() methods for Array access via an index array...
  //
  // The implementation uses recursion to subset the multiple nested Arrays.
  // The methods for induction and anchoring are choosen via template SFINAE (OLD) / constexpr if (NEW): the induction does that by requiring a subsettable type (the operator[] must be applicable) and furthermore last template parameter tests, whether
  //   the desired return type (Scalar) has been reached. This is done via is_same< . ,Scalar> in std::enable_if which does not posess ::type if its first template argument is false (Old).
  // All template types are deduced recursively such that ambigouities are exclude. CV and constness qualifiers are supported: the constness of Scalar must not match the constness of the Array.
  // at() may also access non Array<ComponentSType>-types such as sub-Arrays or std::array since only subsettability is required.
  //  => Since access via an index-array can be provided for any subsettable type (as at(...) does), there might be also stl/boost pendants.


  template<class Scalar,bool noDerefer=false,class Index, Count i=0, class SubSettable> //template<class Scalar,class IndexArray, int i=0, class SubSettable, typename std::enable_if< !std::is_same< SubSettable, Scalar >::value,int >::type=0>
    static inline auto&& at(SubSettable &&array,  Index &&index) noexcept
  {
    //return at<Scalar,IndexArray,i+1>( array[const_cast<const IndexArray>(indexArray)[i]], std::forward<IndexArray>(indexArray) ); //const casting removes the need to define the parameter const which in return allows for template<class T>f(T &&t) universal reference.
    if constexpr( std::is_same_v< remove_cr_t<SubSettable>, remove_cr_t<Scalar> > || is_same_v<void,std::remove_pointer_t<Index>> )
      return std::forward<SubSettable>(array);
    else
      if constexpr( !is_subsettable<Index>::value && std::is_same_v< remove_cr_t<decltype(array[index])>, remove_cr_t<Scalar> > )
         return basicAt<noDerefer>(std::forward<SubSettable>(array),std::forward<Index>(index));
      else
        return at<Scalar,noDerefer,Index,i+1>( basicAt<noDerefer>(std::forward<SubSettable>(array),static_cast<const Index>(index)[i]), std::forward<Index>(index) ); //const casting removes the need to define the parameter const which in return allows for template<class T>f(T &&t) universal reference.
  }

  // Note 1: Public function due to Array<T>::at(array,Index) works for Scalar-typed array, too (in contrast to array.at(Index).) useful for template routines
  //         Example: Key(<=>Index) concept: methods  on arrays are supposed to work on the first level (otherwise method(array(Index)) is efficient.) However, applying the method for certain elements, that a characterised by a second Index, called key, is often necessary. Eg., summing over all elements of column of the matrix.
  //           Then: template<class ScalarUser=void, class T, class Key=int, Scalar=firstNonVoidType( ScalarUser, decltype(Array<Abritrary>::at(T(),Array<int>()) ) > Scalar Func(T a,Key k=0) using internally "at <Scalar> (a[i],key)" instead of pure a[i] provides the generalisation wo runtime costs (thanks to pure template wrapper solutions & dead code elemenination)
  // Note 2: If Index is void* at becomes the identity and should be optimized away. This behaviour allows to implement the trivial case in template algorithm without impacting runtime speed.


  //  An old implementation of at(): it does not make use of std::enable_if and std::is_same but of Array<T> beeing more specialised than Scalar if the Array is still subsettable. However, if Scalar is an Array,
  //                                 the non-Anchor function stays more specialised, nevertheless and recursion still stops when all Arrays are subseted.
  //
  //  template<class Scalar,class ScalarIndex,int i=0>          static Scalar at(      Scalar scalar  ,const Array<ScalarIndex>& indexArray, Scalar &&dummy) //Scalar is expected to be a reference and const if desired - this is not enforced here
  //  {  return scalar; }
  //  template<class Scalar,class ScalarIndex,int i=0,class T>  static Scalar at(      Array<T>& array,const Array<ScalarIndex>& indexArray, Scalar &&dummy)
  //  {  return  at<Scalar,ScalarIndex,i+//>( array[indexArray.pData->pComponents[i]], indexArray, dummy ); }
  //  template<class Scalar,class ScalarIndex,int i=0,class T>  static Scalar at(const Array<T>& array,const Array<ScalarIndex>& indexArray, Scalar &&dummy)

  template<class ScalarUser,class Index> using atScalar_t = first_nonVoid_t<ScalarUser,std::conditional_t<is_subsettable_v<Index>,decltype(getFundamentalType()),ComponentSType>>;

  template<class ScalarUser=void,bool noDerefer=false,class Index>   auto&& at    (Index &&index )       &  noexcept  {  return at<atScalar_t<ScalarUser,Index>,noDerefer>(          *this ,std::forward<Index>(index));  }
  template<class ScalarUser=void,bool noDerefer=false,class Index>   auto&& at    (Index &&index ) const &  noexcept  {  return at<atScalar_t<ScalarUser,Index>,noDerefer>(          *this ,std::forward<Index>(index));  }
  template<class ScalarUser=void,bool noDerefer=false,class Index>   auto&& at    (Index &&index )       && noexcept  {  return at<atScalar_t<ScalarUser,Index>,noDerefer>(std::move(*this),std::forward<Index>(index));  }
  template<class ScalarUser=void,bool noDerefer=false,class Index>   auto&& at    (Index &&index ) const && noexcept  {  return at<atScalar_t<ScalarUser,Index>,noDerefer>(std::move(*this),std::forward<Index>(index));  }
                                                                   // useful beside the operator[](Array<Scalar>) version bc the Scalar template parameter can be used to set anchor's type where the indexing stops such that not all Array<> have to be indexed
                                                                   // in [] recursion ends with Scalar!=Array<...>  bco otherwise it would be choosen as the most specialised type
                                                                   // call syntax: Array<T> a; ... ; a.template at<ST1,ST2>(...); -- the  template keyword helps g++ to evaluate the template list
  template<class ScalarUser=void,                     class Index>   auto&& atFast(Index &&index )       &  noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(          *this ,std::forward<Index>(index));  }// the const versions does not differ from ordinary at but exists in order to have a persistent interface
  template<class ScalarUser=void,                     class Index>   auto&& atFast(Index &&index ) const &  noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(          *this ,std::forward<Index>(index));  } //skips checks for multiple references to same data
  template<class ScalarUser=void,                     class Index>   auto&& atFast(Index &&index )       && noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(std::move(*this),std::forward<Index>(index));  } //skips checks for multiple references to same data
  template<class ScalarUser=void,                     class Index>   auto&& atFast(Index &&index ) const && noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(std::move(*this),std::forward<Index>(index));  }// the const versions does not differ from ordinary at but exists in order to have a persistent interface
  template<class ScalarUser=void,                     class Index>   auto&& atF   (Index &&index )       &  noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(          *this ,std::forward<Index>(index));  }
  template<class ScalarUser=void,                     class Index>   auto&& atF   (Index &&index ) const &  noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(          *this ,std::forward<Index>(index));  }
  template<class ScalarUser=void,                     class Index>   auto&& atF   (Index &&index )       && noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(std::move(*this),std::forward<Index>(index));  } //skips checks for multiple references to same data
  template<class ScalarUser=void,                     class Index>   auto&& atF   (Index &&index ) const && noexcept  {  return at<atScalar_t<ScalarUser,Index>,true     >(std::move(*this),std::forward<Index>(index));  }// the const versions does not differ from ordinary at but exists in order to have a persistent interface

  // if we want to ensure IndexArray to be subsettable is_array<T> might be used:
  // However, there is no danger to confuse with operator[](Scalar) since the call of at() for the fundamental type (which is ComponentSType or even more resolved) needs at least one access to the indexArray.
  //   template<class IndexArray, typename std::enable_if<is_array<IndexArray>::value,int>::type=0> auto&& operator[](IndexArray &&indexArray)
  // In order to prevent expensive return-type deduction (which occurs every time when operator[] would be used!), we will favor that over the simpler
  //   auto&& operator[](auto &&indexArray)
  // Also: use of operator[] within the primary at<.>(.)-template requires the return type to be 'deduced'=>non auto OR defined before the at-template.
  template<class IndexArray, std::enable_if_t<is_subsettable<IndexArray>::value,int>* =nullptr>         ComponentSType&  operator[](IndexArray &&indexArray)       &  noexcept  {  __ARRAY_DEBUG__printV__("Array::[](Subsetable)       & ",2);    return at<ComponentSType>(          *this ,indexArray);  }
  template<class IndexArray, std::enable_if_t<is_subsettable<IndexArray>::value,int>* =nullptr>         ComponentSType&& operator[](IndexArray &&indexArray)       && noexcept  {  __ARRAY_DEBUG__printV__("Array::[](Subsetable)       &&",2);    return at<ComponentSType>(std::move(*this),indexArray);  }
  template<class IndexArray, std::enable_if_t<is_subsettable<IndexArray>::value,int>* =nullptr>   const ComponentSType&  operator[](IndexArray &&indexArray) const &  noexcept  {  __ARRAY_DEBUG__printV__("Array::[](Subsetable) const & ",2);    return at<ComponentSType>(          *this ,indexArray);  }
  template<class IndexArray, std::enable_if_t<is_subsettable<IndexArray>::value,int>* =nullptr>   const ComponentSType&& operator[](IndexArray &&indexArray) const && noexcept  {  __ARRAY_DEBUG__printV__("Array::[](Subsetable) const &&",2);    return at<ComponentSType>(std::move(*this),indexArray);  }


  template<bool periodic=true> Array<ComponentSType>& shift(const Count offset) noexcept
  {
    const Count n=size();
    if constexpr (periodic)
    {
      derefer();
      for(Count i=0; true;) // swap all elements of distance offset in a cyle. Cycle starting at i.
      {
        for(Count j=(i-offset+n)%n; j!=i ; j=(j-offset+n)%n) // swap elements until the cycle is completed <=> element i reached again.
          std::swap(pData->pComponents[i],pData->pComponents[j]);
        ++i;
        if( n%offset==0 ? i==offset : offset%(n%offset)!=0||i==n%offset ) // o%(n%o)==0 <=> cyclic swaping of elements of distance o=offset happens in a cycle that does not cover all elements
          break;
      }
    }
    else
    {
      __ARRAY_DEBUG__assert__(n<-offset,"ERROR: Array.shift<false>("<<offset<<") but size=="<<size());
      Data* pD=new Data(n+offset);
        if(offset<0)
          for(Count i=pD->componentSCount;i-->0;)
            pD->pComponents[i]        = pData->pComponents[i-offset];
        else
          for(Count i=n;i-->0;)
            pD->pComponents[i+offset] = pData->pComponents[i];

      finishDataSUse();
      useData(pD);
    }

    return *this;
  }


  template<class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type>Array& operator*=(const Scalar scalar)
  {
    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]*=scalar;
    return *this;
  }

  template<class ParameterComponentSType>Array& operator*=(const Array<ParameterComponentSType> &array)
  {
    __ARRAY_DEBUG__assert__(array.size()!=size(),"ERROR: Array::*=(&array): size=="<<size()<<"!=array.size=="<<array.size());

    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]*=array.pData->pComponents[i];
    return *this;
  }

  template<class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type>Array& operator/=(const Scalar scalar)
  {
    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]/=scalar;
    return *this;
  }

  template<class ParameterComponentSType>Array& operator/=(const Array<ParameterComponentSType> &array)
  {
    __ARRAY_DEBUG__assert__(array.size()!=size(),"ERROR: Array::/=(&array): size=="<<size()<<"!=array.size=="<<array.size());
    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]/=array.pData->pComponents[i];
    return *this;
  }

  template<class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type>Array& operator+=(const Scalar scalar)
  {
    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]+=scalar;
    return *this;
  }

  template<class ParameterComponentSType>Array& operator+=(const Array<ParameterComponentSType> &array)
  {
    __ARRAY_DEBUG__assert__(array.size()!=size(),"ERROR: Array::+=(&array): size=="<<size()<<"!=array.size=="<<array.size());
    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]+=array.pData->pComponents[i];
    return *this;
  }

  template<class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type>Array& operator-=(const Scalar scalar)
  {
    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]-=scalar;
    return *this;
  }

  template<class ParameterComponentSType>Array& operator-=(const Array<ParameterComponentSType>& array)
  {
    __ARRAY_DEBUG__assert__(array.size()!=size(),"ERROR: Array::-=(&array): size=="<<size()<<"!=array.size=="<<array.size());
    derefer();
    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]-=array.pData->pComponents[i];
    return *this;
  }


  template<class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type>inline Array& operator=(const Scalar scalar) noexcept
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_2");

    if(pData==NULL)//provides compatibility with built-in types in the case of  Array<T> x; x=t; with t in T; T could be double etc
      useData(new Data(1));
    else
      derefer();

    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]=scalar;
    return *this;
  }

  template<class ArrayT, class=typename enable_if< inheritsArray<ArrayT>::value&&!std::is_same<Array,ArrayT>::value >::type>inline Array& operator=(ArrayT &&array) noexcept
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_2b");

    if(size()!=array.size())//provides compatibility with built-in types in the case of  Array<T> x; x=t; with t in T; T could be double etc
      resize(array.size());
    else
      derefer();

    for(Count i=0;i<pData->componentSCount;++i)
      pData->pComponents[i]=std::as_const(array)[i];
    return *this;
  }

  inline Array& operator=(const Array &array) noexcept
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_3");

    if(this!=&array)
    {
      finishDataSUse();
      useData(array.pData);
    }

    return *this;
  }
  inline Array& operator=(Array &&array) noexcept //only marginaly faster, though the operator allows to assign temporaries as from =(a+b); However, thanks to the container concepts operators
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_3b_opr=(Array&&)");

    if(this!=&array)
    {
      finishDataSUse();
      pData=array.pData;
      array.pData=NULL;
    }
    return *this;
  }

  template<class Scalar> inline Array& operator=(const std::initializer_list<Scalar> list) noexcept
  {
    __ARRAY_DEBUG__print__("lineareAlgebraSBloedeKommentare_3c_opr=(init-list); list.size="<<list.size());

    if(size()!=list.size()) // use the default-interface for internal needs, too! (<-is optimised decently)
     resize(list.size());
    auto l=list.begin();
    for(auto i=0; i<pData->componentSCount; l++, ++i)
      pData->pComponents[i]=*l; // only implicit conversions for safety.
    return *this;
  }




  // for implicit type conversions to ComponentSType:
  // allows to mimic ComponentSType (eg. int) such that generalistions of scalar algoritms is easy.

  operator ComponentSType&() noexcept
  {
    derefer();
    return (pData->pComponents[0]);
  }

  operator ComponentSType() const noexcept
  {
    return (pData->pComponents[0]);
  }

  //


  template<int lineBreakDim=0> void print()
  {
    print<lineBreakDim>(*this);
  }
};



/*OPERATOREN:
-parameter darf kein temporary(const &) sein, denn dann wuerde pData->useCount auch 1 sein, wenn eine lokale variable durchgereicht wird. Andererseits besaesse auch ein nur als temporary-Variabel existierender RueckgabeWert die NutzungsAnzahl. Mit diesem System kann man also die Parameter nicht recyclen
-der zuerst angegebene Parameter muss bzgl. der Array-template-Schachtelung (Array< Array<...> >) der komplexere sein, denn aus ihm leitet sich der RueckgabeTyp ab
--mann koente diesbezueglich 2 templateVersionen ueberladen, wo die Rolle des 1. und 2. Parameters vertauscht. Der Compiler kann die passenden Version finden, indem man eine Operation Array-Membern des ComponentSType ausfuehrt. Am RekursionsEnde passt dann hier nur der komplexere Typ. Allerdings ist das fuer gleichKomplexe Parameter uneindeutig!
--SkalarOperationen sollten deshalb primaer als 2. Parameter eingeplant werden
-der zuerst angegebene Parameter muss einen kleineren oder gleichen Wert fuer pData->componentSCount besitzen, arbeiten mit dieser Information
*/
template<class ComponentSTypeA,class ComponentSTypeB> ComponentSTypeA operator*(const Array<ComponentSTypeA> arrayA,const Array<ComponentSTypeB> arrayB)
{
  std::cout<<"depricate Array*Array";
  ComponentSTypeA componentSTypeResult(arrayA.pData->pComponents[0]*arrayB.pData->pComponents[0]);

  for(Count i=1;i<arrayA.pData->componentSCount;++i)
    componentSTypeResult+=arrayA.pData->pComponents[i]*arrayB.pData->pComponents[i];

  return componentSTypeResult;
}
template<class FieldA,class FieldB> FieldA mult(FieldA a,const FieldB b) // elementwise multiplication. Because, here, we need to implement the scalar version too, we do not restrict the parameter fields via Array<.> (useful for templates with support for Scalar types). Nevertheless, the definition is designed around the requirements of Array<.>
{
  FieldA result(a*=b);//derefer() geschieht bei .= Operationen automatisch
  if constexpr( inheritsArray_v<FieldA> )
    a.finishUse();
  return result;
}


template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator*(Array<ComponentSType> array,const Scalar scalar)
{
  Array<ComponentSType> resultArray(array*=scalar);//derefer() geschieht bei .= Operationen automatisch
  array.finishUse();
  return resultArray;
}
template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator*(const Scalar scalar,Array<ComponentSType> array)
{
  Array<ComponentSType> resultArray(array*=scalar);//derefer() geschieht bei .= Operationen automatisch
  array.finishUse();
  return resultArray;
}


template<class ComponentSTypeA,class ComponentSTypeB> Array<ComponentSTypeA> operator/(Array<ComponentSTypeA> arrayA,const Array<ComponentSTypeB> arrayB)
{
  Array<ComponentSTypeA> resultArray(arrayA/=arrayB);
  arrayA.finishUse();
  return resultArray;
}
template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator/(Array<ComponentSType> array,const Scalar scalar)
{
  Array<ComponentSType> resultArray(array/=scalar);//derefer() geschieht bei .= Operationen automatisch
  array.finishUse();
  return resultArray;
}
template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator/(const Scalar scalar,Array<ComponentSType> array)
{
  array.derefer();
  for(Count i=0;i<array.pData->componentSCount;i++)
    array.pData->pComponents[i] = scalar/array.pData->pComponents[i];
  return array;
}


template<class ComponentSTypeA,class ComponentSTypeB> Array<ComponentSTypeA> operator+(Array<ComponentSTypeA> arrayA,const Array<ComponentSTypeB> arrayB)
{
  //Array<ComponentSTypeA> resultArray;//lokale Variablen, die zurueckgegeben werden, werden nicht kopiert und die lokale Instanz destruiert, sondern durchgereicht;, falls dies nicht standard festgelegt ist, dann kann man noch zurueckgegebenen Arrays die useCount==0 mitgeben und per statusVar mitteilen, dass pData nur dann geloescht werden soll, wenn der rueckgabeWert nicht verwendet wurde - dh. useCount==0

  //cout<<"opr+,ASuseCount="<<arrayA.pData->useCount<<",BSuseCount="<<arrayB.pData->useCount<<"; ";
  /*geht nicht weil, der rueckgabeTyp bei umgedrehter Rechnung geandert wird*///if(arrayA.pData->useCount<=1)
  //{
    Array<ComponentSTypeA>  resultArray(arrayA+=arrayB);//wenn A nirgenwo anders verwendet wird, dann koennen wir es hier verwenden, ohne neuen anzulegen.
    arrayA.finishUse();
  //}
  //else
  //{
  //  resultArray=(arrayB+=arrayA);//wenn A nirgenwo anders verwendet wird, dann koennen wir es hier verwenden, ohne neuen anzulegen.
  //  arrayB.finishUse();
  //}

  return resultArray;//derefer() geschieht bei .= Operationen automatisch
}
template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator+(Array<ComponentSType> array,const Scalar scalar)
{
  Array<ComponentSType> resultArray(array+=scalar);//derefer() geschieht bei .= Operationen automatisch
  array.finishUse();
  return resultArray;
}
template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator+(const Scalar scalar,Array<ComponentSType> array)
{
  Array<ComponentSType> resultArray(array+=scalar);//derefer() geschieht bei .= Operationen automatisch
  array.finishUse();
  return resultArray;
}


template<class ComponentSTypeA,class ComponentSTypeB> Array<ComponentSTypeA> operator-(Array<ComponentSTypeA> arrayA,const Array<ComponentSTypeB> arrayB)
{
  //Array<ComponentSTypeA> resultArray;//lokale Variablen, die zurueckgegeben werden, werden nicht kopiert und die lokale Instanz destruiert, sondern durchgereicht;, falls dies nicht standard festgelegt ist, dann kann man noch zurueckgegebenen Arrays die useCount==0 mitgeben und per statusVar mitteilen, dass pData nur dann geloescht werden soll, wenn der rueckgabeWert nicht verwendet wurde - dh. useCount==0

  //cout<<"opr-,ASuseCount="<<arrayA.pData->useCount<<",BSuseCount="<<arrayB.pData->useCount<<"; ";
  /*geht nicht weil, der rueckgabeTyp bei umgedrehter Rechnung geandert wird*///if(arrayA.pData->useCount<=1)||arrayB.pData->useCount>1)
  //{
    Array<ComponentSTypeA> resultArray(arrayA-=arrayB);//wenn A nirgenwo anders verwendet wird, dann koennen wir es hier verwenden, ohne neuen anzulegen.
    arrayA.finishUse();
  //}
  //else
  //{
  //  //printf("HeuPferd; ");
  //  resultArray=((arrayB-=arrayA)*=(-1));//wenn A nirgenwo anders verwendet wird, dann koennen wir es hier verwenden, ohne neuen anzulegen, allerdings muss das vorzeichen aller Elemente nochmal umgedreht werden
  //  arrayB.finishUse();
  //}

  return resultArray;
}
template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator-(Array<ComponentSType> array,const Scalar scalar)
{
  Array<ComponentSType> resultArray(array-=scalar);//derefer() geschieht bei .= Operationen automatisch
  array.finishUse();
  return resultArray;
}
template<class ComponentSType,class Scalar, class=typename enable_if_notInheritsArray<Scalar>::type> Array<ComponentSType> operator-(const Scalar scalar,Array<ComponentSType> array)
{
  Array<ComponentSType> resultArray((array-=scalar)*=(-1));//derefer() geschieht bei .= Operationen automatisch
  array.finishUse();
  return resultArray;
}


template<class ComponentSTypeA,class ComponentSTypeB> bool operator==(const Array<ComponentSTypeA>& arrayA,const Array<ComponentSTypeB>& arrayB) noexcept
{
  if(arrayA.pData==arrayB.pData)
    return(true);
  if(arrayA.size()!=arrayB.size())
    return(false);
  for(Count i=0;i<arrayA.pData->componentSCount;i++)
    if(!(arrayA.pData->pComponents[i]==arrayB.pData->pComponents[i]))
      return(false);

  return(true);
}

template<class ComponentSTypeA,class ComponentSTypeB> bool operator!=(const Array<ComponentSTypeA>& arrayA,const Array<ComponentSTypeB>& arrayB) noexcept
{
  return !(arrayA==arrayB);
}

int arraySStreamOprOutSWidth=0;
int arraySStreamOprOutSLineBreakDimension=2;//to get this working calls on components must be overtaken by a function with same signature but parameter 'dimension' and the linebreak action eg at the beginning of the definition + passing dimension-1 for component calls
template<class ComponentSType> std::ostream& operator<<(std::ostream& ostream,const Array<ComponentSType>& array)
{
  if(array.pData!=NULL)
  {
    for(Count i=0;i<array.pData->componentSCount;++i)
    {
      ostream << (array.pData->pComponents[i]) << " ";
      if(arraySStreamOprOutSWidth>0 && ostream.tellp()<int(i)*(arraySStreamOprOutSWidth+1))
      {
        Count fillN=Count(int(i)*(arraySStreamOprOutSWidth+1)-ostream.tellp());
        for(Count i=0;i<fillN;++i)
          ostream<<" ";
      }
    }
  }
  return ostream;
}
// a specialisation for std::complex is unwanted for file IO...
//
//template<class Scalar> std::ostream& operator<<(std::ostream& ostream,const Array< std::complex<Scalar> >& array)
//{//could be avoide together with Array::print(*) specialisation for std::complex if an << is defined for complex directly
//  if(array.pData!=NULL)
//    for(Count i=0;i<arrayA.pData->componentSCount;++i)
//      ostream << (array.pData->pComponents[i].real()) << "+i*" << (array.pData->pComponents[i].imag()) << " ";

//  return ostream;
//}
