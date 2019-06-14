// Traits & helper structs used by array.h


// trait-like to select 1fst non-void type.
//   Used to overide template dependent type B with user-defined tpl argument A, that occurs for convenience before B in the list of declarations of template arguments
template<class A, class B> struct first_nonVoid           {  typedef A type;  };
template<         class B> struct first_nonVoid <void,B>  {  typedef B type;  };

template<class A,class B> using first_nonVoid_t = typename first_nonVoid<A,B>::type;
template<class A,class B> using nonVoid_t       = typename first_nonVoid<A,B>::type;

#define keyClasses(Key,Scalar,ScalarUser,ComponentSType) class Key=void*, class Scalar=nonVoid_t<std::remove_reference_t<ScalarUser>,ComponentSType> // if the key shold be used by setting Scalar!=ComponentSType and Key=int has not be overridden by providing a Subsetable-typed Argument, Array<>::at<...> will fail via key[i]


// constexpr to forward an lvalue "member", considered to behave as a member of B base, as an rvalue if base is one and as lvalue otherwise...
template <class B,class M,enable_if_t< is_reference_v<B> && is_reference_v<M> >* =nullptr> constexpr static std::remove_reference_t<M>&  forwardMember(M &&member) noexcept //M &&member is universal reference but should be a lvalue in most cases
{  return member;  }                                //due to std::remove_reference<B> && is not universal anymore, ie. its possible to distinguish beween & and &&
template <class B,class M,enable_if_t<!is_reference_v<B> >* =nullptr>                      constexpr static std::remove_reference_t<M>&& forwardMember(M &&member) noexcept //M &&member is universal reference but should be a lvalue in most cases
{  return static_cast<std::remove_reference_t<M>&&>(member);  }



// test for operator[](Count) structue for SFINAE - could be discarded. Reduces time usage for compilation in non-static at(IndexAray) methods.
template<class  , class = void> struct is_subsettable                                                                     : std::false_type  { };
template<class T>               struct is_subsettable<T, std::void_t<decltype(std::declval<T>()[std::declval<Count>()])>> : std::true_type   { };

template<class T> inline static constexpr bool is_subsettable_v = is_subsettable<T>::value;
  //std::decval<T>() could be implemented by *static_cast<T*>(NULL) - this is meaned for use with decltype(.) and similar needs of a value signature. Accessing the value explicitly results in the error or dereferencing a null pointer


template<class T> struct inheritsArray //tests inheritance neglecting ComponentSType
{
  template<template<class> class ArrayType,class ComponentSTypeL> static std::is_base_of<Array<bool>,ArrayType<bool>>  check( ArrayType<ComponentSTypeL> ); //a more specialised template to extract the outer template ArrayType and test for inheritance at equal template Parameter ComponentSType=bool
  template<class TT>                                              static std::false_type                               check( TT ); //const is not necessairy, but it adds an additional reason not to select this definition -- in case the template specialisation measure changes.

  static const bool value = decltype(check(std::declval<typename std::decay<T>::type>()))::value;
};
template<class T> inline static constexpr bool inheritsArray_v = inheritsArray<T>::value;
template<class T> struct enable_if_inheritsArray    : enable_if< inheritsArray_v<T>> {};
template<class T> struct enable_if_notInheritsArray : enable_if<!inheritsArray_v<T>> {};      // if Array is inherited from  a class "D", than function with unconstrained template arguments (often named Scalar) might be prefered over functions with Array<.> arguments:
                                                                                              // the upcast D into Array is disfavoured over using the template to create an method directly for D: template<class Scalar> f(Scalar) with [Scalar=D]
                                                                                              // Implicit conversion Array<ComponentSType> into ComponentSType allow every method with Scalar template arguments to be used that without errors.
template<class T> using enable_if_inheritsArray_t    = typename enable_if_inheritsArray   <T>::type;
template<class T> using enable_if_notInheritsArray_t = typename enable_if_notInheritsArray<T>::type;


template<class T> struct fundamental // this does not takes & && (inherited from T ) into account. In most cases, the projection of the reference-attributes from T onto the fundamentalType is not desired, anyway.
                                     // though, it would be easy to add &/&& support bco Array::opr[] is propagating the reference type.
{
  static auto get()  { if constexpr (inheritsArray_v<T>) return typename fundamental<decltype(std::declval<T>()[0])>::type{}; else return std::remove_reference_t<T>{}; } //T{} <=> T&& which is, what std::declval does. (uuuuuuHowever, declval appears do not by applicable on certain type, as double& -> Scalar&-types ?)
  typedef decltype(get()) type;
};
template<class T> using fundamental_t = typename fundamental<T>::type;


template< class SubSettable, class Scalar=fundamental_t<SubSettable> > struct DimensionCount
{
  template<class SubSettableL=SubSettable,Count d=0> static constexpr Count get()
  { if constexpr ( !is_subsettable_v<SubSettableL> || std::is_same<decay_t<Scalar>,decay_t<SubSettableL>>::value )   return d;  else  return get<decltype(std::declval<SubSettableL>()[0]),d+1>(); }
  static const Count value = get();
};
template< class SubSettable, class Scalar=fundamental_t<SubSettable> > inline static constexpr Count rank = DimensionCount<SubSettable,Scalar>::value;

template<Count rank,template<class... > class Container,class... FundamentalT> struct nested
{
  template<class T> struct id // std::declval<T,,,>() didn't work but would be preferable. Also T...{} didn't work. -> id<T...>::type{}. T... is expected to have only one element
  { typedef T type; };
  template<Count d, class... T> static constexpr auto get()
  { if constexpr (d==0)  return typename id<T...>::type{};  else  if constexpr (d==1) return Container<T...>{}; else return get<d-1,Container<T...>>(); }
  typedef decltype(get<rank,FundamentalT...>()) type;
};
template<Count rank,template<class ...> class Container, class... FundamentalT> using nested_t = typename nested< rank, Container, FundamentalT... >::type;
template<Count rank,                                     class    FundamentalT> using tensor_t = typename nested< rank, Array    , FundamentalT    >::type;
