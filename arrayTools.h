#pragma once
#include<tuple>


namespace ArrayTools
{
  template<class SubSettable,class Scalar=fundamental_t<SubSettable>> inline static constexpr Count rank=Array<bool>::rank<SubSettable,Scalar>;


  template<class Tt,class Ts> struct copy_referenceType //similar to forward, can also handle non-reference types
  {
    template<class=typename std::enable_if< std::is_reference<Ts>::value>::type > static typename remove_reference<Tt>::type &&  helper(typename remove_reference<Ts>::type &&);
    template<class=typename std::enable_if< std::is_reference<Ts>::value>::type > static typename remove_reference<Tt>::type  &  helper(typename remove_reference<Ts>::type & );
    template<class=typename std::enable_if<!std::is_reference<Ts>::value>::type > static typename remove_reference<Tt>::type     helper(                          Ts          );
    //typename remove_reference<Tt>::type     helper(typename remove_reference<Ts>::type   );

    typedef decltype(helper(std::declval<Ts>())) type;
  };


  template<class  ,class=void> struct has_atFast                                                             : std::false_type {  };
  template<class T>            struct has_atFast<T,std::void_t<decltype(std::declval<T>().atFast(Count{}))>> : std::true_type  {  };
  template<class T> inline static constexpr bool has_atFast_v = has_atFast<T>::value;
  template<class  ,class=void> struct has_at                                                                 : std::false_type {  };
  template<class T>            struct has_at    <T,std::void_t<decltype(std::declval<T>().at    (Count{}))>> : std::true_type  {  };
  template<class T> inline static constexpr bool has_at_v     = has_atFast<T>::value;

  //template<class  ,class... P, class=void> struct has_n                                                                          : std::false_type {  };
  //template<class T,class... P>             struct has_n<T,P...,std::void_t<decltype(std::declval<T>().n(std::declval<P>()...))>> : std::true_type  {  };

  //template<class T,class... P> inline static constexpr bool has_n_v = has_n<T,P...>::value;


  template<class ScalarUser=void,class SubSettable,class Index>                                                        inline auto&& at    (SubSettable&& s,Index&& i)
  {
    if constexpr(inheritsArray_v<SubSettable>)
      return std::forward<SubSettable>(s).template at<ScalarUser>(std::forward<Index>(i));
    else
      return std::forward<SubSettable>(s)[std::forward<Index>(i)];
  }
  template<class ScalarUser=void,class SubSettable,class Index>                                                        inline auto&& atFast(SubSettable&& s,Index&& i)
  {
    if constexpr(has_atFast_v<SubSettable>)
      if constexpr(inheritsArray_v<SubSettable>)
        return std::forward<SubSettable>(s).template atFast<ScalarUser>(std::forward<Index>(i));
      else
        return std::forward<SubSettable>(s).         atFast            (std::forward<Index>(i));
    else
      return std::forward<SubSettable>(s)[std::forward<Index>(i)];
  }


  using std::size;
  template<class Array,class Size> inline auto&& resize(Array &&a,Size &&size,bool keepComponents=false)
  {
    if constexpr(inheritsArray_v<Array>)
      a.resize(std::forward<Size>(size),keepComponents);
    else
      a.resize(std::forward<Size>(size));
    return std::forward<Array>(a);
  }

  template<class... P,class Array>            auto&& n(Array&& a)
  {
    //if constexpr(has_n_v<Array>)
      return std::forward<Array>(a).template n<P...>();
    //else
    //  return Count(1);
  }
  template<class... P,class Array,class Size> inline auto&& n(Array&& a,Size&& size, bool keepComponents=false)
  {
    return std::forward<Array>(a).template n<P...>(std::forward<Size>(size),keepComponents);
  }



  namespace // internals
  {

    template<class ComponentSType> struct SubArraySArrayLoop
    {
      Array < ComponentSType > result;
      const Count& dimensionIdRemoveSIndex;

      template<class ScalarType>void process(const Array< ScalarType > & array __attribute__ ((unused)),double result __attribute__ ((unused)),const Count dimensionIdRemaining __attribute__ ((unused)))
      {  }//template recursion anchor

      template<class SubComponentSType>void process(const Array< Array< SubComponentSType > >& array,Array< SubComponentSType >& result,const Count dimensionIdRemaining)
      {
        result.resize(size(array));
        if(dimensionIdRemaining==1)
        {
          for(Count i=0;i<size(array);++i)
            atFast(result,i)=array[i][dimensionIdRemoveSIndex];
        }
        else
          for(Count i=0;i<size(array);++i)
            process(array[i],atFast(result,i),dimensionIdRemaining-1);
      }


      SubArraySArrayLoop(const Array< Array<ComponentSType> >& array,const int& dimensionIdRemove,const Count& i)
        :dimensionIdRemoveSIndex(i)
      {
        process (array,result,dimensionIdRemove);
      }
    };

  }



  template<class Invocable> struct get_args
  {
    template<class R,class... Args>               static std::tuple<Args...>                    helper( R(*)(Args...) );                   // ptr to func
    template<class R,class Closure,class... Args> static std::tuple<Args...>                    helper( R(Closure::*)(Args...) );          // ptr to member func
    template<class R,class Closure,class... Args> static std::tuple<Args...>                    helper( R(Closure::*)(Args...) const);     // ptr to member func const
    template<class Closure>                       static decltype(helper(&Closure::operator())) helper( Closure );                         // lambda
    //template<class Closure,class=typename enable_if<is_class<Closure>::value>::type> static auto                helperA( Closure ) -> decltype(helperA(&Closure::operator()));  // lambda

    typedef decltype(helper(std::declval<Invocable>())) type;
  };
  template<class Invocable> using get_args_t = typename get_args<Invocable>::type;


  template< class Array, class Invocable, class... Parameter > struct Map
  {
  public:
    typedef remove_cr_t < typename std::tuple_element<0,typename get_args<Invocable>::type>::type> Scalar;
    typedef remove_cr_t < typename std::tuple_element<1,typename get_args<Invocable>::type>::type> Index;

  private:
    Invocable                    &&invocable;
    std::tuple<Parameter &&...>    parameter;
    static const bool   subsettableIndex = is_subsettable<Index>::value;

    //template<std::enable_if_t<!subsettableIndex>* =nullptr,class ArrayL=Array> struct IndexContainer             {  };
    template<class ArrayL=Array, bool=subsettableIndex> struct IndexContainer {};
    template<class ArrayL>                              struct IndexContainer<ArrayL,true>
    {
      Index v;
      IndexContainer()
      {
        if constexpr ( inheritsArray<Index>::value )
          v.n(DimensionCount<ArrayL,Scalar>::value);
      }
    };

    IndexContainer<> indexContainer;

    template< int d=0, std::size_t... PI, class ArrayL > // PI:<=>ParameterIndices: needed to transform tupel parameter into a varic list via "std::get<PI>(param)...
       inline void map (ArrayL &&array, std::index_sequence< PI... >* = NULL )
    {
      auto indexLambda = [this]()-> decltype(auto) // is either Index{0} or (Index[]&){indexContainer.v[d]}  => we loop over a plain scalar or use the previously
                         {  if constexpr ( subsettableIndex ) return indexContainer.v[d]=0; else return std::decay_t<Index>{0};  };

      const auto n=size(array);
      if(n>0)
      {
        array[0];// calls either Array<.>::derefer() <= ArrayL inherits Array or gets optimised away.
        for( decltype(auto) i=indexLambda() ;i<n; ++i ) //declt(auto) : replace auto w initialization list & uses rules of decltype => reference detection.
                                                                // REPLACE with iterators ?
          if constexpr( ! std::is_same
                        <
                          std::remove_cv_t< std::remove_reference_t<decltype(array[0])> >,
                          Scalar
                        >::value )
            map< d+1, PI... >(atFast(std::forward<ArrayL>(array),i));
          else
            if constexpr( subsettableIndex )
              invocable( atFast(std::forward<ArrayL>(array),i), indexContainer.v, get<PI> (parameter)... );
            else
              invocable( atFast(std::forward<ArrayL>(array),i), i               , get<PI> (parameter)... );
      }
      //return std::forward<Array>(array);
    }


  public:

    Map(Array &&array, Invocable &&invocable, Parameter &&... parameter)
      : invocable(std::forward<Invocable>(invocable)), parameter(std::forward<Parameter>(parameter)...)
    {
      map(std::forward<Array>(array), static_cast< std::make_index_sequence<std::tuple_size_v<std::tuple<Parameter &&...>>>* >(NULL) );
    }

  };


  template<class Array, class Invocable, class... Parameter>
     inline decltype(auto)
       map (Array &&array,Invocable &&invocable,Parameter &&... parameter)
  // The invocable should have the following signature void (ComponentSType&(&), IntegralType OR Subsettable Integraltype, Parameter...)
  //   whereas Array is resolved into Elements recursively until its ComponentSType equals ComponentSType from Invocable up to cv qualification
  //   Simple one-parameter-Invocables can be mapped directly using the suceeding wrapper.
  // map could be realised with a different-efficient signature: Array map(Array, Invocable&&, Parameters...): for_array(array,i) array[i]=invocable(array[i],i,parameter...);
  //   return value optimisation & parameter replacement theoretically enable copy elusion, st. array=map(array, . , . ); does not result in an array copy.
  //   this ansatz allows might in some cases allow for more copy elusion: pointer-alias-problems.
  //
  // We deal with the constness problem of the indexing operator by the use of atFast(.) and initial call of "array[0];" ( which might perform data copies ("derefer()") for lazy-copy container types or gets optimised away, otherwise.)
  {
    if constexpr( !inheritsArray<Array>::value || std::is_same< remove_cr_t<Array>, typename Map< Array, Invocable, Parameter... >::Scalar >::value )
    {
       typename Map< Array, Invocable, Parameter... >::Index i {0};
       invocable(std::forward<Array>(array),i,parameter...);
    }
    else
      Map< Array, Invocable, Parameter... > map(std::forward<Array>(array),std::forward<Invocable>(invocable),std::forward<Parameter>(parameter)...);

    return std::forward<Array>(array);
  }

  template<class Array, class Invocable, std::enable_if_t< std::tuple_size<get_args_t<Invocable>>::value==1 >* =nullptr >
    inline decltype(auto)
      map (Array &&array,Invocable &&func)
  //Invocables w only one parameter have a different calling signature than Lambda(F,Integral,Parameters...), that is used in the more general "map" above.a
  //  For these simple Invocables, this wrapper of "map" can be invoked directly:
  //   if func has a return value: x=func(x) ; x is an element of array matching the type of f's argument
  //   if func returns void, only func(x) is invoced.
  //     func should than a have a reference typed argument or side effects in order to have a non-trivial map expression
  {
    typedef typename std::tuple_element<0,get_args_t<Invocable>>::type FuncSArg;
    return map( std::forward<Array>(array), [&func](FuncSArg x,Count) // FuncSArg already may contain && or &
                                            {
                                              if constexpr ( std::is_same<void,decltype(func(std::declval<FuncSArg>()))>::value )
                                                func(std::forward<FuncSArg>(x));
                                              else
                                                x=func(std::forward<FuncSArg>(x));
                                            } );
  }


  template<class ComponentSType> Array < ComponentSType > subArray(const Array < Array < ComponentSType > > array,const int  dimensionIdRemove,const Count i)
  {
    if(dimensionIdRemove==0)
      return array[i];
    SubArraySArrayLoop<ComponentSType> arrayLoop(array,dimensionIdRemove,i);

    return arrayLoop.result;
  }


  template<class ComponentSType> Array < ComponentSType > interval(const Array < ComponentSType > array,Count min,Count max)
  {
    if(min<0)  min=0;  if(max>=size(array))  max=size(array)-1;
    Array < ComponentSType > result(max-min+1,NULL);

    for(Count i=0;i<size(result);++i)
      atFast(result,i)=array[i+min];

    return result;
  }


/*  template<class ScalarUser=void,bool local=true,bool boundaryCheck=false,Count d=0,template<class...> class ArrayT,class... ComponentSType, class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,decltype(std::declval<ArrayT<ComponentSType...>>()[0])>>>
    decltype(auto) blockify(ArrayT<ComponentSType...> array, const Array<Count> m)  // transforms a tensor into a double rank tensor by introducing blocks of length m
                                                                                    // m is a scalar or vector of block length.
  {
    if constexpr(is_same_v<Scalar,remove_cr_t<ArrayT<ComponentSType...>>>)
      return std::forward<ArrayT<ComponentSType...>>(array);
    else
    {
      Count b=m[min(d,m.n()-1)];
      nested_t<2*rank<ArrayT<ComponentSType...>,Scalar>,ArrayT,Scalar> r;  // remove_cr_t<decltype(blockify<Scalar>(array[0],0))> r;
        resize(r,size(array)/b);

      for(Count i=size(r);i-->0;)
      {
        auto &&r_i=atFast(r,i);
        resize(r_i,b);
        for(Count j=( !boundaryCheck ? b : min(b,size(array)-i*b) );j-->0;)
          if constexpr(is_same_v<remove_cr_t<decltype(array[0])>,Scalar>)
            atFast(r_i,j)=atFast(array,i*b+j);
          else
            atFast(r_i,j)=blockify<Scalar,boundaryCheck,d+1>(atFast(array,i*b+j),m);
      }
      if constexpr( local&&!is_same_v<remove_cr_t<decltype(array[0])>,Scalar> ) //resort st. that slices from dim d-1 are rearanged (blockwise) wrt to their proximity in dim d:
      {                                                                         //the first step is still
        decltype(r) r_(r.n(),NULL);
          for(Count i)
//could in principle be done simpler by passing a size array - however this would remove compatibility w std::vector and the resize([Array<>]) is more

      }

      return r;
    }
  }
  template<class ScalarUser=void,bool boundaryCheck=false,class ArrayT, class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,ArrayT>>, class Index=Count, std::enable_if_t< !is_subsettable_v<ArrayT> >* =nullptr>
    decltype(auto) blockify(ArrayT array, Count=0)
  {
    if constexpr(is_same_v<Scalar,ArrayT>) // catches wrong scalar specifications
      return std::forward<ArrayT>(array);
  }*/
  template<class ScalarUser=void,template<class...> class ArrayT,class... ComponentSType, class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,decltype(std::declval<ArrayT<ComponentSType...>>()[0])>>,Count rank=rank<ArrayT<ComponentSType...>,Scalar> >
    decltype(auto) blockify(ArrayT<ComponentSType...> array, const Array<Count> m)  // transforms a tensor into a double rank tensor by introducing blocks of length m
                                                                                    // m is a scalar or vector of block length.
  {
    if constexpr(!rank)
      return std::forward<ArrayT<ComponentSType...>>(array);
    else
    {
      Array<Count> size=array.template getSize<Scalar>();
      size.n(rank*2,true);
      for ( Count d=rank;d-->0; )
         {  size[rank+d]=m[d]; size[d]/=m[d];  }

      nested_t<2*rank,ArrayT,Scalar> r;  // remove_cr_t<decltype(blockify<Scalar>(array[0],0))> r;
        resize(r,size);

      map(std::as_const(array), [&r,&m](const Scalar &e,std::array<Count,rank> &i)
      {
        Count blockIndices[rank],intraBlockIndices[rank];
          for(Count d=rank;d-->0;)
          {
                 blockIndices[d]=i[d]/m[d];
            intraBlockIndices[d]=i[d]%m[d];
          }
        atFast<Scalar>(atFast<ArrayT<ComponentSType...>>(r,blockIndices),intraBlockIndices) = e;
      });

      return r;
    }
  }
  template<class ScalarUser=void,class ArrayT, class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,ArrayT>>, class Index=Count, std::enable_if_t< !is_subsettable_v<ArrayT> >* =nullptr>
    decltype(auto) blockify(ArrayT array, Count=0)
  {
    if constexpr(is_same_v<Scalar,ArrayT>) // catches wrong scalar specifications
      return std::forward<ArrayT>(array);
  }

  template<class ScalarUser=void,bool boundaryCheck=true,bool flatten=false,Count d=0,template<class...> class ArrayT, class... ComponentSType, class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,decltype(std::declval<ArrayT<ComponentSType...>>()[0][0])>>>
    decltype(auto) deblockify(ArrayT<ComponentSType...> array,const Array<Count> m=0)
  {
    if constexpr(is_same_v<Scalar,remove_cr_t<ArrayT<ComponentSType...>>>)
      return std::forward<ArrayT<ComponentSType...>>(array);
    else
    {
      nested_t<rank<ArrayT<ComponentSType...>,Scalar>/2,ArrayT,Scalar> r;  // remove_cr_t<decltype(blockify<Scalar>(array[0],0))> r;

      if constexpr (!flatten)
      {
        Count b=0;
        if(m.n()<=d||m[d]==0)
        {
          for(Count i=size(array);i-->0;)
            if(size(atFast(array,i))>b)
              b=size(atFast(array,i));
        }
        else
          b=m[d];
          resize(r,size(array)*b);

        for(Count i=size(array);i-->0;)
        {
          auto &&a_i=atFast(array,i);
          for(Count j=( !boundaryCheck ? b : size(a_i) );j-->0;)
            if constexpr(is_same_v<remove_cr_t<decltype(a_i[0])>,Scalar>)
              atFast(r,i*b+j)=atFast(a_i,j);
            else
              atFast(r,i*b+j)=deblockify<Scalar,boundaryCheck,flatten,d+1>(atFast(a_i,j),m);
        }
      }
      else
      {
        Count n=0;
          for(Count i=size(array);i-->0;)
            n+=size(atFast(array,i)); // if not just every 2nd dim should be eleminated, use AT::map to compute size
          resize(r,n);

        for(Count i=size(array);i-->0;)
        {
          auto &&a_i=atFast(array,i);
          for(Count j=size(a_i);j-->0;)
            if constexpr(is_same_v<remove_cr_t<decltype(a_i[0])>,Scalar>)
              atFast(r,--n)=atFast(a_i,j);
            else
              atFast(r,--n)=deblockify<Scalar,boundaryCheck,flatten,d+1>(atFast(a_i,j));
        }
      }

      return r;
    }
  }
  template<class ScalarUser=void,bool boundaryCheck=true,bool flatten=false,class ArrayT, class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,ArrayT>>, std::enable_if_t< rank<ArrayT><2 >* =nullptr>
    decltype(auto) deblockify(ArrayT&& array,Count=0)
  {
    if constexpr(is_same_v<Scalar,ArrayT>) // catches wrong scalar specifications
      return std::forward<ArrayT>(array);
  }

  template<class ScalarUser=void,class ArrayT> decltype(auto) flatten(ArrayT&& array) // Reduces dimension of a tensor by unrolling subdimensions.
  {  return deblockify<ScalarUser,false,true>(std::forward<ArrayT>(array));  }        // Note that the present version uses deblockify<...>(...) and hence unrolls just every 2nd dimension. For total subdimension elemination, a code-split seems to be a good idea.


  template<bool maximum,bool returnIndex=false, class ScalarUser=void,class SubSettable,class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,decltype(std::declval<SubSettable>()[0])>>>  decltype(auto) extremum(SubSettable &&array)
  {
    typedef std::array<Count,rank<SubSettable,Scalar>> I;
    I r{0};
    Scalar       extrV=std::as_const(array).template at<Scalar>(r);
    map(array, [&r,&extrV](const Scalar &s, const I &i) // later, if opr[] passes rvalues (or: map uses a global at(Array,Index) is used that passes rvalues) Scalar &s might need to become an rvalue if array is one.
    {
      if( (maximum&&s>extrV) || (!maximum&&extrV>s) )
        { extrV=s; for(Count d=0; d< rank<SubSettable,Scalar>; ++d) r[d]=i[d]; }
    } );

    if constexpr (returnIndex)
      if constexpr (rank<SubSettable,Scalar> ==1)
        return r[0]; else return Array<Count>(r);
    else
      return forwardMember<SubSettable>(array.template at<Scalar>(r));
  }

  template<class ComponentSType>  Count maximalComponentId(const Array< ComponentSType > array) //simple but fast.
  {
    Count result=0;
    for(Count i=1;i<size(array);++i)
      if(array[i]>array[result])
    result=i;

    return result;
  }

  template<class ComponentSType>  Count minimalComponentId(const Array< ComponentSType > array)
  {
    Count result=0;
    for(Count i=1;i<size(array);++i)
      if(array[i]<array[result])
    result=i;

    return result;
  }

  template<class ScalarUser=void,class K>  decltype(auto) max(K &&array)
  {  if constexpr (!is_subsettable_v<K>) return array; else return extremum<true ,false,ScalarUser>(std::forward<K>(array));  }
  template<class ScalarUser=void,class K>  decltype(auto) min(K &&array)
  {  if constexpr (!is_subsettable_v<K>) return array; else return extremum<false,false,ScalarUser>(std::forward<K>(array));  }


  template<class ComponentSType>  Array< Array <Count> > histogramClassesSIndexArray(const Array< ComponentSType >& array, Count classN, ComponentSType minimum, ComponentSType maximum) //returns an array of classes whwereas each element is a vector of indices of the records that fall into the class
  {
    Array< Array <Count> > classesSIndexArray(classN,NULL);
    Array <Count> classesSSize(classN,0,NULL);

    struct
    {
      Count get(Count classN,ComponentSType& minimum, ComponentSType& maximum,ComponentSType value)
      {	return Count(double(value-minimum)/double(maximum-minimum)*(double)(classN)); }
    } classId;

    for(Count i=0;i<size(array);++i)
    {
      Count c=classId.get(classN,minimum,maximum,array[i]);
      if(c>=0&&c<classN&&array[i]>=minimum&&array[i]<=maximum)
	      classesSSize[c]++;
    }
    for(Count c=0;c<classN;++c)
      classesSIndexArray[c].resize(classesSSize[c]);
    {
      Array <Count>& classesSFillId=classesSSize; //reuse memory
        classesSFillId=0;

      for(Count i=0;i<size(array);++i)
      {
        Count c=classId.get(classN,minimum,maximum,array[i]);
        if(c>=0&&c<classN&&array[i]>=minimum&&array[i]<=maximum)
        {
      	  classesSIndexArray[c][classesSFillId[c]]=i;
	        ++classesSFillId[c];
        }
      }
    }

    return classesSIndexArray;
  }

  template<class ComponentSType>  Array <Count> histogram(const Array< ComponentSType >& array, Count classN, ComponentSType minimum, ComponentSType maximum)
  {
    const Array<Array<Count>> hCsIA=histogramClassesSIndexArray(array,classN,minimum,maximum);
    Array<Count> r(hCsIA.n(),NULL);
      for_Array(r,c)
        atFast(r,c)=hCsIA[c].n();
    return r;
  }


  template<class ComponentSType>  void quickSortIndexArray(const Array< ComponentSType >& array,Array< Count >& indexArray,Count componentIdMinimal,Count componentIdMaximal)
  {
    std::cout<<"ERROR: template<class ComponentSType> static void quickSortIndexArray(const Array< ComponentSType >& array,Array< Count >& indexArray,Count componentIdMinimal,Count componentIdMaximal) contains a bug in the recursive calls.  -- limits? See. Cf. wikipedia -> quicksort\n";
    ComponentSType pivot   (array[indexArray[componentIdMaximal]]);
    Count          leftId =componentIdMinimal;
    Count          rightId=componentIdMaximal-1;

    std::cout<<"A_1\n";
    do
    {
      std::cout<<"A_2 "<<leftId<<" "<<rightId<<"\n";
      for( ;array[indexArray[leftId]]<=pivot && leftId<componentIdMaximal; ++leftId);
      for( ;array[indexArray[rightId]]>=pivot && rightId>componentIdMinimal; --rightId);

      if(leftId<rightId) //swap ids
      {
        Count tmpId=indexArray[leftId];
        indexArray[leftId]=indexArray[rightId];
        indexArray[rightId]=tmpId;
      }
    }
    while(leftId<rightId);

    if( array[indexArray[leftId]] > pivot )
    {
      Count tmpId=indexArray[leftId];
      indexArray[leftId]=indexArray[componentIdMaximal];
      indexArray[componentIdMaximal]=tmpId;
    }// => all elements on the left are smaller than pivot...

    if(leftId<componentIdMaximal)
      quickSortIndexArray(array,indexArray,leftId,componentIdMaximal); //Error here?
    if(leftId>0&&--leftId>componentIdMinimal)
      quickSortIndexArray(array,indexArray,componentIdMinimal,leftId-1);
  }
  template<class ComponentSType>  Array<Count> quickSortIndexArray(Array< ComponentSType >& array,Count componentIdMinimal=0,Count componentIdMaximal=0)
  {
    if(componentIdMaximal==0)
      componentIdMaximal=size(array)-1;
    if(componentIdMinimal>=componentIdMaximal)
    {
      Array<Count> indexArrayEmpty;
      return indexArrayEmpty;
    }
    Array<Count> indexArray(componentIdMaximal-componentIdMinimal+1,NULL);
    for(Count i=componentIdMinimal;i<=componentIdMaximal;++i)
      indexArray[i]=i;

    quickSortIndexArray(array,indexArray,componentIdMinimal,componentIdMaximal);

    return indexArray;
  }
}


template<class ComponentSType> class MathematicaArray: public Array<ComponentSType>
{
public:
  using Array<ComponentSType>::Array;
  template<class K>MathematicaArray(Array<K> a)
  {
    this->resize(size(a));
    for(Count i=0;i<size(a);++i)
      this->operator[](i)=a[i];  //this->operator=(a); // results in recursive construction.
  }
};

template<class K> MathematicaArray(Array<K>) -> MathematicaArray< nested_t<DimensionCount<K,fundamental_t<K>>::value,MathematicaArray,fundamental_t<K>> >; // deduction guide


template<class ComponentSType> std::ostream& operator<<(std::ostream& ostream,const MathematicaArray<ComponentSType>& array)
{
  ostream<<"{";
  if(array.pData!=NULL)
  {
    for(Count i=0;i<array.size();++i)
    {
                            ostream << (array[i]);
      if(i+1!=array.size()) ostream << ",";
      if(arraySStreamOprOutSWidth>0 && ostream.tellp()<int(i)*(arraySStreamOprOutSWidth+1))
      {
        Count fillN=Count(int(i)*(arraySStreamOprOutSWidth+1)-ostream.tellp());
        for(Count i=0;i<fillN;++i)
          ostream<<" ";
      }
    }
  }
  ostream<<"}";
  return ostream;
}

/* C++11 implemantation without "if constexpr" but SFINAE
     template<class Array, class Lambda, class Parameter...>
       inline static typename enable_if< inheritsArray<Array>::value && !is_same< Array, std::tuple_elements<0,lambda_args<Lambda>::type>::type> >::value > >::value, decltype(auto)>::type
         map (Array&& array,Lambda&& lambda,Parameter... parameter)
 {
      for(auto i=0;i<array.size();++i)
        map(array[i],i,parameter...);
        lambda(array[i],i,parameter...);
    }
*/

    /*struct ArrayLoop//bzw. mapping-template
    {
      virtual void evaluate(Real& sum,Real& squareSum,Real& count,const int& dimensionId) = 0;
      virtual void evaluate(Real& real,const int& dimensionId) = 0;

      void dimensionSLoop(Real& sum,Real& squareSum,Real& count,const int& dimensionId)
      {  evaluate(sum,squareSum,count,dimensionId);  }
      void dimensionSLoop(Real& real,const int dimensionId)
      {  evaluate(real,dimensionId);  }

      template<class ComponentSType>void dimensionSLoop(Array<ComponentSType>& sumArray,Array<ComponentSType>& squareSumArray,Array<ComponentSType>& countArray,const int dimensionId=0)
      {
        for(int i=0;i<sumArray.pData->componentSCount;++i)
          dimensionSLoop(sumArray[i],squareSumArray[i],countArray[i],dimensionId+1);
      }
      template<class ComponentSType>void dimensionSLoop(Array<ComponentSType>& array,const int dimensionId=0)
      {
        for(int i=0;i<array.pData->componentSCount;++i)
          dimensionSLoop(array[i],dimensionId+1);
      }
    };
    typedef ArrayLoop Map;*/
