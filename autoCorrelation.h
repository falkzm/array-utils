//nach Prof. U. Wolf: Monte Carlo errors with less errors
//alle TextVerweise beziehen sich auf obige Literatur

class AutoCorrelation
{//nach U. Wolf: Monte Carlo errors with less errors
public:
  Count n;//#Values in observablesSSetsSValues
  Vector< Vector<Real> > observablesSSetsSMeanValue;
  Vector<Real> observablesSMeanValues;

  Vector< Vector<Real> > fProjectedSetsSDeviancesToMeanValue;//f ist die auszuwertende (abgeleitete) Observable. Die Projektion meint projektion der ObservablenDimension auf den GardientVonF, wobei der Gradient nun selbst durch ableitung von f nach den Obversablen geschehen ist. Die Projektion lässt also die ObservablenAbhaengigkeit verschwinden. 

  Real fSMeanValue,fSError;
  Vector<Real> setsSFSMeanValue;

  Vector<Real> correlationLengthsSAutoCorrelation;
  Vector<Real> windowsSIntegratedAutoCorrelationTime;
  Real integratedAutoCorrelationTime;


  template<class ComponentSType> ComponentSType calculateMeanValue(Vector<ComponentSType> values)
  {
    ComponentSType sumComponentSType=0;
    Real normalisationReal=1/Real(values.pData->componentSCount);

    for(int i=0;i<values.pData->componentSCount;++i)
      sumComponentSType+=values[i]*normalisationReal;

    return sumComponentSType;
  }


  Real componentSSum(Real real) {return real;}
  template<class ComponentSType> Real componentSSum(Array<ComponentSType> array,int componentIdMinimal,int componentIdMaximal)
  {
    cout<<"6.1\n";
    Real result=0;

    cout<<"6.2\n";
    for(int i=componentIdMinimal;i<array.pData->componentSCount&&i<componentIdMaximal;++i)
      result=componentSSum(array[i]);//so werden auch unterVektoren einbezogen, per Ueberladung wird das RekursionsEnde ggf. ausgewaehlt.

    cout<<"6.3\n";
    return result;
  }
  template<class ComponentSType> Real componentSSum(Array<ComponentSType> array)
  {
    return componentSSum(array,0,array.pData->componentSCount);
  }


  virtual Real f/*derived Observable*/(Vector<Real> observablesSValues)
  {
    cout<<"7\n";
    return observablesSValues[0];
  }

private://da zu viele SeitenEffekte

  void calculateFProjectedSetsSDeviancesToMeanValue(Vector< Vector < Vector<Real> > > observablesSSetsSValues,int evaluationObservableId)
  {//based on n, means, fSError
    cout<<"3.1\n";
    if(evaluationObservableId>=0)
      fProjectedSetsSDeviancesToMeanValue=observablesSSetsSValues[evaluationObservableId]-observablesSMeanValues[evaluationObservableId];
    else
    {
      cout<<"3.2\n";
      Vector<Real> fSGradient(observablesSSetsSValues.pData->componentSCount,NULL);

      for(int o=0;o<observablesSSetsSValues.pData->componentSCount;++o)
      {
        cout<<"3.3\n";
        Vector< Vector<Real> > vectorTmp=observablesSSetsSValues[o]-observablesSMeanValues[o];//Abweichung der SetsValues vom MittelWert
        cout<<"3.3.1\n";
        vectorTmp*=vectorTmp;//komponentenWeise multiplikation
        cout<<"3.3.2\n";
        Real fSGradientSMeasureDistanceReal=sqrt(componentSSum(vectorTmp))/Real(n);//wie bei StandardAbweichung
        cout<<"3.4\n";

        if(fSGradientSMeasureDistanceReal==0)
          fSGradient[o]=0;//da fSGradientSMeasureDistance die StdAbweichung ist, fluktuiert f für dieses o nicht
        else
        {//DifferenzenQuotient...
          cout<<"3.5\n";
          observablesSMeanValues[o]+=fSGradientSMeasureDistanceReal;
          fSGradient[o]=f(observablesSMeanValues);
          observablesSMeanValues[o]-=2*fSGradientSMeasureDistanceReal;
          fSGradient[o]-=f(observablesSMeanValues);
          observablesSMeanValues[o]+=fSGradientSMeasureDistanceReal;//AusgangsZustand wiederhergestellt.
          fSGradient[o]/=fSGradientSMeasureDistanceReal;
        }
        cout<<"3.6\n";
      }

      cout<<"3.7\n";
      fProjectedSetsSDeviancesToMeanValue=(observablesSSetsSValues-observablesSMeanValues)*fSGradient;
    }
    cout<<"3.8\n";
  }


  Count setsSValuesSCountMinimal(Vector< Vector< Vector<Real> > > observablesSSetsSValues)
  {
    Count resultCount=observablesSSetsSValues[0][0].pData->componentSCount;

    int o=0;//Annahme: alle Observablen wurden gleichoft vermessen.
      for(int s=0;s<observablesSSetsSValues[o].pData->componentSCount;++s)
        if(observablesSSetsSValues[o][s].pData->componentSCount<resultCount)
          resultCount=observablesSSetsSValues[o][s].pData->componentSCount;

    return resultCount;
  }

  void calcFErrosAndAutoCorrelations(Vector< Vector < Vector<Real> > > observablesSSetsSValues,Real autoCorellationTimesSPropornationalityConstantReal)
  {//based on n, means, fprojected deviances
    int windowMaximal=setsSValuesSCountMinimal(observablesSSetsSValues)/2; //window meint hier die # der aufeinanderfolgende Messwerte, ueber die die AutoKorrelation ausgewertet wird.
    int windowOptimal;//der Wert, wo die systematischen Fehler auf Grund des begrenzten Fesnters + die statistischen minimal sind. (Vgl. Paper v. Wolf , S.10)
    Vector<Real> windowsSFProjectedAutoCorrelation(windowMaximal+1,NULL);
    Real fProjectedIntegratedAutoCorrelationTimeUnnormalized=1;//der window=0 - Anteil
    bool aborted=true;

    for(int window=0;window<=windowMaximal;++window)
    {///eigentlich hier ist window hier vielmehr die Spanne t mit der die AutoKorrelationsFunktion aufgerufen wird
      //windowsSFProjectedAutoCorrelation ausrechnen:
      Real& windowSFProjectedAutoCorrelation=windowsSFProjectedAutoCorrelation[window];//eigentlich hier ist window hier vielmehr die Spanne t mit der die AutoKorrelationsFunktion aufgerufen wird
      {
        const int o=0;//alle Observablen wurden gleichoft vermessen

        windowSFProjectedAutoCorrelation=0;//es folgt (31) mit (33) - die Projektion wurde schon vorher eingebaut
        for(int s=0;s<observablesSSetsSValues[o].pData->componentSCount;++s)
        {
          for(int v=0;v<observablesSSetsSValues[o][s].pData->componentSCount-window;++v)
            windowSFProjectedAutoCorrelation+=fProjectedSetsSDeviancesToMeanValue[s][v]*fProjectedSetsSDeviancesToMeanValue[s][v+window];
        }
        windowSFProjectedAutoCorrelation/=Real(n-observablesSSetsSValues[o].pData->componentSCount*window);
      }

      //feststellen, ob dies das "optimale" window ist
      if(window>0)
      {
        Real autoCorrelationTime=0.000001;//"tiny Value", falls nachfolgendes fehlschlägt: Annahme: autoCorrelationTime=autoCorellationTimesSPropornationalityConstantReal*fProjectedIntegratedAutoCorrelationTime //Vgl U.Wolff S.12

        fProjectedIntegratedAutoCorrelationTimeUnnormalized+=2*windowSFProjectedAutoCorrelation/windowsSFProjectedAutoCorrelation[0];//die AutocorrelationsZeit ergibt sich eben aus der Summe ueber alle AutoKorrelationen und der Normalisierung; nach (35)
        if(fProjectedIntegratedAutoCorrelationTimeUnnormalized>1)
          autoCorrelationTime=autoCorellationTimesSPropornationalityConstantReal/log( (fProjectedIntegratedAutoCorrelationTimeUnnormalized+1)/(fProjectedIntegratedAutoCorrelationTimeUnnormalized-1) );

        const Real autoCorrelationErrorDerivationUnnormalized=exp(-window/autoCorrelationTime)-autoCorrelationTime/sqrt(Real(window*n));//1. Term systematischer Fehler der AutokorrelationsZeit wegen FensterBeschneidung, 2. statistischer Fehler. ; Vgl U.Wolff S.12
        cout<<"suppenSchildKröte = "<<autoCorrelationErrorDerivationUnnormalized<<"\nWindowMaximal = "<<windowMaximal<<"\n";
        if(autoCorrelationErrorDerivationUnnormalized<0)//VorzeichenWechsel -> 1.Ableitung hier ~ gesetzt => Minimum des Fehlers gefunden => wOptimal gefunden
        {
          windowOptimal=window;
          windowMaximal=( 2*window<windowMaximal ? 2*window : windowMaximal ); //nicht mehr so lange suchen
          aborted=false;
        }
      }
    }
    windowsSFProjectedAutoCorrelation.resize(windowMaximal+1,true);

    //primaere Auswertung
    if(aborted)
    {
      printf("WARNING: no window with error minimisation found up to window=%i\n",windowMaximal);
      fProjectedIntegratedAutoCorrelationTimeUnnormalized*=windowsSFProjectedAutoCorrelation[0];//das haben wir vorhin nur hineindividiert;fProjectedIntegratedAutoCorrelationTimeUnnormalized==CFbb (opt) nach (41) nach (35) mit (33) nach (25)
    }
    else
      fProjectedIntegratedAutoCorrelationTimeUnnormalized=windowsSFProjectedAutoCorrelation[0]+2*componentSSum(windowsSFProjectedAutoCorrelation,1,windowOptimal);//def CFbb;&& fProjectedIntegratedAutoCorrelationTimeUnnormalized==CFbb (opt) nach (41) nach (35) mit (33) nach (25) (wie oben)
    //U.Wolff: GammaFbb=GammaFbb+CFbbopt/N;                       % bias in Gamma corrected
    windowsSFProjectedAutoCorrelation+=fProjectedIntegratedAutoCorrelationTimeUnnormalized/Real(n); //Vgl Seite 8 unten oder (49)
    //U.Wolff: CFbbopt=GammaFbb(1) + 2*sum(GammaFbb(2:Wopt+1));   % refined estimate
    fProjectedIntegratedAutoCorrelationTimeUnnormalized=windowsSFProjectedAutoCorrelation[0]+2*componentSSum(windowsSFProjectedAutoCorrelation,1,windowOptimal);//wie oben - wird nur mit bias-korrigiertem windowSFProjectedAutoCorrelation neuberechnet;//hier kein wOpt+1 weil wOpt im Wollf-Alg in zusammenhang mit t+1-Groessen gefunden wurde
    //U.Wolff: sigmaF=sqrt(CFbbopt/N);                            % error of F
    fSError=sqrt(fProjectedIntegratedAutoCorrelationTimeUnnormalized/Real(n));//Vgl (21) bzw. (44)
    //U.Wolff: rho=GammaFbb/GammaFbb(1);                          % normalized autocorr.
    correlationLengthsSAutoCorrelation=windowsSFProjectedAutoCorrelation*(1/windowsSFProjectedAutoCorrelation[0]);//fast (41) - es fehlt 1/2; zudem vür verschiedene KorellationsLaengen t=="window"
    //U.Wolff: tauintFbb=cumsum(rho)-0.5;
    windowsSIntegratedAutoCorrelationTime.resize(correlationLengthsSAutoCorrelation.pData->componentSCount);
    windowsSIntegratedAutoCorrelationTime[0]=correlationLengthsSAutoCorrelation[0];
    for(int t=1;t<correlationLengthsSAutoCorrelation.pData->componentSCount;++t)
      windowsSIntegratedAutoCorrelationTime[t]=windowsSIntegratedAutoCorrelationTime[t-1]+correlationLengthsSAutoCorrelation[t];//cumulativ sum
    windowsSIntegratedAutoCorrelationTime-=0.5;//ähnlich wie fProjectedIntegratedAutoCorrelationTimeUnnormalized; nur werden hier NORMALISIERTE AutoKorrelationen aller laengen gleichoft+einfach gezählt und die t=0te Korrelation geht haber nur halbsooft ein (sum(kor(t) von t=-inf bis +inf =>Korr(0) einmal enthalten) und  die 0. eben zu 1 normalisiert wurde, muss nur 0.5 abgezogen werden =1*Gewicht=0.5  ( hier wird nicht wie oben *2 gerechnet)
      //tauint  = tauintFbb(Wopt+1);
      integratedAutoCorrelationTime=windowsSIntegratedAutoCorrelationTime[windowOptimal];//hier kein wOpt+1 weil wOpt im Wollf-Alg in zusammenhang mit t+1-Groessen gefunden wurde
    //% changed subleading terms (prev. versions!) here:
    //dtauint = tauint*2*sqrt((Wopt-tauint+0.5)/N);
    //sigmaSigmaF = sigmaF*sqrt((Wopt+0.5)/N);

  }


  void calcFMeans(Vector< Vector< Vector<Real> > > observablesSSetsSValues,int evaluationObservableId)
  {
    cout<<"5.1\n";
    if(evaluationObservableId>0)
    {
      cout<<"5.2\n";
      fSMeanValue     =observablesSMeanValues    [evaluationObservableId];
      setsSFSMeanValue=observablesSSetsSMeanValue[evaluationObservableId];
    }
    else
    {
      cout<<"5.3\n";
      fSMeanValue=f(observablesSMeanValues);

      setsSFSMeanValue.resize(observablesSSetsSValues[0].pData->componentSCount);
      {
        cout<<"5.4\n";
        Vector<Real> setSObservablesSMeanValues(observablesSSetsSMeanValue.pData->componentSCount,NULL); 

        for(int s=0;s<observablesSSetsSValues[0].pData->componentSCount;++s)
        {
          cout<<"5.5\n";
          for(int o=0;o<observablesSSetsSValues.pData->componentSCount;++o)
            setSObservablesSMeanValues[o]=observablesSSetsSMeanValue[o][s];
          setsSFSMeanValue[s]=f(setSObservablesSMeanValues);
          cout<<"5.6\n";
        }
      }
    }
    cout<<"5.7\n";

    if(observablesSSetsSValues[0].pData->componentSCount>1)//Bias von fSMean nur behandelbar, falls mehere DatenSets je zu einer Observable vorliegen
    {
      cout<<"5.8\n";
      int   o=0;
      Real  fSMeanValueBySetsSFSMeanValue=0;
      Count nSets;
      {
        Real normalisationReal=1/Real(n);

        nSets=observablesSSetsSValues[o].pData->componentSCount;
        for(int s=0;s<nSets;++s)
          fSMeanValueBySetsSFSMeanValue+=setsSFSMeanValue[s]*Real(observablesSSetsSValues[o][s].pData->componentSCount)*normalisationReal;
      }
      cout<<"5.9\n";

      //U. Wolff: bias cancellation for the mean value
      Real bias=(fSMeanValueBySetsSFSMeanValue-fSMeanValue)/(nSets-1);
      fSMeanValue-=bias;
      if(absReal(bias)>fSError/4)
        cout<<"WARNING: removed bias="<<bias<<" > Error/4="<<fSError/4<<"\n";
      cout<<"5.10\n";

      for(int s=0;s<nSets;++s)
        setsSFSMeanValue[s]-=bias*Real(n)/Real(observablesSSetsSValues[o][s].pData->componentSCount);
    }
    cout<<"5.11\n";
    /*U. Wollf:
      if R >= 2
        bF = (Fb-Fbb)/(R-1);
        Fbb=Fbb - bF;
        if abs(bF) > sigmaF/4
          beep
          fprintf('\n WARNING: a %.1f sigma bias of the mean has been cancelled \n',bF/sigmaF)
        end
        Fbr = Fbr - bF*N./Nrep;
        Fb  = Fb  - bF*R;
      end
    */
  }

public:

  void calc(Vector< Vector < Vector<Real> > > observablesSSetsSValues/*observablesVectorSSetsVectorSValuesVector - *Vector weglassen -> ergibt FeldTypNotation*/,int evaluationObservableId=-1/*<0 heisst  es wird die agbeleitet ObservablenFkt betrachtet*/,Real autoCorellationTimesSPropornationalityConstantReal=1.5)
  {
    cout<<"1\n";
    observablesSSetsSMeanValue.resize(observablesSSetsSValues.pData->componentSCount);//Dim o,s
    {
      for(int o=0;o<observablesSSetsSMeanValue.pData->componentSCount;++o)
      {
        observablesSSetsSMeanValue[o].resize(observablesSSetsSValues[o].pData->componentSCount);//einstellen der setAnzahl an der Stelle Observable=o -- eigentlich eine konstante uerber alle o
        for(int s=0;s<observablesSSetsSMeanValue[o].pData->componentSCount;++s)
          observablesSSetsSMeanValue[o][s]=calculateMeanValue(observablesSSetsSValues[o][s]);
      }
    }

    cout<<"2\n";
    observablesSMeanValues.resize(observablesSSetsSMeanValue.pData->componentSCount);
    n=0;//#Values in observablesSSetsSValues
    {
      for(int s=0;s<observablesSSetsSValues[0].pData->componentSCount;++s)//Bedingung: jede Observable wurde in gleichmächtigen Sets vermessen.
        n+=observablesSSetsSValues[0][s].pData->componentSCount;
      Real normalisationReal=1/Real(n);
      n*=observablesSSetsSValues.pData->componentSCount;

      for(int o=0;o<observablesSSetsSMeanValue.pData->componentSCount;++o)
      {
        observablesSMeanValues[o]=observablesSSetsSMeanValue[o][0]*(Real(observablesSSetsSValues[o][0].pData->componentSCount)*normalisationReal);
        for(int s=1;s<observablesSSetsSMeanValue[o].pData->componentSCount;++s)
          observablesSMeanValues[o]=observablesSSetsSMeanValue[o][s]*(Real(observablesSSetsSValues[o][s].pData->componentSCount)*normalisationReal);
      }
    }

    cout<<"3\n";
    calculateFProjectedSetsSDeviancesToMeanValue(observablesSSetsSValues,evaluationObservableId);
    cout<<"4\n";
    calcFErrosAndAutoCorrelations(observablesSSetsSValues,autoCorellationTimesSPropornationalityConstantReal);
    cout<<"5\n";
    calcFMeans(observablesSSetsSValues,evaluationObservableId);
  }


  void print(bool printObservablesSMeanValues=true,bool printSetsSFSMeanValue=true,bool printCorrelationLengthsSAutoCorrelation=true,bool printWindowsSIntegratedAutoCorrelationTime=true)
  {
    cout<<"=== gamma-method ===\n\n";

    cout<<"  #Values n="<<n<<"\n";//#Values in observablesSSetsSValues
    if(printObservablesSMeanValues)
    {
      cout<<"  observablesSMeanValues:\n";
      for(int o=0;o<observablesSMeanValues.pData->componentSCount;++o)
          cout<<"    o="<<o<<": "<<observablesSMeanValues[o]<<"\n";
    }
    cout<<"\n";

     Vector< Vector<Real> > fProjectedSetsSDeviancesToMeanValue;//f ist die auszuwertende (abgeleitete) Observable. Die Projektion meint projektion der ObservablenDimension auf den GardientVonF, wobei der Gradient nun selbst durch ableitung von f nach den Obversablen geschehen ist. Die Projektion lässt also die ObservablenAbhaengigkeit verschwinden. 
    cout<<"  derived Observable f:\n";
    cout<<"    meanValue:"<<fSMeanValue<<"\n";
    cout<<"    error:"<<fSError<<"\n";
    if(printSetsSFSMeanValue)
    {
      cout<<"    setsSMeanValue:\n";
      for(int s=0;s<setsSFSMeanValue.pData->componentSCount;++s)
       cout<<"      s="<<s<<": "<<setsSFSMeanValue[s]<<"\n";
    }
    cout<<"\n";

    cout<<"  autoCorrelation:\n";
    if(printCorrelationLengthsSAutoCorrelation)
    {
      cout<<"    correlationLengthsSAutoCorrelation:\n";
      for(int t=0;t<correlationLengthsSAutoCorrelation.pData->componentSCount;++t)
       cout<<"      t="<<t<<": "<<correlationLengthsSAutoCorrelation[t]<<"\n";
    }
    if(printWindowsSIntegratedAutoCorrelationTime)
    {
      cout<<"    windowsSIntegratedAutoCorrelationTime:\n";
      for(int w=0;w<windowsSIntegratedAutoCorrelationTime.pData->componentSCount;++w)
       cout<<"      w="<<w<<": "<<windowsSIntegratedAutoCorrelationTime[w]<<"\n";
    }
    cout<<"    integratedAutoCorrelationTime(from optimal Window)="<<integratedAutoCorrelationTime<<"\n";
    cout<<"\n";

    cout<<"\n";
  }
};
