Array library
=============

Hier ist übrigens mal die ganze Library: Falls du mal sehen willst, wo ich immer diese eigenartigen C++-Probleme herbekomme. 
-> siehe dazu array*.h
Eine Anwendung ist zB. linearAlgebra.h. Musst du dir aber nicht genauer ansehen. 

Ganz besonders stolz bin ich aber auf arrayTools.h: ArrayTools::map(Invocable, Param....). Damit lässt sich wunderbar n-Dimensionales programmieren. 
```c++
  int n=0; map ( Tensor, [&n](T& t) { t=n++; } );
```
Läuft beispielsweise über den Tensor "Tensor", wobei die Elemente bis Typ T aufgelöst werden.
```c++
  int n=0; 
  map ( Tensor, [&n](Vec<T>& v) { v=Vec<T>{1,2}*n++; } )
```
macht dann Equivalentes, aber für vektorartige Elemente Vec<T>. Es wird also eine Dimension weniger aufgelöst.

Gpp
===

Dann gibt es noch das bash-Skript gpp:
Es führt c++17-Code aus der Konsole heraus aus (hab es sogar ins /bin/ VZ verlinkt, sodass ich es überall nutzen kann.)
```bash
  > gpp 3.3+7
  10.3
  > gpp "sqrt(9)"
  3
```
Das Skript ist hilfreich um einfache Datenverarbeitung aus der Konsole heraus zu erledigen, für die Bash nicht mehr praktikabel ist -- oder einfach wenn man lieber die c++-Syntax haben möchte.
Ich benutze das regelmäßig.
Hier nun noch ein Beispiel um die Konstante Pi über Monte-Carlo-Integration recht schlecht zu berechnen:
```bash
  ./gpp "pow(DA::mean( Vec<double>(10000000,[&](Count){ return exp(-pow(rand(),2)/2); }) )*L,2)*2"   "double L=99; auto rand=[L](){ return double(std::rand()%200000)/200000*L; }"
```
Oder der Code für das map-Tensor-Beispiel von oben:
```bash
  ./gpp "AT::map( Tensor<3,double>(Vec<int>{3,3,3},NULL), [&n](double &v){ v=n++; }  )" "int n=0"
```
