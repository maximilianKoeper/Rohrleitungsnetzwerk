#include <iostream>
#include <cstdlib>
#include <cmath>
#include "hdnum.hh"
#include <stdexcept>

template<class NumberType>
int flussMatrix( hdnum::DenseMatrix<NumberType> &A )
{
  /** 
  * @brief Erstellt das Gleichungssystem für NxN Knoten
  *
  * Die Matrix \f$ A \f$ wird wie folgt erstellt:\n
  * - Auf der Hauptdiagonalen, also für \f$ \forall i,j \text{ mit } i=j\f$, steht die Summe der
  * Leitwerte aller Zweige, die mit dem Knoten  \f$ i \f$ verbunden sind.\n
  * - An den anderen Stellen  \f$ \forall i,j \text { mit } i \neq j \f$ steht die negative Summe der
  * Leitwerte (hier jeweils 1) zwischen Knoten  \f$ i \f$ und  \f$ j \f$.
  *  
  * @pre Die Matrix A muss quadratisch sein!
  * @pre N > 2
  * @post Matrix A muss symetrisch sein
  * @return 0 bei erfolgreicher Erstellung der Matrix
  * @return 1 falls keine nxn Matrix
  * @return 2 falls Matrix zu klein (n<3)
  *
  */ 
  int M( A.rowsize() );
  int N( A.colsize() );
  if(M != N){
    return 1;
  }
  int n = sqrt(N+1);
    if(n< 3){  // --> double free or corruption (out) Aborted (core dumped)
    return 2;
  }
  int x[n][n];    // Speichert Numerierung der Knoten
  for (std::size_t i=0; i<n; ++i)
  {
    for (std::size_t j=0; j<n; ++j)
    {
      x[i][j] = n*i + j -1;
    }
  }
  
  // Rohleitungsgleichungssystems

  //**********************************************************************
  // Diagonaleinträge auffüllen (Summe der Verbindungen zu anderen Knoten)
  //**********************************************************************
  for (int i = 1; i <= N-1; i++)
  {
    if(i < n){
      A[i-1][i-1] = 3;
    }
    else if(i % n == 0){
      A[i-1][i-1] = 3;
    }
    else if(N-i <= n){
      A[i-1][i-1] = 3;
    }
    else if(i % n == n-1){
      A[i-1][i-1] = 3;
    }
    else{
      A[i-1][i-1] = 4;
    }
  }

  //**********************************************************************
  // Eck-Knoten auffüllen (in unserem Beispiel immer gleich)
  //**********************************************************************
  A[N-1][N-1]=2;
  A[n-2][n-2]=2;
  A[N-n][N-n]=2;

  //**********************************************************************
  // Verbindungen auffüllen (Leitwerte zwischen benachbarten Knoten i,j)
  //**********************************************************************
  for (int i = 0; i < n; i++)
  {
    //***********************
    // vertikale Verbindungen
    //***********************
    if(i not_eq 0 && i not_eq (n-1)){ // Mittlere Zeilen
      for (int j = 0; j < n; j++)
      {
          if(x[i-1][j] > -1){
            A[x[i][j]][x[i-1][j]] = -1;
          }
          A[x[i][j]][x[i+1][j]] = -1;
      }
    }
    else if(i == 0) // Oberste Zeile (j=1 da erster Knoten Referenzknoten)
    {
      for (int j = 1; j < n; j++)
      {
        A[x[i][j]][x[i+1][j]] = -1;
      }    
    }
    else if(i == n-1) // Unterste Zeile
    {
      for (int j = 0; j < n; j++)
      {
        A[x[i][j]][x[i-1][j]] = -1;
      }    
    }

    //*************************
    // horizontale Verbindungen
    //*************************
    if(i == 0){ // 1. Spalte
      for (int j = 1; j < n; j++) // (j=1 da erster Knoten Referenzknoten)
      {
        A[x[j][i]][x[j][i+1]] = -1;
      }      
    }
    else if(i == n-1){ // letzte Spalte
      for (int j = 0; j < n; j++)
      {
        A[x[j][i]][x[j][i-1]] = -1;
      } 
    }
    else if (i not_eq 0 && i not_eq (n-1)){ // Mittlere Spalten
      for (int j = 0; j < n; j++)
      {
        if(x[i-1][j] > -1){
          A[x[j][i]][x[j][i-1]] = -1;
        }
        A[x[j][i]][x[j][i+1]] = -1;
      } 
    }
  }
  return 0;
}


template<class NumberType>
NumberType frobeniusNorm(const hdnum::DenseMatrix<NumberType> &A)
{
  /** 
  *
  * @brief Funktion zur Berechnung der Frobenius-Norm einer Matrix
  * 
  * Diese Funktion berechnet Folgende Norm:
  * \f$||A||_F = \left(\sum_{i,j = 0}^n |a_{ij}|^2\right)^{\frac{1}{2}}\f$
  *
  * @pre  Die Matrix A muss quadratisch sein!
  * @return Gibt berechneten Wert als double zurück
  * @return -1 Falls keine nxn Matrix
  * 
  */ 
  int M(A.rowsize());
  int N(A.colsize());
  if(M != N){
    return -1;
  }

  double result=0.0;

  // Frobenius-Norm
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      result += pow(A[i][j], 2);
    }    
  }
  result = sqrt(result);

  std::cout << "Frobenius-Norm: " << std::setprecision(10) << result << std::endl;
  return result;
}


template<class NumberType>
NumberType maxEigenwert(const hdnum::DenseMatrix<NumberType> &A, int n_iteration_max = 1000, double conv_threshold = 0.00000000001)
{
  /** 
  * @brief Funktion zur Berechnung des betragsgrößten Eigenwertes mit Potenzmethode
  * 
  * - Der Eigenvektor wird wie folgt berechnet 
  * \f$y_{tmp} = \frac{Ay}{||Ay||}\f$ \n
  * - Der dazugehörige Eigenwert: 
  * \f$l = \frac{(y,Ay)_2}{(y,y)_2}\f$\n
  * - Über den Eigenwert wird überprüft ob das Ergebnis konvergiert. Ist \f$l_i-l_{i-1}\f$ < conv_threshold oder i > n_iteration_max wird die
  * Berechnung abgebrochen und das Ergebnis ausgegeben.
  * 
  * @pre  Die Matrix A muss quadratisch sein!
  * @param n_iteration_max die maximale Anzahl der Itartionsschritte (default == 1000)
  * @param conv_threshold Schwelle für die Konvergenz (default == 0.00000000001)
  * @return 0 bei Konvergenz
  * @return 1 falls nicht Konvergiert
  * @return 2 falls keine nxn Matrix
  *
  * @note     Für große N kann die Berechnung auch für kleine n_iteration_max sehr lange dauern
  * @note     Liefert auch ein Ergebnis falls die Berechnung noch nicht Konvergiert ist.
  * @warning  Prüft nicht auf A*y ≠ 0
  */ 
  int M(A.rowsize());
  int N(A.colsize());
  if(M != N){
    return 2;
  }

  // Potenzmethode
  hdnum :: Vector < double > y (N ,1); 
  hdnum :: Vector < double > y_tmp (N);
  double l,l_tmp = 1;
  
  for (int i = 0; i < n_iteration_max; i++)
  {
    // Berechnung laut Aufgabenblatt
    A.mv(y_tmp, y);
    y_tmp /= norm(y_tmp);
    y = y_tmp;

    // Rayleigh-Quotient
    A.mv(y_tmp , y);
    l = (y_tmp*y) / (y*y);

    if(l-l_tmp < conv_threshold && i != 0) // Prüfe auf Konvergenz (Erster Druchgang wird ausgelassen da noch kein l_tmp existiert)
    {    
      std::cout << std::string(50, '-') << std::endl << "Benötigte Iterationsschritte: " << i << std::endl
                << "Betragsmäßig größter Eigenwert: " << l << std::endl
                << "Eigenvektor: " << y << std::endl;
      return 0; 
    }
    l_tmp = l;
  }

  // Nur falls Berechnung nicht konvergiert
  std::cout << std::string(50, '-')
            << "\nACHTUNG Maximale Anzahl der Iteration erreicht! \nNachfolgende Ergebnisse sind eventuell nicht korrekt!" << std::endl
            << std::string(50, '-') << std::endl
            << "Betragsmäßig größter Eigenwert: " << l << std::endl
            << "Eigenvektor: " << y << std::endl;
  return 1;
}



int main(int argc, char ** argv)
{
  /** 
  * @brief Hauptmethode
  *
  * @param argv[1] Falls angegeben Wert für Maximale Iterationsschritte
  * @param argv[2] Falls angegeben Wert für Threshold
  * @return 0 bei erfolgreichem Beenden
  * @return 1 Falls falsche Eingabe
  * @return 2 Falls N < 3
  * @note Stellt automatisch sicher dass Matrix die Größe nxn hat. Alle aufgerufenen Funtionen erwarten eine nxn Matrix!
  *
  */
  std::cout << "Strömung in Rohrleitungsnetzwerken\n" << std::string(50, '-') << std::endl;
  // Anzahl der Knoten
  int N;
  std::cout << "Knotenanzahl N: ";
  std::cin >> N;
  std::cout << std::endl;
  if(std::cin.fail()){  // Eingabe muss vom Typ Int sein
    std::cin.clear();
    std::cin.ignore();
    std::cout << "ERROR: N muss vom Typ Int sein!" << std::endl;
    return 1;
  }
  if(N <= 2){  // Eingabe muss > 2 sein
    std::cout << "ERROR: N muss > 2 sein!" << std::endl;
    return 2;
  }

  // Größe der Matrix
  const int n(N*N-1);

  // Datentyp für die Matrix
  typedef double REAL;

  // Matrix initialisieren (Stellt sicher dass Matrix nxn groß ist!)
  hdnum::DenseMatrix<REAL> A(n,n);

  // Pretty-printing einmal setzen für alle Matrizen
  A.scientific(false);
  A.width(15);
  
  flussMatrix(A);
  if (N<=4){
    std::cout << A << std::endl;
  }

  // Ausgabe der Frobenius-Norm sowie des maximalen Eigenwertes
  frobeniusNorm(A);

  if (argc != 1) // Für angepasste conv_threshold und n_iteration_max Werte.
    {  
      try {
        int arg1 = std::stoi(argv[1], NULL, 10);
        std::cout << std::string(50, '-') << "\nMaximale Iterationen: " << arg1 << std::endl;
        double arg2 = std::stod(argv[2], NULL);
        std::cout << "Konvergenz Threshold: " << std::setprecision(16) << arg2 << std::endl;
        maxEigenwert(A, arg1, arg2);
      }
      catch (...) {
        std::cout << std::string(50, '-') << "\nERROR: Ungültige Argumente in argv[]. \nFallback zu den default Werten" << std::endl;
        maxEigenwert(A);
      }
    }
  else // Für default Werte
    {
        maxEigenwert(A);
    }

  return 0;
}