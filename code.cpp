#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <random>
#include <fstream>
using namespace std;
double a,b; 
void weightedLinearRegression(const vector<double>& x, const vector<double>& y, const vector<double>& sigma_y, double& a, double& b) {
    // Verifica che le dimensioni dei vettori siano uguali
    if (x.size() != y.size() || x.size() != sigma_y.size()) {
        cerr << "Errore: i vettori devono avere la stessa dimensione!" << endl;
        return;
    }

    int n = x.size();
    double sum_w = 0.0;
    double sum_wx = 0.0;
    double sum_wy = 0.0;
    double sum_wxy = 0.0;
    double sum_wx2 = 0.0;

    // Calcolo delle somme pesate
    for (int i = 0; i < n; ++i) {
        double w = 1.0 / (sigma_y[i] * sigma_y[i]); // Peso inversamente proporzionale al quadrato dell'errore
        sum_w += w;
        sum_wx += w * x[i];
        sum_wy += w * y[i];
        sum_wxy += w * x[i] * y[i];
        sum_wx2 += w * x[i] * x[i];
    }

    // Calcolo dei parametri della retta di regressione
    double denominator = sum_w * sum_wx2 - sum_wx * sum_wx;
    if (denominator == 0) {
        cerr << "Errore: divisione per zero nel calcolo dei parametri della regressione!" << endl;
        return;
    }

    b = (sum_w * sum_wxy - sum_wx * sum_wy) / denominator;
    a = (sum_wx2 * sum_wy - sum_wx * sum_wxy) / denominator;
}
int main ()
{
    ofstream out("isto.txt"); 

    normal_distribution<double> n(0,1); 
    default_random_engine e; // generatore pseudo-casuale 
    random_device rd;
    e.seed(rd()); 

    int N = 0;
    double q = 0;  

    cout << "Quante iterazioni vuoi? "; 
    cin >> q; 
    cout << "Inserire numero di GDL (N-2): "; 
    cin >> N; 

    vector <double> y (N,0); // vettore contenente N numeri che stanno sulla retta y = 3x + 2 per esempio  
    vector <double> ey (N,0); // errori gaussiani sulle y 
    vector <double> x (N,0); // var ind. 
    

    vector<double> chiq (q,0); 

    for (int z = 0; z < q; z ++)
    {

        double t = 0; 
        for (int i = 0; i < N; i ++ )
        {
            x.at(i) = i; 
            y.at(i) = (3*x.at(i)) + 2; // adesso le y e le x seguono esattamente il grafico 
            ey.at(i) = 1.0; // errore standard noto (costante o variabile, ma mai negativo o vicino a zero)
            t = n(e) * ey.at(i); // fluttuazione statistica
            y.at(i) += t; // aggiungo il rumore gaussiano

        }

        // adesso abbiamo i dati "sperimentali"

        // facciamo la regressionen lineare per stimare i parametri a e b 

        weightedLinearRegression(x, y, ey, a, b);

        // calcolo del chi-quadro 

        for (int i = 0; i < N; i ++ )
        {
            chiq.at(z) += (pow(y.at(i) - (a+(b*x.at(i))),2))/(pow(ey.at(i),2));
        }
    }

    // trovo quindi la matrice dei chi quadri

    // divido i chi quadri in bin di 0,1

    double bin = 0.1; 

    int num = 0; 
    double max = 0; 
    double min = 0; 
    max = chiq.at(0); 
    min = chiq.at(0); 
    for (int i =0; i < q; i ++ )
    {
        if (max < chiq.at(i))
            max = chiq.at(i); 

        if (min > chiq.at(i))
            min = chiq.at(i);
    }

    num = (max - min)/bin; 

    vector <double> cont (num,0); // conteggi per ogni chiq  

    double h = 0; // counter 

    for (int i = 0; i < q; i ++ ) // passiamo tutti i chi q
    {
        for (int g = 0; g < num; g ++) // passiamo tutti i bin 
        {
            if ( (chiq.at(i) >= (min + (bin*g))) && (chiq.at(i) < (min + (bin*(g+1)))) )
                cont.at(g) += 1;
            // aggiungo uno se il chi quadro è all'interno del bin 
        }
    } 

    // aggiungiamo all'interno di un file 

    for (int i = 0; i < num; i ++ )
    {
        out << (min + (i*bin) + (min + ((i+1)*bin)))/2 << " " << cont.at(i)/(q*bin) << "\n "; 
    }


    return 0; 
}

// su gnuplot 
// set grid 
// set xlabel "chiq"
// set ylabel "f(x)"
// k = GDL 
// f(x) = (x > 0) ? (1.0 / (2**(k/2.0) * gamma(k/2.0))) * x**(k/2.0 - 1) * exp(-x/2.0) : 0
// plot f(x) title "Chi-squared PDF (k=100)", "isto.txt" using 1:2 with points ps 0.1 title "Simulazione"


