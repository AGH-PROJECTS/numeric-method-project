#include "metody.h"
#include "ui_metody.h"
#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <iostream>
#include <QMessageBox>
#include <math.h>
#include <QtCharts>
#include <QBrush>
#include <ctime>
#include <cstdio>
#include <chrono>

using namespace std;

class Punkt {
public:
    double x, y;
};

class Wielomian {
public:
    Punkt* tab;
    int n;

    Wielomian(int N) {
        this->n = N;
        tab = new Punkt[n];
    }
};

double Legendre(int k, double x) {
    if (k == 0) {
        return 1;
    }
    if (k == 1) {
        return x;
    }
    double K = (double)k;
    return ((((2 * K) - 1) / K) * x * Legendre(k - 1, x)) - (((K - 1) / K)* Legendre(k - 2, x));
}

double elementOptymalny(Wielomian a, double x, int n) {
    double returnValue = 0.0;

    for (int i = 0; i < n; i++) {
        double licznik = 0.0;
        double mianownik = 0.0;
        for (int j = 0; j < a.n; j++) {
            double legendrei = Legendre(i, a.tab[j].x);
            licznik += a.tab[j].y*legendrei;
            mianownik += legendrei * legendrei;
        }
        returnValue += (licznik / mianownik)*Legendre(i, x);
    }
    return returnValue;

}

double fun(double x) {
    //return (3*x + 3*x*x + 3*x*x*x) * 3*x*x*x*x;
    //return sin(3*x)*2*cos((4/5)*x);
    return (x*(x + 1) - 0.1);
    //return 2*x;
    //return (1 / (x - 3)) - 6;
}

double f(double x)
{
    return (x+3)*(x-1)*(x-1);
   // return (x+3)*(x-1)*(x-1);
}

double brents_fun(double dolna_granica, double gorna_granica, double TOL, double MAX_ITER)
{
    double a = dolna_granica;
    double b = gorna_granica;
    double fa = f(a);   // obliczenie funkcji dla a
    double fb = f(b);   // obliczenie funkcji dla b
    double fs = 0;

    if (!(fa * fb < 0))
    {
        std::cout << "Wartość f(dolna_granica) i f(gorna_granica) muszą być przeciwnych znaków" << std::endl;
        return 0;
    }

    if (std::abs(fa) < std::abs(b)) // jesli wartosc fa jest mniejsza od gornej granicy
    {
        std::swap(a,b);
        std::swap(fa,fb);
    }

    double c = a;
    double fc = fa;
    bool mflag = true;      // uzywam do oceny warunkow
    double s = 0;           // zwracana wartosc
    double d = 0;           // uzywam jesli mglaga jest nieustawiona

    for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
    {

        if (std::abs(b-a) < TOL)
        {
            std::cout << "Po " << iter << " iteracjach pierwiastek to: " << s << std::endl;
            return s;
        }

        if (fa != fc && fb != fc)
        {
            // metoda interpolacji odwrotnej
            s =   ( a * fb * fc / ((fa - fb) * (fa - fc)) )
                + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
                + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
        }
        else
        {
            // metoda siecznych
            s = b - fb * (b - a) / (fb - fa);
        }

        if (    ( (s < (3 * a + b) * 0.25) || (s > b) ) ||
                ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
                ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
                ( mflag && (std::abs(b-c) < TOL) ) ||
                ( !mflag && (std::abs(c-d) < TOL))  )
        {
            //metoda bisekcji
            s = (a+b)*0.5;
            mflag = true;
        }
        else
        {
            mflag = false;
        }

        fs = f(s);
        d = c;
        c = b;
        fc = fb;

        if ( fa * fs < 0)   // fa i fs maja przeciwny znak
        {
            b = s;
            fb = fs;    // jesli tak
        }
        else
        {
            a = s;
            fa = fs;    // jesli nie
        }

        if (std::abs(fa) < std::abs(fb)) // jesli wartosc fa jest mniejsza od wartosc fb
        {
            std::swap(a,b);
            std::swap(fa,fb);
        }

    }

    std::cout<< "Rozwiązanie nie jest zbieżne lub przekroczono ilosc iteracji" << std::endl;

}

metody::metody(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::metody)
{
    ui->setupUi(this);
}

metody::~metody()
{
    delete ui;
}

void metody::on_pushButton_clicked() //obliczenie aproksymacji
{
    vector<double> tabX;
    vector<double> tabY;


    QFile file("aproksymacja.txt");

    if(file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QTextStream in(&file);
        while(!in.atEnd())
        {
            QStringList lineList=in.readLine().split(" ");
            tabX.push_back(lineList[0].toDouble());
            tabY.push_back(lineList[1].toDouble());
        }

    }

        int stopienL=3;
        int i = 0;
        cout.precision(2);

        for (double j = -1; j <= 1; j += 0.125) {
            i++;
        }

        Wielomian W(i);

        i = 0;
        for (double j = -1; j <= 1; j += 0.125) {

            W.tab[i].x = j;
            W.tab[i].y = fun(j);
            i++;
        }

        auto start=std::chrono::steady_clock::now();
      for (double j = -1; j <= 1; j += 0.25)
        {

           elementOptymalny(W, j, stopienL);


            //cout << j << "\t " << fun(j) << "\t " << elementOptymalny(W, j, stopienL) << endl;
        }
      auto end = std::chrono::steady_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
          std::cout << "It took me " << elapsed.count() << " microseconds." << std::endl;
        //wyswietlanie wynikow ->
/*
        for (i = 0; i < W.n; i++) {
            cout << "X" << i << ": " << W.tab[i].x;
            cout << " Y" << i << ": " << W.tab[i].y << endl;
        }

        cout << "wyniki: " << endl;
        cout << "X" << "\t " << "Y" << "\t " << "element optymalny" << endl;
        for (double j = -1; j <= 1; j += 0.25) {
            //elementOptymalny(W, j, stopienL);
          //  cout << j << "\t " << fun(j) << "\t " << elementOptymalny(W, j, stopienL) << endl;
        }
*/
        //rysowanie funkcji ->

        QScatterSeries *series=new QScatterSeries();
        for(double i =-1;i<1;i+=0.25)
        {
            series->append(i,fun(i));
        }
        series->setMarkerSize(10.0);
        series->setName("Punkty aproksymowane");

        QChart *chart = new QChart();
       // chart->legend()->hide();
        chart->addSeries(series);
        chart->createDefaultAxes();
        chart->setAnimationOptions(QChart::AllAnimations);

        QLineSeries *seriesLine=new QLineSeries();
        for(double i =-1;i<1;i+=0.25)
        {
            seriesLine->append(i, elementOptymalny(W,i, stopienL));
        }
        seriesLine->setName("Wielomian st. 3");

        QLineSeries *seriesLine2=new QLineSeries();
        for(double i =-1;i<1;i+=0.25)
        {
            seriesLine2->append(i, elementOptymalny(W,i, 2));
        }
        seriesLine2->setColor(Qt::red);
        seriesLine2->setName("Wielomian st. 2");
        QLineSeries *seriesLine1=new QLineSeries();
        for(double i =-1;i<1;i+=0.25)
        {
            seriesLine1->append(i, elementOptymalny(W,i, 4));
        }
        seriesLine1->setColor(Qt::blue);
        seriesLine1->setName("Wielomian st. 4");

        chart->addSeries(seriesLine);
        //chart->addSeries(seriesLine1);
        //chart->addSeries(seriesLine2);
        QChartView *chartView=ui->tabWidget->widget(0)->findChild<QChartView *>();
        chartView->setChart(chart);

}

void metody::on_pushButton_2_clicked()
{
    double a;               // lower bound
    double b;               // upper bound
    double TOL = 0.0001;    // tolerance
    double MAX_ITER = 1000; // maximum number of iterations

    //przedzialy
    a=-2.999;
    b=-3.0001;

    //test czasu
    auto start=std::chrono::steady_clock::now();
    brents_fun(a,b,TOL,MAX_ITER);
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "It took me " << elapsed.count() << " microseconds." << std::endl;

    QLineSeries *xAxis=new QLineSeries();
    for(double i=-10;i<10;i+=0.25)
    {
        xAxis->append(i,0);
    }
    xAxis->setColor(Qt::black);

    QLineSeries *yAxis=new QLineSeries();
    for(double i=-150;i<200;i+=0.25)
    {
        yAxis->append(0,i);
    }
    yAxis->setColor(Qt::black);

    QScatterSeries *seriesW=new QScatterSeries();
    seriesW->append(brents_fun(a,b,TOL,MAX_ITER),f(brents_fun(a,b,TOL,MAX_ITER)));
    seriesW->setMarkerSize(10.0);
    seriesW->setColor(Qt::red);

    QLineSeries *seriesLine=new QLineSeries();
    for(double i=a-2;i<b+3;i+=0.25)
    {
        seriesLine->append(i, f(i));
    }
    QChart * chartB= new QChart();
    chartB->legend()->hide();
    chartB->addSeries(seriesW);
    chartB->addSeries(xAxis);
    chartB->addSeries(yAxis);
    chartB->addSeries(seriesLine);
    chartB->createDefaultAxes();
    chartB->setAnimationOptions(QChart::AllAnimations);

    chartB->addSeries(seriesW);
    QChartView *chartView=ui->tabWidget->widget(1)->findChild<QChartView *>();
    chartView->setChart(chartB);
}
