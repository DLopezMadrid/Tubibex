#include <QtCore/QCoreApplication>
#include <iostream>
#include <QDateTime>

#include "ibex.h"

#include "sivia.h"
#include "init.h"

using namespace std;

int main(int argc, char *argv[]){
    cout<<"Tubibex RC 1"<<endl;

    init_functions();

    Data data;

    //change the variables inside this function to modify the parameters of the program
    config_parameters(data);

    if (data.realTimeDraw){
        init_scene(data, false);
    }

    Sivia siviaOuter, siviaInner;

    QTime t;
    t.start();

    siviaOuter = Sivia(data, false);//execute sivia for the outer approximation of the capture tube

    if (data.calcInner){
        siviaInner = Sivia(data, true);}    //execute sivia for the inner approximation


    double timeSiviaCalculations=t.elapsed()/1000.0;
    cout<<"Elapsed time = "<<timeSiviaCalculations<<endl;
}
