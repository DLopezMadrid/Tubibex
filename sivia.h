#ifndef SIVIA_H
#define SIVIA_H

#define __PREC__ 1e-11
#define __METH__ RK4

#include "ibex.h"
#include "init.h"

using namespace ibex;
using namespace std;

class Sivia
{
public:
    Function *f;
    Function *g;

    Sivia(){}                               // Empty constructor

    Sivia(Data &data, bool calcInner);                      // Initializes constraints and prepares data to be processed

    void do_Sivia(Ctc &tubeConstraints, Data &data, Function gdot, bool calcInner);     // Processes the data using contractors and bissections. Classifies the boxes in outside (grey), back_in(yellow) and unsafe (red)

    void print_results(Data data);          // Uses Vibes to represent the results
};

#endif // SIVIA_H
