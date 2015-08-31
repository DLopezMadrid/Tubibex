#ifndef INIT_H
#define INIT_H

#include "ibex.h"
#include <QList>
#include <fstream>
#include <string>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

//for the sleep command in windows
//#include <windows.h>


using namespace ibex;
using namespace std;

struct Data{
    IntervalVector initialBox;
    QList<IntervalVector> boxes;
    QList<IntervalVector> boxesOutside;
    QList<IntervalVector> boxesUnsafe;
    QList<IntervalVector> boxesInside;
    QList<IntervalVector> boxesInsideUnsafe;
    QList<IntervalVector> boxesInsideBackIn;
    QList<IntervalVector> boxesBackIn;
    QList<IntervalVector> boxesUnsafeFuture;
    QList<IntervalVector> boxesBackInFuture;

    QList<double> errorMax;
    vector<double> epsilons;

    double diskR;
    bool showDisk;
    double imageSize;

    double dt;
    double tMin;
    double tMax;


    double delayMs;
    bool delayActivated;
    bool realTimeDraw;
    bool saveAllPossibleFigures;

    int maxNumUnsafeBoxes;
    bool maxNumUnsafeBoxesActivated;
    int numFuturePos;
    bool calcInner;
    bool drawFutureBoxesUnion;
    bool enableBackIn;

    Function *f, *g;

    int numVarF;
    int numVarG;

    double gMax;

    int plane;
    int intMethod;

    double VFepsilonX;
    double VFepsilonY;

    bool myDebug;

};

bool file_exists(const string& filename);                       // Using stat, checks if there is already a file with the current filename

bool init_functions();                                          // Creates f.txt and g.txt in the build directory with the default functions

void config_parameters(Data &data);                             // Contains all the user configurable parameters and initializes some of the variables of data s.t. the initial box

void init_scene(Data data, bool showFuturePos);                 // Initializes the scene in vibes for a certain plane

void draw_box(IntervalVector myInt, int color, int dim);        // Draws a box in vibes in a given plane and using a given color

void draw_update(Data data, IntervalVector box, IntervalVector contractedBox);      // Using VIBes, draws the current boxes stored in data excluding the future ones. A delay can be applied to see how sivia works in real time

void draw_update(Data data, bool showFuturePos);                // Same as previous one but this one also includes the future boxes and removes the time delay

void draw_vector_field(Data data);                              // Draws the general vector field

void draw_vector_field_boxes(QList<IntervalVector> posList, int color);     // Draws the vector field on the center of the provided list of boxes

void draw_vector_field_boxes(QList<IntervalVector> posList);                // Same as the previous one but draws the vectors in the default color

void draw_arrow(vector<double> arrowX, vector<double> arrowY, int color);           // Calculate the correspondient vectors

void calculate_vector(IntervalVector& box, double t, int color);            // Draw arrow using a line and a box from the given coordinates

void export_all_images(Data& data);








#endif // INIT_H
