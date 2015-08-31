#include "init.h"
#include "vibes.h"
#include "vibes.cpp"



/// Contains all the user configurable parameters and initializes some of the variables of data s.t. the initial box
void config_parameters(Data &data){

    Variable y(data.numVarF),t;

    //start of configurable parameters
    data.diskR=1;           //draw a circle of radius diskR
    data.showDisk=true;    //show the previous circle

    data.dt=0.1;            //time step size
    data.tMin=0;            // start time
    data.tMax=0.1;          // end time

    data.maxNumUnsafeBoxes=100;                 //max num of unsafe boxes that is needed to be reached to stop the loop and show results if the next variable is true
    data.maxNumUnsafeBoxesActivated=false;      //to activate the counting of unsafe boxes
    data.numFuturePos=3;                        //number of future positions taken into account by euler method and guaranteed integration
    data.calcInner=true;                       //calculate also inner approximation
    data.drawFutureBoxesUnion=false;             //draw the future boxes as a union or each one separately
    data.enableBackIn=false;                    //try to use euler integration to determine if the boxes come back to the tube in a given time t

    //for the pendulum
    double myx[2][3]={
        {-3, 3, 0.1} ,                 //x[0] min, max, epsilon
        {-3, 3, 0.1}                  //x[1] min, max, epsilon
    };

    //non holonomic vehicle
    //    double myx[3][3]={
    //        {-2, 2, 0.05} ,                 //x[0] min, max, epsilon
    //        {-2, 2, 0.05} ,                 //x[1] min, max, epsilon
    //        {-1, 1, 0.05}
    //    };


    data.realTimeDraw=false;            //draw the paving in real time
    data.delayActivated=true;           //add a delay to the real time drawing to see it better
    data.delayMs=100;                    //how big is the delay in ms

    data.plane=0;                       //  0=XY  1=XZ   2=YZ

    data.VFepsilonX=0.25;                //mesh precision for the vector field for X axis
    data.VFepsilonY=0.25;                //same for Y axis

    data.intMethod=1;                   // 0=guaranteed integration, 1=euler method

    data.myDebug=false;                  //include cout statements with some useful information

    data.saveAllPossibleFigures=true;         //saves a png file for every possible combination of inner approx, outer approx and vector field
    data.imageSize=800;                         //image export size X by X pixels

    //end of configurable parameters




    //store functions f and g from the text files
    data.f = new Function ("f.txt");
    data.g = new Function ("g.txt");

    //extract the number of variables for f and g
    data.numVarF=data.f->nb_var()-1;
    data.numVarG=data.g->nb_var()-1;

    //initialize errorMax
    for(int i = 0; i<data.numVarG+1; i++){
        data.errorMax.push_back(0);
    }

    //Define the initial box
    double auxBox[data.numVarG+1][2];

    for (int i=0;i<data.numVarG;i++) {
        auxBox[i][0]=myx[i][0];        //to simulate a box from -inf to inf
        auxBox[i][1]=myx[i][1];
        data.epsilons.push_back(myx[i][2]);
    }

    auxBox[data.numVarG][0]=data.tMin;
    auxBox[data.numVarG][1]=data.tMax;

    IntervalVector myBox(data.numVarG+1, auxBox);

    data.initialBox=myBox;

}



/// Creates f.txt and g.txt in the build directory with the default functions
bool init_functions(){

    if (!file_exists("f.txt")){
        std::string input;
        std::string line1, line2, line3, line4, line5, line6;
        line1="function f(x[2],t)\n";
        line2="         return(\n";
        line3="             x[1],\n";
        line4="             -sin(x[0])-0.15*x[1]\n";
        line5="         );\n";
        line6="     end\n";
        input=line1+line2+line3+line4+line5+line6;
        std::ofstream out("f.txt");
        out << input;
        out.close();
        cout<<"f.txt created"<<endl;
    }

    if (!file_exists("g.txt")){
        std::string gInput;
        std::string gLine1, gLine2, gLine3, gLine4, gLine5;
        gLine1="function g(x[2],t)\n";
        gLine2="         return(\n";
        gLine3="                 ((x[0])^2+(x[1])^2-1),\n";
        gLine4="         );\n";
        gLine5="     end\n";
        gInput=gLine1+gLine2+gLine3+gLine4+gLine5;
        std::ofstream gOut("g.txt");
        gOut << gInput;
        gOut.close();
        cout<<"g.txt created"<<endl;
    }

    if (!file_exists("f_inner.txt")){
        std::string input;
        std::string line1, line2, line3, line4, line5, line6;
        line1="function f(x[2],t)\n";
        line2="         return(\n";
        line3="             x[0]*0,\n";
        line4="             x[0]*0,\n";
        line5="         );\n";
        line6="     end\n";
        input=line1+line2+line3+line4+line5+line6;
        std::ofstream out("f_inner.txt");
        out << input;
        out.close();
        cout<<"f_inner.txt created"<<endl;
    }


    //checks function syntaxis
    Function *f, *g, *f_inner;
    try{
        f = new Function ("f.txt");
    }catch(SyntaxError&){
        cout<<"Syntax error in f"<<endl;
        return false;
    }


    try{
        g = new Function ("g.txt");
    }catch(SyntaxError&){
        cout<<"Syntax error in g"<<endl;
        return false;
    }

    try{
        f_inner = new Function ("f_inner.txt");
    }catch(SyntaxError&){
        cout<<"Syntax error in f_inner"<<endl;
        return false;
    }

    return true;        //everything is OK
}



/// Using stat, checks if there is already a file with the current filename
bool file_exists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}






/// Initializes the scene in vibes for a certain plane
void init_scene(Data data, bool showFuturePos){
    vibes::beginDrawing();
    string planeName, futureName;
    int posx=220;

    if(data.plane==0){
        planeName="XY";
    }
    if (data.plane==1){
        planeName="XZ";
    }
    if (data.plane==2){
        planeName="YZ";
    }

    if (showFuturePos==true){
        futureName=" w/ predictions";
        posx+=580;                      //shows the vibes with the predictions side by side with the first one
    }
    else
    {
        futureName="";
    }
    string myTitle="Tubibex "+planeName+futureName;
    vibes::newFigure(myTitle);

    vibes::setFigureProperties(myTitle,vibesParams("x",posx,"y",220,"width",data.imageSize,"height",data.imageSize));

    vibes::axisAuto(myTitle);
    vibes::axisLabels("x","y");
    vibes::drawBox(data.initialBox[0].lb(),data.initialBox[0].ub(),data.initialBox[1].lb(),data.initialBox[1].ub());    //draws initial box
    vibes::axisLimits(-2,2,-2,2); //sets axis limits
}

void init_scene(Data data, bool showFuturePos, string myTitle){
    vibes::beginDrawing();
    int posx=120;
    if (showFuturePos==true){
        posx+=600;                      //shows the vibes with the predictions side by side with the first one
    }

    vibes::newFigure(myTitle);

    vibes::setFigureProperties(myTitle,vibesParams("x",posx,"y",220,"width",data.imageSize,"height",data.imageSize));

    vibes::axisAuto(myTitle);
    vibes::axisLabels("x","y");
    vibes::drawBox(data.initialBox[0].lb(),data.initialBox[0].ub(),data.initialBox[1].lb(),data.initialBox[1].ub());    //draws initial box
    //    vibes::axisLimits(1.0*data.initialBox[0].lb(),1.0*data.initialBox[0].ub(),1.0*data.initialBox[1].lb(),1.0*data.initialBox[1].ub()); //sets axis limits
    vibes::axisLimits(-2,2,-2,2); //sets axis limits


}


/// Draws a box in vibes in a given plane and using a given color
void draw_box(IntervalVector myInt, int color, int plane){

    //plane select
    int dim1=0;
    int dim2=1;

    if(plane==0){
        dim1=0;
        dim2=1;
    }
    if (plane==1){
        dim1=0;
        dim2=2;
    }
    if (plane==2){
        dim1=1;
        dim2=2;
    }

    //comment description is no longer accurate for the colors
    if (color==0) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub());     //current boxes data.boxes
    }
    else if (color==1) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[lightGray]black");              //future boxes from back in boxes
    }
    else if (color==2) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[darkGray]black");              //future boxes from back in boxes
    }
    else if (color==3) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[gray]black");              //future boxes from back in boxes
    }
    else if (color==4) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"b[yellow]"); //current box after contractor
    }
    else if (color==5) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"b[red]");                 //current box before contractor
    }
    else if (color==6) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[green]black");              //future boxes from unsafe boxes
    }
    else if (color==7) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[darkRed]black");              //future boxes from back in boxes
    }
    else if (color==8) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[cyan]black");              //future boxes from back in boxes
    }
    else if (color==9) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[blue]black");              //future boxes from back in boxes
    }
    else if (color==10) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[darkGray]black");              //future boxes from back in boxes
    }
    else if (color==11) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[lightGray]black");              //future boxes from back in boxes
    }
    else if (color==20) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[Green]black");              //future boxes from back in boxes
    }
    else if (color==21) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[lightGreen]black");              //future boxes from back in boxes
    }
    else if (color==22) {
        vibes::drawBox(myInt[dim1].lb(),myInt[dim1].ub(),myInt[dim2].lb(),myInt[dim2].ub(),"[lightGray]black");              //future boxes from back in boxes
    }

}


/// Using VIBes, draws the current boxes stored in data excluding the future ones. A delay can be applied to see how sivia works in real time
void draw_update(Data data, IntervalVector box, IntervalVector contractedBox){
    vibes::clearFigure("Tubibex");

    for (int i = 0; i < data.boxesOutside.size(); ++i) {
        draw_box(data.boxesOutside.at(i), 0, data.plane);
    }

    for (int i = 0; i < data.boxes.size(); ++i) {
        draw_box(data.boxes.at(i),1, data.plane);
    }

    for (int i = 0; i < data.boxesBackIn.size(); ++i) {
        draw_box(data.boxesBackIn.at(i), 4, data.plane);
    }

    for (int i = 0; i < data.boxesUnsafe.size(); ++i) {
        draw_box(data.boxesUnsafe.at(i), 0, data.plane);
    }

    if (data.calcInner){
        for (int i = 0; i < data.boxesInside.size(); ++i) {
            draw_box(data.boxesInside.at(i),8,data.plane);
        }

        for (int i = 0; i < data.boxesInsideUnsafe.size(); ++i) {
            draw_box(data.boxesInsideUnsafe.at(i),9,data.plane);
        }
        for (int i = 0; i < data.boxesInsideBackIn.size(); ++i) {
            draw_box(data.boxesInsideBackIn.at(i),6,data.plane);
        }
    }

        draw_box(box, 8, data.plane);

        draw_box(contractedBox, 8, data.plane);

    if (data.delayActivated)
    {
        //if using windows comment this
        usleep(data.delayMs*1000);

        //and uncomment this
        //Sleep(sleepMs);j
    }


}


/// Same as previous one but this one also includes the future boxes and removes the time delay
void draw_update(Data data, bool showFuturePos){

    //draw disk
    if (data.showDisk) {
        double cx=0;
        double cy=0;
        vibes::drawEllipse(cx,cy,data.diskR,data.diskR,0,"[gray]k");
        //        vibes::drawEllipse(cx,cy,1.75,0.7,0,"[gray]k");
    }

        for (int i = 0; i < data.boxes.size(); ++i) {
            draw_box(data.boxes.at(i),0, data.plane);
        }

        for (int i = 0; i < data.boxesOutside.size(); ++i) {
            draw_box(data.boxesOutside.at(i), 0, data.plane);
        }

    if (data.calcInner){
        for (int i = 0; i < data.boxesInside.size(); ++i) {
            draw_box(data.boxesInside.at(i),8,data.plane);
        }

        for (int i = 0; i < data.boxesInsideUnsafe.size(); ++i) {
            draw_box(data.boxesInsideUnsafe.at(i),9,data.plane);
        }
        for (int i = 0; i < data.boxesInsideBackIn.size(); ++i) {
            draw_box(data.boxesInsideBackIn.at(i),9,data.plane);
        }
    }


        for (int i = 0; i < data.boxesBackIn.size(); ++i) {
            draw_box(data.boxesBackIn.at(i), 4, data.plane);
        }

        for (int i = 0; i < data.boxesUnsafe.size(); ++i) {
            draw_box(data.boxesUnsafe.at(i), 5, data.plane);
        }



        vector< vector<double> >  boxVectorUnsafe;
        vector< vector<double> >  boxVectorBackIn;

        if (showFuturePos==true){
            for (int i = 0; i < data.boxesUnsafeFuture.size(); ++i) {
                if (data.drawFutureBoxesUnion){
                    vector<double>  box;
                    box.push_back(data.boxesUnsafeFuture.at(i)[0].lb());
                    box.push_back(data.boxesUnsafeFuture.at(i)[0].ub());
                    box.push_back(data.boxesUnsafeFuture.at(i)[1].lb());
                    box.push_back(data.boxesUnsafeFuture.at(i)[1].ub());
                    boxVectorUnsafe.push_back(box);
                }
                else {
                    draw_box(data.boxesUnsafeFuture.at(i), 10, data.plane);
                }
            }

            for (int i = 0; i < data.boxesBackInFuture.size(); ++i) {
                if (data.drawFutureBoxesUnion){
                    vector<double>  box;
                    box.push_back(data.boxesBackInFuture.at(i)[0].lb());
                    box.push_back(data.boxesBackInFuture.at(i)[0].ub());
                    box.push_back(data.boxesBackInFuture.at(i)[1].lb());
                    box.push_back(data.boxesBackInFuture.at(i)[1].ub());
                    boxVectorBackIn.push_back(box);
                }
                else {
                    draw_box(data.boxesBackInFuture.at(i), 10, data.plane);
                }
            }
        }

        if (data.drawFutureBoxesUnion){
            vibes::drawBoxesUnion(boxVectorUnsafe, "[darkGray]black");
            vibes::drawBoxesUnion(boxVectorBackIn, "[darkGray]black");
        }



}


/// Draws the general vector field
void draw_vector_field(Data data){

    //calculate the boundaries of the mesh
    double myXMin, myXMax, myYMin, myYMax;
    myXMin=5*data.initialBox[0].lb();
    myYMin=5*data.initialBox[1].lb();
    myXMax=5*data.initialBox[0].ub();
    myYMax=5*data.initialBox[1].ub();

    QList<IntervalVector> posVector;
    Interval dummyTime(0,0);

    Function *f = new Function("f.txt");
    int numVar=f->nb_var();

    //generates an interval Qlist with all the boxes generated from the mesh
    for (double i = myXMin; i < myXMax; i=i+data.VFepsilonX) {
        for (double j = myYMin; j < myYMax; j=j+data.VFepsilonY) {
            double pos[2][2]={
                {i-data.VFepsilonX/2.0, i+data.VFepsilonX/2.0} ,            //individual mesh box position
                {j-data.VFepsilonY/2.0, j+data.VFepsilonY/2.0}
            };
            //populate the Qlist
            IntervalVector posInt(numVar, pos);
            posInt[numVar-1]=dummyTime;
            posVector.append(posInt);

        }
    }

    draw_vector_field_boxes(posVector, 1);             //draw the generated list

}

/// Draws the vector field on the center of the provided list of boxes, can specify the color of the vectors
void draw_vector_field_boxes(QList<IntervalVector> posList, int color){
    Function *f = new Function("f.txt");
    int numVar=f->nb_var();

    for(int i = 0; i<posList.size(); i++){
        IntervalVector box = posList.at(i);
        calculate_vector(box, box[numVar-1].mid(), color);
    }

    delete f;
}


/// Same as the previous one but draws the vectors in the default color
void draw_vector_field_boxes(QList<IntervalVector> posList){
    Function *f = new Function("f.txt");
    int numVar=f->nb_var();

    for(int i = 0; i<posList.size(); i++){
        IntervalVector box = posList.at(i);
        calculate_vector(box, box[numVar-1].mid(), 1);
    }

    delete f;
}

/// Calculate the correspondient vectors
void calculate_vector(IntervalVector& box, double t, int color){

    double x1, x2, x3;
    Function *f = new Function("f.txt");
    int numVar = f->nb_var()-1;

    //origin of the vector
    x1 = box[0].mid();
    x2 = box[1].mid();

    //creation of the mesh
    double eps=0.5;  //precision of the mesh

    for(x3=box[2].lb();x3<=box[2].ub();x3+=eps){

        IntervalVector evaluationBox(numVar+1);
        for (int i=0;i<numVar;i++) evaluationBox[i]=box[i].mid();
        evaluationBox[2] = x3;
        evaluationBox[numVar]=t;                        //add the time dimension
        IntervalVector fResults = f-> eval_vector(evaluationBox);
        double f1=fResults[0].mid();
        double f2=fResults[1].mid();

        //calculate the norm and resize it
        double nn=5*sqrt(pow(f1,2)+pow(f2,2));

        //arrow coordinates
        vector<double> arrowX, arrowY;
        arrowX.push_back(x1);
        arrowX.push_back(x1+eps*f1/nn);

        arrowY.push_back(x2);
        arrowY.push_back(x2+eps*f2/nn);

        draw_arrow(arrowX, arrowY, color);

    }
}


/// Draw arrow using a line and a box from the given coordinates
void draw_arrow(vector<double> arrowX, vector<double> arrowY, int color){
    if (color==0) {
        vibes::drawLine(arrowX, arrowY, "[cyan]");
        vibes::drawBox(arrowX[1]-0.005,arrowX[1]+0.005,arrowY[1]-0.005,arrowY[1]+0.005,"[black]");
    }
    if (color==1) {
        vibes::drawLine(arrowX, arrowY, "[cyan]");
        vibes::drawBox(arrowX[1]-0.005,arrowX[1]+0.005,arrowY[1]-0.005,arrowY[1]+0.005,"[black]");
    }
}

void export_all_images(Data& data){

    string myTitle, myFilename;

    //    myTitle="Tube_noInner_noOuter_noVField";
    //    init_scene(data, 0, myTitle);
    //    draw_update(data, 0);
    //    //draw_vector_field(data);
    //    myFilename=myTitle+".png";
    //    vibes::saveImage(myFilename);
    //    vibes::endDrawing();

    //    myTitle="Tube_noInner_noOuter_VField";
    //    init_scene(data, 0, myTitle);
    //    draw_update(data, 0);
    //    draw_vector_field(data);
    //    myFilename=myTitle+".png";
    //    vibes::saveImage(myFilename);
    //    vibes::endDrawing();

    //    myTitle="Tube_noInner_Outer_noVField";
    //    init_scene(data, 1, myTitle);
    //    draw_update(data, 1);
    //    //draw_vector_field(data);
    //    myFilename=myTitle+".png";
    //    vibes::saveImage(myFilename);
    //    vibes::endDrawing();

    data.calcInner=false;
    myTitle="Tube_noInner_Outer_VField";
    init_scene(data, 1, myTitle);
    draw_update(data, 1);
    draw_vector_field(data);
    myFilename=myTitle+".png";
    vibes::saveImage(myFilename);
    vibes::endDrawing();

    //    data.showDisk=true;
    //    myTitle="Tube_Inner_noOuter_noVField";
    //    init_scene(data, 0, myTitle);
    //    draw_update(data, 0);
    //    //draw_vector_field(data);
    //    myFilename=myTitle+".png";
    //    vibes::saveImage(myFilename);
    //    vibes::endDrawing();

    data.calcInner=true;
    myTitle="Tube_Inner_noOuter_VField";
    init_scene(data, 0, myTitle);
    draw_update(data, 0);
    draw_vector_field(data);
    myFilename=myTitle+".png";
    vibes::saveImage(myFilename);
    vibes::endDrawing();

    //    myTitle="Tube_Inner_Outer_noVField";
    //    init_scene(data, 1, myTitle);
    //    draw_update(data, 1);
    //    //draw_vector_field(data);
    //    myFilename=myTitle+".png";
    //    vibes::saveImage(myFilename);
    //    vibes::endDrawing();

    myTitle="Tube_Inner_Outer_VField";
    init_scene(data, 1, myTitle);
    draw_update(data, 1);
    draw_vector_field(data);
    myFilename=myTitle+".png";
    vibes::saveImage(myFilename);
    vibes::endDrawing();
}
