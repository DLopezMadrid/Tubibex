#include "sivia.h"
#include <QList>
#include "vibes.h"
#include <QDateTime>

#define __PREC__ 1e-11
#define __METH__ RK4
#define __DURATION__ 0.5

/// Initializes constraints and prepares data to be processed
Sivia::Sivia(Data &data, bool calcInner)
{
    //constraints to check to see if a trajectory belongs to a tube

    //c1=   gdot= d/dx(gi)(x,t)*f(x,t)+d/dt(gi)(x,t)>=0
    //c2=   gi(x,t)=0j
    //c3=   g(x,t)<=0

    Function g("g.txt");
    Function dg(g, Function::DIFF);                                     //  d/dx(gi)(x,t)

    Variable x(data.numVarF),t;                                         //we have x[] and t as variables for our fns

    // initialize auxMat and auxVector to the correct sizes and fill with zeros
    IntervalMatrix auxMat(data.numVarF+1, data.numVarF,Interval::ZERO);
    IntervalVector auxVector(data.numVarF+1,Interval::ZERO);

    //put 1 in the diagonal of auxMat
    for (int i=0; i<data.numVarF; i++){
        auxMat[i][i]=1;}

    auxVector[data.numVarF]=1;

    if (calcInner){         //for the inner approximation of the tube, we set a new function f and its correspondent constraints
        cout<<endl<<"Start inner approx calculation"<<endl;
        Function f=("f.txt");
        Function gdot(x,t,dg(x,t)*(auxMat*transpose(f(x,t))+auxVector));


        NumConstraint c3(g, EQ);                                           //LEQ means less or equal 0

        Array<Ctc> individualTubeConstraints;                               //put together the constraints in an array

        individualTubeConstraints.add(*new CtcHC4(Array<NumConstraint>(c3)));

        CtcUnion unionTubeConstraints(individualTubeConstraints);           //calculate the union

        Ctc3BCid tubeConstraints(unionTubeConstraints);                     //this contracts the whole union to its final form

        data.boxes.push_back(data.initialBox);                              //initialize the boxes

        do_Sivia(tubeConstraints, data, gdot, calcInner);

        print_results(data);
    }

    else{                       //for the outer approximation
        Function f("f.txt");
        Function gdot(x,t,dg(x,t)*(auxMat*transpose(f(x,t))+auxVector));
        //c1 & c2
        Array<NumConstraint> c1, c2;

        int numConstraints = data.g->expr().dim.max_index()+1;              //find how many gi we have

        for (int i = 0; i < numConstraints; ++i) {                          //create constraints based on the dimensions of g
            c1.add(*new NumConstraint(x,t,gdot(x,t)[i] >= 0));
            c2.add(*new NumConstraint(x,t,g(x,t)[i] = 0));
        }

        NumConstraint c3(g, LEQ);                                           //LEQ means less or equal 0

        Array<Ctc> individualTubeConstraints;                               //put together the constraints in an array

        for (int i=0;i<numConstraints ;i++) {
            individualTubeConstraints.add(*new CtcHC4(Array<NumConstraint>(c1[i],c2[i],c3)));}

        CtcUnion unionTubeConstraints(individualTubeConstraints);           //calculate the union

        Ctc3BCid tubeConstraints(unionTubeConstraints);                     //this contracts the whole union to its final form

        data.boxes.push_back(data.initialBox);                              //initialize the boxes

        do_Sivia(tubeConstraints, data, gdot, calcInner);
        if (!data.calcInner){
            print_results(data);
        }
    }
}



/// Processes the data using contractors and bissections. Classifies the boxes in outside (grey), back_in(yellow) and unsafe (red)
void Sivia::do_Sivia(Ctc& tubeConstraints, Data &data, Function gdot, bool calcInner){

    QTime tSivia;
    tSivia.start();

    if (calcInner)                  //inner approximation calculation
    {
        int count=0;
        while (!data.boxes.empty()) {
            IntervalVector currentBox = data.boxes.front();                 //start from the first one
            data.boxes.pop_front();                                         //once it has been copied remove the first box

            IntervalVector auxBox=currentBox;                               //store it in aux variable to compare later

            tubeConstraints.contract(currentBox);                           //contract the current box using the previously calculated constraints
            if (currentBox!=auxBox){                                        //if the box has been contracted
                IntervalVector* removedByContractorInner;
                int setDiff=auxBox.diff(currentBox, removedByContractorInner);   //set difference between the contracted box and the original box
                for (int i = 0; i < setDiff; ++i) {

                    bool testInside=true;
                    IntervalVector gg=data.g->eval_vector(removedByContractorInner[i]);

                    for(int j = 0; j<gg.size(); j++){
                        testInside = testInside && (gg[j].ub()<=0);
                    }
                    if (testInside) {
                        data.boxesInside.append(removedByContractorInner[i]);
                    }
                }
                delete[] removedByContractorInner;
            }

            if(data.realTimeDraw){                                          //draw the boxes processing in real time
                draw_update(data, auxBox, currentBox);
            }


            bool allBoxesLessEpsilon=true;                                                                              //check if all the boxes are smaler than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox.diam()[i])<=data.epsilons[i]));
            }
            allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox[currentBox.size()-1].diam())<=data.dt));         //check the time box also


            bool boxesLessEpsilon=false;                                                                                //check if at least one box is smaller than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                boxesLessEpsilon = boxesLessEpsilon||((currentBox[i].diam())<=data.epsilons[i]);
            }
            boxesLessEpsilon = boxesLessEpsilon&&((currentBox[currentBox.size()-1].diam())<=data.dt);                   //check time box

            if (allBoxesLessEpsilon) {                                                          //if allBoxesLessEpsilon = true the box is unsafe and I continue my loop
                (data.boxesInsideUnsafe).push_back(currentBox);
                count++;
                if (count >=data.maxNumUnsafeBoxes && data.maxNumUnsafeBoxesActivated){         //If I have more boxes than nbPerhaps I stop the loop and I display the results
                    break;
                }
            }
            else {                                                                              //Otherwise we bissect following the widest diameter
                double l = 0;
                double l_temp = 0;
                int v = -1;
                for(int i = 0; i<currentBox.size()-1; i++){                                     //test that the diameter of the boxes doesnt depend on time
                    if(currentBox[i].is_bisectable()||!(currentBox[i].is_degenerated())){
                        l_temp = currentBox[i].diam();
                        if(l_temp>=data.epsilons[i] && l_temp/(data.epsilons[i]) > l){
                            l = l_temp/(data.epsilons[i]);
                            v = i;
                        }
                    }
                }

                l_temp = currentBox[currentBox.size()-1].diam();                                //test the time interval
                if(l_temp>=data.dt && l_temp/(data.dt) > l){
                    v = currentBox.size()-1;
                }
                if(v != -1 && currentBox[v].is_bisectable()){                                   // then the test interval of the state variables, and then it bisects the interval which has the largest diameter
                    pair<IntervalVector,IntervalVector> boxes=currentBox.bisect(v, 0.5);
                    (data.boxes).push_back(boxes.first);
                    (data.boxes).push_back(boxes.second);
                }
                else{
                    if (data.myDebug){
                        std::cout<<"Cannot be bisected \n";
                    }
                }
            }
        }
    }


    else                            //outer approximation
    {
        int count=0;
        //process all the boxes in data
        while (!data.boxes.empty()) {
            IntervalVector currentBox = data.boxes.front();                 //start from the first one
            data.boxes.pop_front();                                         //once it has been copied remove the first box

            IntervalVector auxBox=currentBox;                               //store it in aux variable to compare later

            tubeConstraints.contract(currentBox);                           //contract the current box using the previously calculated constraints
            if (currentBox!=auxBox){                                        //if the box has been contracted
                IntervalVector* removedByContractor;
                int setDiff=auxBox.diff(currentBox, removedByContractor);   //set difference between the contracted box and the original box
                for (int i = 0; i < setDiff; ++i) {
                    data.boxesOutside.push_back(removedByContractor[i]);    //add the areas removed by the contractor to the outside set
                }
                delete[] removedByContractor;
            }

            if(data.realTimeDraw){                                          //draw the boxes processing in real time
                draw_update(data, auxBox, currentBox);
            }


            bool allBoxesLessEpsilon=true;                                                                              //check if all the boxes are smaler than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox.diam()[i])<=data.epsilons[i]));
            }
            allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox[currentBox.size()-1].diam())<=data.dt));         //check the time box also

            bool boxesLessEpsilon=false;                                                                                //check if at least one box is smaller than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                boxesLessEpsilon = boxesLessEpsilon||((currentBox[i].diam())<=data.epsilons[i]);
            }
            boxesLessEpsilon = boxesLessEpsilon&&((currentBox[currentBox.size()-1].diam())<=data.dt);                   //check time box

            if (boxesLessEpsilon && !allBoxesLessEpsilon){
                IntervalVector xnext = currentBox.subvector(0, data.numVarF-1).mid();   //using the middle point of the box calculate the future positions using euler method
                IntervalVector x = currentBox.mid();
                bool testBackIn;
                for (int i = 0;i<data.numFuturePos;i++){                                // Euler method: x(n+1)=x(n)+dt*fx
                    x[data.numVarF]= x[data.numVarF].mid();
                    testBackIn = true;
                    xnext=xnext+(data.dt)*data.f->eval_vector(x);
                    x.put(0, xnext);
                    x[data.numVarF] = x[data.numVarF]+(data.dt);
                    IntervalVector gg=data.g->eval_vector(x);
                    for(int j = 0; j<gg.size(); j++){
                        testBackIn = testBackIn && (gg[j].ub()<0);                      //test if it comes back to the bubble

                    }
                    if(testBackIn == true){
                        break;
                    }
                }

                if(testBackIn == true && data.enableBackIn){                                                 //If my box was back in the bubble after integration, I store it in boxesbackin
                    (data.boxesBackIn).append(currentBox);

                    continue;
                }
            }


            if (allBoxesLessEpsilon) {                                                          //if allBoxesLessEpsilon = true the box is unsafe and I continue my loop
                (data.boxesUnsafe).push_back(currentBox);
                count++;
                if (count >=data.maxNumUnsafeBoxes && data.maxNumUnsafeBoxesActivated){         //If I have more boxes than nbPerhaps I stop the loop and I display the results
                    break;
                }
            }
            else {                                                                              //Otherwise we bissect following the widest diameter
                double l = 0;
                double l_temp = 0;
                int v = -1;
                for(int i = 0; i<currentBox.size()-1; i++){                                     //test that the diameter of the boxes doesnt depend on time
                    if(currentBox[i].is_bisectable()||!(currentBox[i].is_degenerated())){
                        l_temp = currentBox[i].diam();
                        if(l_temp>=data.epsilons[i] && l_temp/(data.epsilons[i]) > l){
                            l = l_temp/(data.epsilons[i]);
                            v = i;
                        }
                    }
                }

                l_temp = currentBox[currentBox.size()-1].diam();                                //test the time interval
                if(l_temp>=data.dt && l_temp/(data.dt) > l){
                    v = currentBox.size()-1;
                }
                if(v != -1 && currentBox[v].is_bisectable()){                                   // then the test interval of the state variables, and then it bisects the interval which has the largest diameter
                    pair<IntervalVector,IntervalVector> boxes=currentBox.bisect(v, 0.5);
                    (data.boxes).push_back(boxes.first);
                    (data.boxes).push_back(boxes.second);
                }
                else{
                    if (data.myDebug){
                        std::cout<<"Can not be bisected \n";
                    }
                }
            }
        }


        double maxGValues[data.numVarG-1];                              //init vector to store the max values of G

        for (int i = 0; i < data.numVarG-1; ++i) {
            maxGValues[i]=0; }


        for(int i=0; i<data.boxesUnsafe.size();i++) {                   //process unsafe boxes

            IntervalVector currentBox=data.boxesUnsafe.at(i);
            IntervalVector nextBox = currentBox.subvector(0, data.numVarF-1);

            if (data.intMethod==0){                                     //Guaranteed integration

                // State variables
                Variable y(data.numVarF);

                // Initial conditions
                IntervalVector yinit(data.numVarF);

                for (int i = 0; i < data.numVarF; ++i) {
                    yinit[i] = currentBox[i];
                    cout<<currentBox[i]<<endl;
                }

                // system fn has to be re entered here, cannot be loaded directly from text file

                //pendulum
                Function ydot = Function (y,Return (y[1], -sin(y[0])-0.15*y[1]));

                //non holonomic
                //                Interval t = currentBox[data.numVarF];
                //                Interval xd = 7*t;
                //                Interval xdd = 7;
                //                Interval yd = sin(0.1*t);
                //                Interval ydd = 0.1*cos(0.1*t);
                //                Interval xdiff = (xd-y[0]+xdd);
                //                Interval ydiff = (yd-y[1]+ydd);
                //                Interval norm =  ( sqrt((xdiff)^2 +(ydiff)^2) );

                //                Function ydot = Function (y,Return (( sqrt((xd-y[0]+xdd)*(xd-y[0]+xdd) +((yd-y[1]+ydd))*(yd-y[1]+ydd)) )*cos(y[2]), ( sqrt(((xd-y[0]+xdd))*(xd-y[0]+xdd) +((yd-y[1]+ydd))*(yd-y[1]+ydd)) )*sin(y[2]), 10*(cos(y[2])*((yd-y[1]+ydd))-sin(y[2])*((xd-y[0]+xdd)))/( sqrt(((xd-y[0]+xdd))*(xd-y[0]+xdd) +((yd-y[1]+ydd))*(yd-y[1]+ydd)) )));           // Ivp contruction (initial time is 0.0)

                QTime t1;
                t1.start();

                ivp_ode problem = ivp_ode (ydot, 0.0 , yinit);

                // Simulation construction and run
                simulation simu = simulation (&problem,data.dt*data.numFuturePos, __METH__, __PREC__);          //uses Runge-kutta4 method
                data.boxesUnsafeFuture.append(simu.run_simulation());                                           //modified ibex_simulation.h to make it return a list with all the solutions, not just the last one

                double timeSiviaCalculations1=t1.elapsed()/1000.0;
                double timeSiviaCalculationsTotal=tSivia.elapsed()/1000.0;
                cout<<endl<<"Unsafe # "<<i<<"  , Box time = "<<timeSiviaCalculations1<<" , Total elapsed time = "<<timeSiviaCalculationsTotal<<endl;
            }

            if (data.intMethod==1){                                     //euler method

                for (int i = 0;i<data.numFuturePos;i++){

                    IntervalVector gdotValue=gdot.eval_vector(currentBox);                                      //evaluate the g and gdot functions to inspect the constraints
                    IntervalVector gValue=data.g->eval_vector(currentBox);

                    if (data.myDebug){
                        cout<<"box = "<<currentBox<<endl;

                        for(int j = 0; j<gValue.size(); j++){
                            cout<<"gdot"<<j<<" = "<<gdotValue<<"  /  "<<(gdotValue[j].lb()>0) <<endl;                       //gdot i values
                        }
                    }

                    for(int j = 0; j<gValue.size(); j++){
                        if (data.myDebug){
                            cout<<"g"<<j<<" = "<<gValue[j]<<"  /  "<<((gValue[j].ub()>0)&&(gValue[j].lb()<0)) <<endl;}      //print gi values

                        if((gValue[j].ub()>maxGValues[j])&&(gValue[j].ub()<999999)){                                        //check max values for each gi, ignore if system goes to infinity
                            maxGValues[j]=gValue[j].ub();}
                    }

                    nextBox=nextBox+(data.dt)*data.f->eval_vector(currentBox);              //euler method
                    data.boxesUnsafeFuture.append(nextBox);

                    currentBox.put(0, nextBox);
                    currentBox[data.numVarF] = currentBox[data.numVarF]+(data.dt);          //increase time for the next step
                }
            }
        }


        for(int i=0; i<data.boxesBackIn.size();i++){                                        //process back_in boxes

            IntervalVector currentBox=data.boxesBackIn.at(i);
            IntervalVector nextBox = currentBox.subvector(0, data.numVarF-1);

            if (data.intMethod==0){                 //guaranteed integration                                                       //Guaranteed integration

                // State variables
                Variable y(data.numVarF);

                // Initial conditions
                IntervalVector yinit(data.numVarF);

                for (int i = 0; i < data.numVarF; ++i) {
                    yinit[i] = currentBox[i];
                    cout<<currentBox[i]<<endl;
                }

                QTime t2;
                t2.start();

                // system fn has to be re entered here, cannot be loaded directly from text file

                //pendulum
                Function ydot = Function (y,Return (y[1], -sin(y[0])-0.15*y[1]));

                //non holonomic
                //                Interval t = currentBox[data.numVarF];
                //                Interval xd = 7*t;
                //                Interval xdd = 7;
                //                Interval yd = sin(0.1*t);
                //                Interval ydd = 0.1*cos(0.1*t);
                //                Interval xdiff = (xd-y[0]+xdd);
                //                Interval ydiff = (yd-y[1]+ydd);
                //                Interval norm =  ( sqrt((xdiff)^2 +(ydiff)^2) );

                //                Function ydot = Function (y,Return (( sqrt((xd-y[0]+xdd)*(xd-y[0]+xdd) +((yd-y[1]+ydd))*(yd-y[1]+ydd)) )*cos(y[2]), ( sqrt(((xd-y[0]+xdd))*(xd-y[0]+xdd) +((yd-y[1]+ydd))*(yd-y[1]+ydd)) )*sin(y[2]), 10*(cos(y[2])*((yd-y[1]+ydd))-sin(y[2])*((xd-y[0]+xdd)))/( sqrt(((xd-y[0]+xdd))*(xd-y[0]+xdd) +((yd-y[1]+ydd))*(yd-y[1]+ydd)) )));           // Ivp contruction (initial time is 0.0)


                ivp_ode problem = ivp_ode (ydot, currentBox[data.numVarF].lb() , yinit);

                // Simulation construction and run
                simulation simu = simulation (&problem,data.dt*data.numFuturePos, __METH__, __PREC__);          //uses Runge-kutta4 method
                data.boxesUnsafeFuture.append(simu.run_simulation());                       //modified ibex_simulation.h to make it return a list with all the solutions, not just the last one

                double timeSiviaCalculations2=t2.elapsed()/1000.0;
                double timeSiviaCalculationsTotal=tSivia.elapsed()/1000.0;

                cout<<endl<<"Back_in # "<<i<<"  , Box time = "<<timeSiviaCalculations2<<" , Total elapsed time = "<<timeSiviaCalculationsTotal<<endl;
            }

            for (int i = 0;i<data.numFuturePos;i++){                                        //euler method

                if (data.intMethod==1){
                    IntervalVector gdotValue=gdot.eval_vector(currentBox);                  //evaluate the g and gdot functions to inspect the constraints
                    IntervalVector gValue=data.g->eval_vector(currentBox);

                    if (data.myDebug){
                        cout<<"box = "<<currentBox<<endl;

                        for(int j = 0; j<gValue.size(); j++){
                            cout<<"gdot"<<j<<" = "<<gdotValue<<"  /  "<<(gdotValue[j].lb()>0) <<endl;               //gdoti values
                        }
                    }

                    bool testBackIn= true;
                    for(int j = 0; j<gValue.size(); j++){
                        if (data.myDebug){
                            cout<<"g"<<j<<" = "<<gValue[j]<<"  /  "<<((gValue[j].ub()>0)&&(gValue[j].lb()<0)) <<endl;}          //print gi values

                        testBackIn = testBackIn && (gValue[j].ub()<0);                                              //test if it comes back to the bubble

                        if((gValue[j].ub()>maxGValues[j])&&(gValue[j].ub()<999999)){                                //check max values for each gi, ignore if system goes to infinity
                            maxGValues[j]=gValue[j].ub();
                        }
                    }


                    nextBox=nextBox+(data.dt)*data.f->eval_vector(currentBox);                                      // euler method
                    if (!testBackIn){
                        data.boxesBackInFuture.append(nextBox);
                    }
                    currentBox.put(0, nextBox);
                    currentBox[data.numVarF] = currentBox[data.numVarF]+(data.dt);                                  //increase time for the next step
                }
            }
        }


        if (data.myDebug){
            for (int i = 0; i < data.numVarG-1; ++i) {
                cout<<"Max G"<<i<<" = "<<maxGValues[i]<<endl;}
        }
    }
}

/// Uses Vibes to represent the results
void Sivia::print_results(Data data){
    cout<<"Boxes outer approximation = "<<data.boxesBackIn.size()+data.boxesOutside.size()+data.boxesUnsafe.size()<<endl;
    cout<<"Outside boxes = "<<data.boxesOutside.size()<<endl;
    cout<<"Unsafe boxes = "<<data.boxesUnsafe.size()<<endl;
    cout<<"Back in boxes = "<<data.boxesBackIn.size()<<endl;
    cout<<"Boxes inner approximation = "<<data.boxesInside.size()+data.boxesInsideUnsafe.size()+data.boxesInsideBackIn.size()<<endl;
    cout<<"Inside boxes = "<<data.boxesInside.size()<<endl;
    cout<<"Inside unsafe boxes = "<<data.boxesInsideUnsafe.size()<<endl;

    //use vibes to show the results
    //    bool showFuturePos=false;
    //    init_scene(data, showFuturePos);
    //    draw_update(data, showFuturePos);
    //    draw_vector_field(data);
    //    vibes::endDrawing();

    //    showFuturePos=true;
    //    init_scene(data, showFuturePos);
    //    draw_update(data, showFuturePos);
    //    draw_vector_field(data);
    //    vibes::endDrawing();

    if (data.saveAllPossibleFigures){
        export_all_images(data);}

}

