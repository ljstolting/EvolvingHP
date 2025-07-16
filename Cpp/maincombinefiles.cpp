//if I artificially parallelized code by making it output a slice
//in several pieces, load those pieces together into one file
//without having to use copy/paste

#include "random.h"
#include "VectorMatrix.h"

const int num_rows1 = 2905050; //first four files same length
const int num_rows2 = 2963151; //last longer because has to include end row
const int num_cols = 3;

int main(){
    ifstream slice1;
    ifstream slice2;
    ifstream slice3;
    ifstream slice4;
    ifstream slice5;

    slice1.open("./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15slice1.dat");
    slice2.open("./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15slice2.dat");
    slice3.open("./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15slice3.dat");
    slice4.open("./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15slice4.dat");
    slice5.open("./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15slice5.dat");

    ofstream slice;
    slice.open("./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15.dat");

    TVector<double> row(1,num_cols);
    for (int i=0;i<num_rows1;i++){
        slice1 >> row;
        // for (int j=1;j<=num_cols;j++){
        //     row(j) = row(j) * 3; //the data was messed up and divided three times as many steps as it took, so rescale
        // }
        slice << row << endl;
    }
    for (int i=0;i<num_rows1;i++){
        slice2 >> row;
        // for (int j=1;j<=num_cols;j++){
        //     row(j) = row(j) * 3; 
        // }
        slice << row << endl;
    }
    for (int i=0;i<num_rows1;i++){
        slice3 >> row;
        // for (int j=1;j<=num_cols;j++){
        //     row(j) = row(j) * 3; 
        // }
        slice << row << endl;
    }
    for (int i=0;i<num_rows1;i++){
        slice4 >> row;
        // for (int j=1;j<=num_cols;j++){
        //     row(j) = row(j) * 3; 
        // }
        slice << row << endl;
    }
    for (int i=0;i<num_rows2;i++){
        slice5 >> row;
        // for (int j=1;j<=num_cols;j++){
        //     row(j) = row(j) * 3; //the data was messed up and divided by all three cycles instead of just the first one
        // }
        slice << row << endl;
    }

    slice1.close();
    slice2.close();
    slice3.close();
    slice4.close();
    slice5.close();
    slice.close();

    return 0;
}