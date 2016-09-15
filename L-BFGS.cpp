//
//
// TODO: dimension d
// TODO: add angle transformation




#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include "meta.h"
#include "problem.h"
#include "lbfgssolver.h"

#define PI 3.1415926

using namespace std;


void openFile (ifstream & inputfile, string name){
    inputfile.open(name.c_str());
    if (inputfile.fail()) {
        cout << "Error opening input data file\n";
        exit(1);
    }
}

string ParseControlFile(ifstream & inputfile, int & dim, int &numpts, double & s, int & c, int & max_neighbor, int & numFile, int & numIteration, bool & infile){
    
    
    infile = false;
    string filename = "";
    
    string line;
    int lineNumber = 0;
    while (! inputfile.eof()) {
        lineNumber++;
        getline(inputfile, line);
        stringstream tmp(line);
        string k;
        if (lineNumber == 4) {
            tmp >> k >> k >> s;
        }else if (lineNumber == 5) {
            tmp >> k >> k >> dim;
        }else if (lineNumber == 6) {
            tmp >> k >> k >> c;
        }else if (lineNumber == 7) {
            tmp >> k >> k >> infile;
        }else if (lineNumber == 8) {
            tmp >> k >> k >> numpts;
        }else if (lineNumber == 9) {
            tmp >> k >> k >> numIteration;
        }else if (lineNumber == 10) {
            tmp >> k >> k >> numFile;
        }else if (lineNumber == 11) {
            tmp >> k >> k >> max_neighbor;
        }else if (lineNumber == 12) {
            tmp >> k >> k >> filename;
        }
    }
    cout << "\nSummary of the control file:\n\n";
    cout << "S value: " << s << "\n";
    cout << "Dimension: " << dim << "\n";
    cout << "C value: " << c << "\n";
    cout << "Infile request: " << infile << "\n";
    cout << "Number of points: " << numpts << "\n";
    cout << "Number of iterations: " << numIteration << "\n";
    cout << "Number of output files: " << numFile << "\n";
    cout << "Max neighbor: " << max_neighbor << "\n";
    if (infile) cout << "Input filename: " << filename << "\n\n";
    else cout << "No input file request; program will generate a random configuartion.\n\n";
    return filename;
}


void ToVector(const cppoptlib::Matrix<double> & M, cppoptlib::Vector<double> & V )
{
    int c = M.cols();
    for (int i=0; i<M.rows(); ++i )
        V.segment(i*c,c) = M.row(i);
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////



// NOTE!!!!
// the first n-2 angles range from 0 to PI, while the last angle range from 0 to 2 PI,
// in 3d the coordinate is permuted in the order of z, x, y

void ToND(const cppoptlib::Vector<double> & angles, cppoptlib::Vector<double> & coords) {
    double tmp = 1;
    coords[0] = cos(angles[0]);
    for (int i = 1; i < coords.size(); i++) {
	   	tmp *= sin(angles[i-1]);
        coords[i] = cos(angles[i])*tmp;
//        cout << coords[i] << endl;
    }
    coords[coords.size()-1] = tmp*sin(angles[angles.size()-1]);
}
















double dist_squared(const cppoptlib::Vector<double> & angles1,const cppoptlib::Vector<double> & angles2)
{
    return (2-2*(sin(angles1(1))*sin(angles2(1))*cos(angles1(0)-angles2(0))+cos(angles1(1))*cos(angles2(1))));
}



void To3D(const cppoptlib::Vector<double> & angles, cppoptlib::Vector<double> & coords)
{
    coords(0) = cos(angles(0)) * sin(angles(1));
    coords(1) = sin(angles(0)) * sin(angles(1));
    coords(2) = cos(angles(1));
}



void ComputeJacobian(const double & theta, const double & phi, cppoptlib::Matrix<double> & temp){
    //x = sin(phi) cos(theta)
    //y = sin(phi) sin(theta)
    //z = cos(phi)
    //
    // derivatives w.r.t. theta:
    temp(0,0) = -sin(phi) * sin(theta); // x
    temp(1,0) =  sin(phi) * cos(theta); // y
    temp(2,0) =  0;                      // z
    // derivatives w.r.t. phi:
    temp(0,1) =  cos(phi) * cos(theta); // x
    temp(1,1) =  cos(phi) * sin(theta); // y
    temp(2,1) = -sin(phi);             // z
}



void ToAngles(cppoptlib::Matrix<double> & all_points, cppoptlib::Matrix<double> & all_angles)
{
    //cppoptlib::Matrix<double> angles(numpts, dim-1);
    for (int i = 0; i < all_points.rows(); i++) {
        double x = all_points(i,0);
        double y = all_points(i,1);
        double z = all_points(i,2);
        double r = sqrt(x*x+y*y+z*z);
        all_angles(i,0) = atan2(y,x);
        all_angles(i,1) = acos(z/r);
    }
}



// do cutoff need to change?

double cutoff(double distance, double cutoff_radius){
    return pow(1-pow(distance/cutoff_radius, 4), 3);
}

double cutoffPrime(double distance, double cutoff_radius){
    double t = distance / cutoff_radius;
    return 12*pow(1-pow(t,4),2)*pow(t,3)/distance/cutoff_radius;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////





double Energy(const cppoptlib::Vector<double> & V, const double s_power, const int dim)
{
    // V contains spherical coordinates
    double e = 0;
    for (int i=0; i<V.size()/dim; ++i)
    {
        for (int j=0; j<i; ++j)
        {
            e += pow(dist_squared(V.segment(i*dim, dim), V.segment(j*dim, dim)), -s_power/2.0);
        }
    }
    return 2.0 * e;
}




void writeFile (ofstream & outputfile, string name, cppoptlib::Vector<double> V, int dim){
    outputfile.open(name.c_str());
    if (outputfile.fail()) {
        cout << "Error writing output data file" << endl;
//        exit(1);
    }
    int c = dim-1;
    
    outputfile << setprecision(6);
    outputfile << fixed;
    
    cppoptlib::Vector<double> tmp2(dim);
    cppoptlib::Vector<double> tmp(c);
    
    for (int i =0; i < V.rows()/c; i++) {
        tmp = V.segment(i*c, c);
        To3D(tmp, tmp2);
        
        
        outputfile << tmp2(0);
        for (int j = 1; j < dim; j++) {
            outputfile << "\t" << tmp2(j);
        }
        outputfile << "\n";
    }
}







class minimizeEnergy : public cppoptlib::Problem<double> {
    
    double cutoff_radius;
    double s;
    int cubes_per_side;
    int dim;
    int numpts;
    cppoptlib::Matrix<double> Cubes;
    cppoptlib::Matrix<double> ptsND;
    
    
public:
    minimizeEnergy(double r, double s_value, int d, int n, int c, int c_cube, int max_neighbor)
    :cutoff_radius(r), s(s_value), cubes_per_side(c), dim(d), numpts(n), Cubes(max_neighbor+1, c_cube), ptsND(dim, numpts){};
    
    
    double TruncatedEnergy(const cppoptlib::Vector<double> & V, const cppoptlib::Vector<int> & neighbors,
                           const int index, int d);
    
    int Point2Cube(const cppoptlib::Vector<double> & thisPt);
    
    
    void FindNeighborCubes(int cube_index, cppoptlib::Vector<int> & neighbors);
    
    
    void BuildIndex(const cppoptlib::Vector<double> &);
    
    
    double value(const cppoptlib::Vector<double> &x) 
    {
        double total_energy = 0;
        
        //int max_neighbor = Cubes.rows();
        cppoptlib::Vector<int> neighbor_cube_indices (pow(3, dim));
        
        // get the 3d points matrix -- pts
        // assign points to cubes
        // the first field in each column is the number of points inside that column.
        BuildIndex(x);
        
        // find the neighbor_cube_indices for each cube and calculate the energy.
        for (int i = 0; i < Cubes.cols(); ++i) 
        {
            neighbor_cube_indices = cppoptlib::Vector<int>::Ones(pow(3, dim))*(-1);
            FindNeighborCubes(i, neighbor_cube_indices);
            //cout << neighbor_cube_indices.transpose() << endl << endl;
            int points_in_cube_number = Cubes(0, i);
            // j goes over all points in the i-th cube
            for (int j = 1; j <= points_in_cube_number; ++j) {
                int point_index = Cubes(j, i); // absolute point index
                total_energy += TruncatedEnergy(x, neighbor_cube_indices, point_index, dim);
            }
        }
//        cout << fixed;
//        cout << "Energy: " << total_energy << endl;
        
        return total_energy;
    }
    
    

    void gradient(const cppoptlib::Vector<double> &x, cppoptlib::Vector<double> &grad) {
        int dim_angle = dim-1;
        int max_neighbor = Cubes.rows();
        double distance;
        cppoptlib::Vector<int> neighbor_cube_indices (pow(3, dim));
        cppoptlib::Vector<double> temp_sum(dim), temp(dim);
        cppoptlib::Matrix<double> temp_jacobian(dim, dim_angle);
        temp_sum.setZero();
        BuildIndex(x);
        
        // find the neighbor_cube_indices for each cubes and calculate the energy.
        for (int index_cube = 0; index_cube < Cubes.cols(); ++index_cube)  
        {
            neighbor_cube_indices = cppoptlib::Vector<int>::Ones(pow(3, dim))*(-1);
            FindNeighborCubes(index_cube, neighbor_cube_indices);

            int points_in_cube = Cubes(0, index_cube);
            for (int j = 1; j <= points_in_cube; ++j) 
            {
                temp_sum.setZero();
                int point_index = Cubes(j, index_cube);
                for (int k = 0; k < neighbor_cube_indices.size(); ++k) 
                {
                    int tmp = neighbor_cube_indices(k);
                    if (tmp != -1) 
                    {
                        int points_in_other_cube = Cubes(0, tmp);
                        for (int l = 1; l <= points_in_other_cube; l++) 
                        {
                            int other_point_index = Cubes(l, tmp);
                            temp = ptsND.col(point_index) - ptsND.col(other_point_index);
                            distance = sqrt(temp.dot(temp));
                            if (other_point_index != point_index && distance < cutoff_radius)
                            {
                                temp_sum += (cutoffPrime(distance, cutoff_radius)*pow(distance, -s) +
                                (-s)*cutoff(distance,cutoff_radius)*pow(distance, -s-2))*temp;
                            }
                        }
                    }
                }
                //
                
                ComputeJacobian(x(point_index*dim_angle+0), x(point_index*dim_angle+1), temp_jacobian);
                grad.segment(point_index*dim_angle, dim_angle) = temp_sum.transpose() * temp_jacobian;
            }
        }

    }
};


double minimizeEnergy::TruncatedEnergy(const cppoptlib::Vector<double> & V,
                       const cppoptlib::Vector<int> & neighbors, const int index, int dim){
    double energy = 0;
    int c = dim -1;
    for (int i = 0; i < neighbors.size(); ++i) 
    {
        int tmp = neighbors(i); // goes over all neighbor cubes
        if (tmp != -1) // check if there is a neighbor in this direction
        { 
            int points_in_cube_number = Cubes(0, tmp);
            for (int j = 1; j <= points_in_cube_number; j++) {
                int point_index = Cubes(j, tmp);
                double distance = sqrt(dist_squared(V.segment(index*c, c), V.segment(point_index*(c), c)));
                if (point_index != index && distance < cutoff_radius)
                    energy += pow(distance, -s)*cutoff(distance, cutoff_radius);
            }
        }
    }
    return energy;
}


int minimizeEnergy::Point2Cube(const cppoptlib::Vector<double> & thisPt){
    int tmp = 0, output = 0;
    for(int i = 0; i < dim; i++){
        tmp = floor((thisPt(i) +1.0) / cutoff_radius);
        output += tmp * pow(cubes_per_side, i);
    }
    return output;
}


void minimizeEnergy::FindNeighborCubes(int cube_index, cppoptlib::Vector<int> & neighbors){
    
    cppoptlib::Vector<int> current_index(dim);
    //neighbors = Eigen::MatrixXd::Constant(neighbors.size(), 1, -1);
    int tmp = cube_index;
    
    for (int i = 0; i < dim; ++i){
        current_index(i) = tmp % cubes_per_side;
        tmp = tmp / cubes_per_side;
    }
    
    for (int i = 0; i < pow(3, dim); i++){
        int step = i;
        int flag = 0, index = 0;
        // index is the index of neighboring cube in the matrix of cubes;
        
        for (int j = 0; j < dim && flag == 0; j++){
            int k = step%3 - 1 + current_index[j];
            step = step/3;
            if (k < 0 || k >= cubes_per_side) flag=1;
            index += k * pow(cubes_per_side, j);
        }
        if (flag == 0) neighbors[i] = index;
    }
}


void minimizeEnergy::BuildIndex(const cppoptlib::Vector<double> & x){
// get the 3d points matrix -- pts
// assign points to cubes
// the first field in each column is the number of points inside that column.
    Cubes.row(0).setZero();
    
    cppoptlib::Vector<double> temp_vector(dim);
    for (int i=0; i<numpts; ++i)
    {
        To3D(x.segment((dim-1)*i, dim-1), temp_vector);
        ptsND.col(i) = temp_vector.transpose();
        int max_neighbor = Cubes.rows();
        int cube_index = Point2Cube(ptsND.col(i));
        int tmp_numpts = Cubes(0, cube_index) + 1;
        if (tmp_numpts >= max_neighbor) {
            cout << "Warning: exceeding maximum neighbor; ignore point " << i << endl;
        }else{
            Cubes(tmp_numpts, cube_index) = i;
            Cubes(0, cube_index) = tmp_numpts;
        }
    }
    //cout << Cubes.transpose() << endl;
    //cout << Point2Cube(ptsND.col(4)) << endl;
}






// generate random sphere configuration
void randptSphere(double coordinates[], int dim){
    
    double z;
    double norm;
    double normsq=2;
    
    while(normsq>1 || normsq==0){
        normsq=0;
        
        for(int i=0;i<dim;i++){
            z=1-(2*(double)rand()/(double)RAND_MAX);
            normsq += z*z;
            coordinates[i] = z;
        }
    }
    
    norm=sqrt(normsq);
    
    for(int i=0;i<dim;i++){
        coordinates[i] = coordinates[i]/norm;
    }
    
}



/////////////////////////////////////////////////////////////////////////////////
int main() {
    
    
    cppoptlib::Vector<double> test1(2), test2(3);
    test1[0] = PI/2;
    test1[1] = PI/2;
    ToND(test1, test2);
    cout << test2 << endl;
    
    
    
    
    
    // declare all the parameter used.
//    double s, radius;
//    int dim = 0, numpts=0, c=0, cubes_per_side=0, max_neighbor=0, numFile=0, numIteration=0;
//    bool infile;
//    ifstream inputfile, pointfile;
//    ofstream outputfile;
//
//    
//    openFile(inputfile, "control.inp");
//    string filename = ParseControlFile(inputfile, dim, numpts, s, c, max_neighbor, numFile, numIteration, infile);
//    inputfile.close();
//    
//    
//    
//    radius = c*pow(numpts,-1.0/(dim-1));
//    cubes_per_side = ceil(2/radius);
//    
//    
//    
//    
//    cppoptlib::Matrix<double> X(numpts, dim), A(numpts, dim-1);
//    cppoptlib::Vector<double> V(A.size()), G(A.size()), GFull(A.size());
//    
//    // read points
//    if (infile)
//    {
//        openFile(pointfile, filename);
//        int lineNumber = 0;
//        while (!pointfile.eof() && lineNumber < numpts)
//        {
//            for (int i=0; i<dim;  ++i) pointfile >> X(lineNumber, i);
//            lineNumber++;
//        }
//        pointfile.close();
//    }
//    // generate random configuration
//    else
//    {
//        srand(time(0));
//        double apoint[dim];
//        for (int i = 0; i < numpts; i++) {
//            randptSphere(apoint, dim);
//            for (int j = 0; j < dim; j++) {
//                X(i, j) = apoint[j];
//            }
//        }
//    }
//    
//    
//    
//
//    ToAngles(X,A);
//    ToVector(A,V);
//    
//    minimizeEnergy f(radius, s, dim, numpts, cubes_per_side, pow(cubes_per_side, dim), max_neighbor);
//    cout << "Energy without cutoff:      " << Energy(V, s, dim-1) << endl;
//    cout << "Energy with cutoff:         " << f(V) << endl;
//
//    
//    cppoptlib::LbfgsSolver<double> solver;
//    
//    solver.setNumFile(numFile);
//    solver.setNumIteration(numIteration);
////    solver.setFileName(filename);
//    
//    solver.minimize(f, V);
//    
//    
//    cout << "Energy now: " << f(V) << endl;
//    cout << "Full energy now: " << Energy(V, s, dim-1) << endl;
//    
//    
//    
//    writeFile(outputfile, "output.txt", V, dim);
//    return 0;
}











