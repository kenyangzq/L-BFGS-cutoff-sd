// dist_squared
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

#define PI 3.141592653589793

using namespace std;


void openFile (ifstream & inputfile, string name){
    inputfile.open(name.c_str());
    if (inputfile.fail()) {
        cout << "Error opening input data file\n";
        exit(1);
    }
}

string ParseControlFile(ifstream & inputfile, int & dim, int &numpts, double & s, int & c, 
		int & max_neighbor, int & numFile, int & numIteration, bool & infile){
    
    
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
    for (int i = 1; i < coords.size()-1; i++) {
	   	tmp *= sin(angles[i-1]);
        coords[i] = cos(angles[i])*tmp;
//        cout << coords[i] << endl;
    }
    coords[coords.size()-1] = tmp*sin(angles[angles.size()-1]);
}


// angels are indexed 0 to dim-2
void ToAngle(cppoptlib::Vector<double> & angels, const cppoptlib::Vector<double> & coords) {
    int dim = coords.size();
    double squaresum = pow(coords[dim - 1],2);
    
    for (int i = dim - 2; i >= 0; i--) {
        squaresum += pow(coords[i], 2);
        angels[i] = acos(coords[i]/sqrt(squaresum));
    }
    
    if (coords[dim-1] < 0) {
        angels[dim-2] = 2 * PI - angels[dim-2];
    }
}

void ComputeJacobianN(const cppoptlib::Vector<double> &angles, cppoptlib::Matrix<double> &temp) {
	// all angles range from 0 to PI except for the last one, which range from 0 to 2PI.
	// The Jacobian Matrix should be size n by n-1;
	// Multiply by a -tan value if it's differentiating cos
	// Multiply by a cot value if it's differentiating sin
	//
	// Note that the Jacobian Matrix in this case is lower diagonal
	//
	int n = angles.size() + 1;
	temp = cppoptlib::Matrix<double>::Zero(n,n-1);
	cppoptlib::Vector<double> ndvector(n);
	ToND(angles,ndvector);
	// main loop
	for (int i = 0; i < n-1; i++) 
    {
		// calculate the diagonal value
		temp(i,i) = ndvector(i)*(-tan(angles(i)));
		for(int j = 0; j < i; j++) 
        {
			// calculate the part in the lower diagonal matrix
			// except the last row of the matrix
			temp(i,j) = ndvector(i)/tan(angles(j));
		}
		// do the last row of the matrix here, separately
		temp(n-1,i) = ndvector(n-1)/tan(angles(i));
	}
}
		


double dist_squared(const cppoptlib::Vector<double> & angles1,const cppoptlib::Vector<double> & angles2)
{
    cppoptlib::Vector<double> temp2(angles2.size()+1), temp1(angles1.size()+1);
    ToND(angles2, temp2);
    ToND(angles1, temp1);
    return (temp2- temp1).squaredNorm();
}



//void To3D(const cppoptlib::Vector<double> & angles, cppoptlib::Vector<double> & coords)
//{
//    coords(0) = cos(angles(0)) * sin(angles(1));
//    coords(1) = sin(angles(0)) * sin(angles(1));
//    coords(2) = cos(angles(1));
//}
//void ComputeJacobian(const double & theta, const double & phi, cppoptlib::Matrix<double> & temp){
//    //x = sin(phi) cos(theta)
//    //y = sin(phi) sin(theta)
//    //z = cos(phi)
//    //
//    // derivatives w.r.t. theta:
//    temp(0,0) = -sin(phi) * sin(theta); // x
//    temp(1,0) =  sin(phi) * cos(theta); // y
//    temp(2,0) =  0;                      // z
//    // derivatives w.r.t. phi:
//    temp(0,1) =  cos(phi) * cos(theta); // x
//    temp(1,1) =  cos(phi) * sin(theta); // y
//    temp(2,1) = -sin(phi);             // z
//}









void ToAngles(cppoptlib::Matrix<double> & all_points, cppoptlib::Matrix<double> & all_angles, int dim)
{
 
	cppoptlib::Vector<double> apoint(dim);
    cppoptlib::Vector<double> aangle(dim-1);
	for (int i = 0; i < all_points.rows(); i++) {
        apoint = all_points.row(i);
        ToAngle(aangle, apoint);
        all_angles.row(i) = aangle;
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
        ToND(tmp, tmp2);
        
        
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
    
    
	// need to debug gradient here
	// get to infinite when execute
	// doubt that cutoff prime or cutoff should change as well
	//
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
                temp_jacobian.setZero();
                ComputeJacobianN(x.segment(point_index*dim_angle,dim_angle), temp_jacobian);
                
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
        ToND(x.segment((dim-1)*i, dim-1), temp_vector);
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
    
//	int  dim = 20; 
//    cppoptlib::Vector<double> test1(dim), test2(dim);
//    cppoptlib::Vector<double> vec(dim+1);
//    test1 = PI/4.0 * cppoptlib::Vector<double>::Ones(20);
//    test2 =  3.0*PI/4.0 * cppoptlib::Vector<double>::Ones(20);
//    test2[dim-1] = 5.0*PI/4.0;
//    cout << test1.transpose() << endl;
//    cout << test2.transpose() << endl;
//    cout << dist_squared(test1, test2) << endl;

//    ToND(test2, vec);
//    ToAngle(test1, vec);
//    cout << (test2-test1).norm() << endl;
    Eigen::MatrixXf m(4,4);
    m <<  1, 2, 3, 4,
        5, 6, 7, 8,
        9,10,11,12,
       13,14,15,16;
//    cout << m << endl;
    cout << m(3,2) << endl;
    
	// test the ComputeJacobianN
	cppoptlib::Vector<double> ang1(5);
	ang1 << PI/2, PI/5, PI/7, PI/3, 3*PI/8;
	cppoptlib::Matrix<double> jac(6,5);
	ComputeJacobianN(ang1,jac);
	cout << jac << endl;
	cout << "\n\n\n\n"; 

	//test ends
    
    
    
    // declare all the parameter used.
    double s, radius;
    int dim = 0, numpts=0, c=0, cubes_per_side=0, max_neighbor=0, numFile=0, numIteration=0;
    bool infile;
    ifstream inputfile, pointfile;
    ofstream outputfile;

    
    openFile(inputfile, "control.inp");
    string filename = ParseControlFile(inputfile, dim, numpts, s, c, max_neighbor, numFile, numIteration, infile);
    inputfile.close();
    
    
    
    radius = c*pow(numpts,-1.0/(dim-1));
    cubes_per_side = ceil(2/radius);
    
    
    
    
    cppoptlib::Matrix<double> X(numpts, dim), A(numpts, dim-1);
    cppoptlib::Vector<double> V(A.size()), G(A.size()), GFull(A.size());
    
    // read points
    if (infile)
    {
        openFile(pointfile, filename);
        int lineNumber = 0;
        while (!pointfile.eof() && lineNumber < numpts)
        {
            for (int i=0; i<dim;  ++i) pointfile >> X(lineNumber, i);
            lineNumber++;
        }
        pointfile.close();
    }
    // generate random configuration
    else
    {
        srand(time(0));
        double apoint[dim];
        for (int i = 0; i < numpts; i++) {
            randptSphere(apoint, dim);
            for (int j = 0; j < dim; j++) {
                X(i, j) = apoint[j];
            }
        }
    }
    
    ToAngles(X, A, dim);
    
    
    ToVector(A,V);
    
    minimizeEnergy f(radius, s, dim, numpts, cubes_per_side, pow(cubes_per_side, dim), max_neighbor);
    cout << "Energy without cutoff:      " << Energy(V, s, dim-1) << endl;
    cout << "Energy with cutoff:         " << f(V) << endl;

    
    cppoptlib::LbfgsSolver<double> solver;
    
    solver.setNumFile(numFile);
    solver.setNumIteration(numIteration);
//    solver.setFileName(filename);
    
    solver.minimize(f, V);
    
    
    cout << "Energy now: " << f(V) << endl;
    cout << "Full energy now: " << Energy(V, s, dim-1) << endl;
    
    
    
    writeFile(outputfile, "output.txt", V, dim);
    return 0;
}











