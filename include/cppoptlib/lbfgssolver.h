// CppNumericalSolver
#include <iostream>
#include <Eigen/LU>
#include <iomanip>
#include <fstream>
#include "isolver.h"
#include "./linesearch/armijo.h"
#include "./linesearch/morethuente.h"
#include "./linesearch/wolfeheuristic.h"
#include <string>

#ifndef LBFGSSOLVER_H_
#define LBFGSSOLVER_H_




// Need to fix writeFiles 
// Need to convert to n d
// Fix search Direction or q







namespace std {
    void writeFiles (ofstream & outputfile, string name, cppoptlib::Vector<double> V, int dim){
        
        int c = dim-1;
        
        outputfile << setprecision(6);
        outputfile << fixed;
        
        cppoptlib::Vector<double> tmp2(dim);
        cppoptlib::Vector<double> tmp(c);
        
        for (int i =0; i < V.rows()/c; i++) {
            tmp = V.segment(i*c, c);

            tmp2(0) = cos(tmp(0)) * sin(tmp(1));
            tmp2(1) = sin(tmp(0)) * sin(tmp(1));
            tmp2(2) = cos(tmp(1));
            
            outputfile << tmp2(0) << "\t" << tmp2(1) << "\t" << tmp2(2) << "\n";
        }
    }
}







namespace cppoptlib {

template<typename T>
class LbfgsSolver : public ISolver<T, 1> {
  public:
    
    size_t numFile = 20;
    
    size_t numIteration = 1000;
    
    size_t m = 1;
    
//    std::string filename = "config";
    
    void setHistorySize(const size_t hs) { m = hs; }
    
    void setNumFile(const size_t nf) { numFile = nf; }
    
    void setNumIteration(const size_t i) { numIteration = i; }
    
//    void setFileName(const std::string fn) { filename = fn; }
    
    
    void minimize(Problem<T> &objFunc, Vector<T> & x0) {
        
        printf("===============================\n");
        printf("Minimization start: \n");
        
        std::ofstream outputfile;
    
        const size_t DIM = x0.rows();

        Matrix<T> sVector = Matrix<T>::Zero(DIM, m);
        Matrix<T> yVector = Matrix<T>::Zero(DIM, m);

        Vector<T> alpha = Vector<T>::Zero(m);
        Vector<T> grad(DIM), q(DIM), grad_old(DIM), s(DIM), y(DIM);
        objFunc.gradient(x0, grad);
		
         printf("grad's norm is: %f\n",grad.norm()); 
        Vector<T> x_old = x0;

        size_t iter = 0;
        double H0k = 1;
        int tmp = numIteration/numFile;

        this->m_current.reset();
        do {
            const T relativeEpsilon = static_cast<T>(0.0001) * std::max(static_cast<T>(1.0), x0.norm());

            if (grad.norm() < relativeEpsilon) {
                printf("Minimize complete.\n");
                break;
            }
            //Algorithm 7.4 (L-BFGS two-loop recursion)
			// Need to check why q fails 
			//
			q = grad;
            const int k = std::min(m, iter);
            
            // for i = k − 1, k − 2, . . . , k − m
            for (int i = k - 1; i >= 0; i--) {
                // alpha_i <- rho_i*s_i^T*q
                double product = static_cast<Vector<T>>(sVector.col(i))
                .dot(static_cast<Vector<T>>(yVector.col(i)));
			   	const double rho = 1.0 / product;
                alpha(i) = rho * static_cast<Vector<T>>(sVector.col(i)).dot(q);
                // q <- q - alpha_i*y_i
                q = q - alpha(i) * yVector.col(i);
            }
            // r <- H_k^0*q
            
            q = H0k * q;
            
            //for i k − m, k − m + 1, . . . , k − 1
            for (int i = 0; i < k; i++) {
                // beta <- rho_i * y_i^T * r
                double product = static_cast<Vector<T>>(sVector.col(i))
                .dot(static_cast<Vector<T>>(yVector.col(i)));
                const double rho = 1.0 / product;
                const double beta = rho * static_cast<Vector<T>>(yVector.col(i)).dot(q);
                // r <- r + s_i * ( alpha_i - beta)
                q = q + sVector.col(i) * (alpha(i) - beta);
            }
            // stop with result "H_k*f_f'=q"

            // any issues with the descent direction ?
            double descent = -grad.dot(q);
            
            double alpha_init =  1.0 / grad.norm();
            
            if (descent > -0.0001 * relativeEpsilon) {
                q = -1 * grad;
                iter = 0;
                alpha_init = 1.0;
            }
	
            // find steplength
            WolfeHeuristic<T, decltype(objFunc), 1>::linesearch(x0, -q,  objFunc, alpha_init) ;
            // update guess
            
//            x0 = x0 - rate * q;
            

            grad_old = grad;
            objFunc.gradient(x0, grad);

            s = x0 - x_old;
            y = grad - grad_old;

            // update the history
            if (iter < m) {
                sVector.col(iter) = s;
                yVector.col(iter) = y;
            } else {
                sVector.leftCols(m - 1) = sVector.rightCols(m - 1).eval();
                sVector.rightCols(1) = s;
                yVector.leftCols(m - 1) = yVector.rightCols(m - 1).eval();
                yVector.rightCols(1) = y;
            }
            
            
            
            // update the scaling factor
            double dot =  static_cast<double>(y.dot(y));
            if (dot <= 1e-7) {
                
                printf("\nMinimize complete.\n");
                printf("Norm between two consecutive iterations is %.10f\n", s.norm());
                printf("Now the norm of gradient is: %f\n", grad.norm());
                printf("===============================\n");
                break;
            }else
                H0k = y.dot(s) / dot;
            
            
            
            x_old = x0;
            
            
            // std::cout << "iter: "<<globIter<< ", f = " <<  objFunc.value(x0) << ", ||g||_inf "
            // <<gradNorm  << std::endl;
            iter++;
            
//            std::cout << "iter: " << iter << std::endl;
            
            
            if (iter % tmp == 0) {
                std::string name = "./output/output" + std::to_string(iter/tmp) + ".txt";
                outputfile.open(name.c_str());
                
                printf("Energy: %f\n", objFunc.value(x0));
                // writeFiles(outputfile, name, x0, 3);
                
                outputfile.close();
            }
            printf("Point 4 \n");
            
            ++this->m_current.iterations;
            this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
            this->m_status = checkConvergence(this->m_stop, this->m_current);
        } while (this->m_status == Status::Continue && iter < numIteration);

    }

};

}
/* namespace cppoptlib */

#endif /* LBFGSSOLVER_H_ */
