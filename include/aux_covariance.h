#include<Eigen/Dense>
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

void print_covariance(Vector3d r, Vector3d v, MatrixXd P, ostream &os);
void print_DVs(double t, Vector3d r, Vector3d DV, ostream &os);