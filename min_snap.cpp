#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "mosek.h"
#include <nlohmann/json.hpp>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

// (time, position)
struct Waypoint1D { 
  double t, x; 
};

// (time, x, y)
struct Waypoint2D { 
  double t, x, y; 
};

static double factorial(int k) {
  double f = 1.0;
  for(int i = 1; i <= k; ++i) 
    f *= i;
  return f;
}

static void MSKAPI printstr(void*, const char *msg) {
  printf("%s", msg);
  fflush(stdout);
}

vector<double> solveMinSnap1D(const vector<Waypoint1D>& W) {
  int Mseg = int(W.size()) - 1;      // number of segments
  int d    = 7;                      // polynomial degree
  int n    = d + 1;                  // coefficients per segment
  int N    = Mseg * n;               // total decision variables
  int C    = 3 + 3 + 3*(Mseg-1) + (Mseg-1);

  // Build the Q matrix 
  MatrixXd Q = MatrixXd::Zero(N,N);
  for(int i = 0; i < Mseg; ++i) {
    double dt = W[i+1].t - W[i].t;
    MatrixXd Qi = MatrixXd::Zero(n,n);
    for(int p = 4; p <= d; ++p)
      for(int q = 4; q <= d; ++q) {
        double v = factorial(p)/factorial(p-4)
                 * factorial(q)/factorial(q-4)
                 / double(p+q-7)
                 * pow(dt, p+q-7);
        Qi(p,q) = v;
      }
    Q.block(i*n, i*n, n, n) = Qi;
  }

  // Build A x = b constraints
  MatrixXd A = MatrixXd::Zero(C,N);
  VectorXd b = VectorXd::Zero(C);
  int row = 0;

  // Boundary at t=0
  A(row,0) = 1;  b(row) = W[0].x; ++row;
  A(row,1) = 1;  b(row) = 0.0;    ++row;
  A(row,2) = 2;  b(row) = 0.0;    ++row;

  // Boundary at final time
  double DT = W.back().t - W[Mseg-1].t;
  for(int k = 0; k < 3; ++k) {
    for(int j = k; j <= d; ++j)
      A(row,(Mseg-1)*n + j) = factorial(j)/factorial(j-k) * pow(DT, j-k);
    b(row) = (k == 0 ? W.back().x : 0.0);
    ++row;
  }

  // 3) Continuity + position match at interior waypoints
  for(int i = 1; i < Mseg; ++i) {
    double dti = W[i].t - W[i-1].t;
    for(int k = 0; k < 3; ++k) {
      for(int j = k; j <= d; ++j)
        A(row,(i-1)*n + j) = factorial(j)/factorial(j-k) * pow(dti, j-k);
      A(row,i*n + k) = -factorial(k);
      b(row) = 0.0;
      ++row;
    }
    // position continuity
    A(row,i*n) = 1.0;
    b(row) = W[i].x;
    ++row;
  }

  // Quick feasibility check 
  {
    VectorXd sol = A.fullPivLu().solve(b);
    if((A*sol - b).norm() > 1e-8)
      throw runtime_error("1D constraints infeasible");
  }

  // Create MOSEK environment and task
  MSKenv_t  env  = nullptr;
  MSKtask_t task = nullptr;
  if(MSK_makeenv(&env, nullptr) != MSK_RES_OK ||
     MSK_maketask(env, C, N, &task) != MSK_RES_OK)
    throw runtime_error("Failed to create MOSEK task");

  MSK_putintparam(task, MSK_IPAR_LOG, 0);
  MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, nullptr, printstr);

  MSK_appendvars(task, N);
  MSK_appendcons(task, C);
  for(int j = 0; j < N; ++j)
    MSK_putvarbound(task, j, MSK_BK_FR, 0.0, 0.0);

  // Objective: minimize 1/2 xᵀ Q x
  MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  vector<double> cz(N, 0.0);
  MSK_putcslice(task, 0, N, cz.data());
  vector<MSKint32t> qi, qj;
  vector<double> qv;
  for(int i = 0; i < N; ++i)
    for(int j = 0; j <= i; ++j) {
      double v = Q(i,j);
      if(fabs(v) > 1e-12) {
        qi.push_back(i);
        qj.push_back(j);
        qv.push_back(v);
      }
    }
  MSK_putqobj(task, (MSKint64t)qv.size(), qi.data(), qj.data(), qv.data());

  // Load A x = b
  for(int rr = 0; rr < C; ++rr) {
    vector<MSKint32t> idx;
    vector<double>    aval;
    for(int j = 0; j < N; ++j) {
      if(fabs(A(rr,j)) > 1e-12) {
        idx.push_back(j);
        aval.push_back(A(rr,j));
      }
    }
    MSK_putarow(task, rr, (MSKint32t)idx.size(), idx.data(), aval.data());
    MSK_putconbound(task, rr, MSK_BK_FX, b(rr), b(rr));
  }

  // Solve QP
  MSKrescodee trmcode;
  if(MSK_optimizetrm(task, &trmcode) != MSK_RES_OK)
    throw runtime_error("MOSEK solve error");

  // Check solution status
  MSKsolstae solsta;
  MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
  if(solsta != MSK_SOL_STA_OPTIMAL)
    throw runtime_error("Solution not optimal");

  // Retrieve solution
  vector<double> xOpt(N);
  MSK_getxx(task, MSK_SOL_ITR, xOpt.data());

  // Cleanup
  MSK_deletetask(&task);
  MSK_deleteenv(&env);

  return xOpt;
}

int main() {
  vector<Waypoint2D> W2d = {
    {0.0, 0.0, 0.0},
    {1.0, 1.2, 2.1},
    {2.5, -0.5, 1.0},
    {4.0, 2.0, 0.5}
  };

  int Mseg = int(W2d.size()) - 1;
  int d    = 7, n = d + 1;

  vector<Waypoint1D> Wx, Wy;
  for(auto &wp : W2d) {
    Wx.push_back({wp.t, wp.x});
    Wy.push_back({wp.t, wp.y});
  }


  auto cx = solveMinSnap1D(Wx);
  auto cy = solveMinSnap1D(Wy);

  // Save coefficients to JSON
  json j;
  j["times"]    = vector<double>{0.0, 1.0, 2.5, 4.0};
  j["coeffs_x"] = cx;
  j["coeffs_y"] = cy;

  ofstream o("min_snap_2d_coeffs.json");
  o << j.dump(2) << endl;

  cout << "Wrote JSON to min_snap_2d_coeffs.json\n";
  return 0;
}
