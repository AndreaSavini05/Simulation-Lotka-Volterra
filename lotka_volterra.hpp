#ifndef PF_LOTKA_VOLTERRA_HPP
#define PF_LOTKA_VOLTERRA_HPP
//#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>  //per moduli (abs)
#include <fstream>
#include <iostream>
#include <numeric>//std::accumultate
#include <vector>
#include <SFML/Graphics.hpp>

namespace pf {

struct DataSim {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> H;
};

struct DataPopulRel {
  std::vector<double> x_rel;
  std::vector<double> y_rel;
  //  std::vector<double> h_wsgn;
};
 struct Stats {
    double mean_x;
    double mean_y;
    double mean_H;
    double min_x, max_x;
    double min_y, max_y;
    double min_H, max_H;
  };


class Simulation {
  double A{};
  double B{};
  double C{};
  double D{};
  double x_0{};
  double y_0{};
  const double dt{0.001};

  DataSim data{};  // creazione variabile menbro privata chiamata data di tipo
                   // DataSim
  DataPopulRel datarel{};
  Stats stats{};
  std::vector<double> t = {};

 public:
  Simulation(double newA, double newB, double newC, double newD, double newx_0,
             double newy_0, double newdt);

  /*void setValueA(double newA) { A = newA; }
  void setValueB(double newB) { B = newB; }
  void setValueC(double newC) { C = newC; }
  void setValueD(double newD) { D = newD; }
  void setx_0(double newx_0) { x_0 = newx_0; }
  void sety_0(double newy_0) { y_0 = newy_0; }*/

  // metodo per coordinate punto di equiloibrio non banale
  double e2_x();
  double e2_y();

  void evolve();

  // metodo per iniziallizare i vettoro con x_0, y_0, H
  void startingVectors();

  /*metodo vector relativi
    void startingRelVectors();*/

  // metodo per far contunuare la simulazione
  void nextstep(int n);

  // methods to get elements inside the vectors (useful for testing and graphic)
  const std::vector<double>& getValuex() const;  // capire questa parte, & const
  const std::vector<double>& getValuey() const;  // tipo di ritorno: const std:....
  const std::vector<double>& getValueH() const;  // metodo getValue() è const
  const std::vector<double>& gett() const;

  void statistical_data();
  // scrivere su txt
  void parameters_txt();
  void simulation_txt() const;
  void stats_txt() const;
  void createOrbit();
  void sinusoidal();
};
}  // namespace pf

#endif
