
#include <cstdlib>
#include <exception>

#include "lotka_volterra.hpp"

int main() {
  std::cout << "Digit 4 positive parameters A, B, C, D separeted by a space"
            << "\n";
  double newA, newB, newC, newD;
  std::cin >> newA >> newB >> newC >> newD;
  if (!std::cin) {
    std::cout << "Wrong data type" << "\n";
    std::exit(1);
  }
  if (newA <= 0 || newB <= 0 || newC <= 0 || newD <= 0) {
    std::cout << "Not all positive" << "\n";
    std::exit(0);
  }

  std::cout
      << "Choose the initial value of the preys (x_0) and the "
         "predators (y_0) separeted by a space. The must be positive numbers"
      << "\n";
  double x_0, y_0;
  std::cin >> x_0 >> y_0;
  if (!std::cin) {
    std::cout << "Wrong data type" << "\n";
    std::exit(1);
  }
  if (x_0 <= 0 || y_0 <= 0) {
    std::cout << "Not all positive or wrong data type" << "\n";
    std::exit(0);
  }

  std::cout
      << "Choose the duration (s) of the simulation. It must be a positive "
         "number"
      << "\n";
  double time;
  std::cin >> time;
  if (!std::cin) {
    std::cout << "Wrong data type" << "\n";
    std::exit(1);
  }
  if (time <= 0) {
    std::cout << "Not all positive or wrong data type" << "\n";
    std::exit(0);
  };
  pf::Simulation simulation(newA, newB, newC, newD, x_0, y_0, 0.001);
  simulation.e2_x();
  simulation.e2_y();

  int steps = static_cast<int>(
      time / 0.001);  // static_cast per convertire da un tipo all altro
  simulation.startingVectors();
  simulation.nextstep(steps);
  simulation.statistical_data();
  simulation.parameters_txt();
  simulation.simulation_txt();
  simulation.stats_txt();
  std::cout << "Simulation ended, data are in statistic_txt, "
               "value_iteration.txt and parameters_and_Coordinates.txt"
            << "\n";
  simulation.createOrbit();
  simulation.sinusoidal();
}
