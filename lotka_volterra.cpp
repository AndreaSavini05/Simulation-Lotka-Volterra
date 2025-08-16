#include "lotka_volterra.hpp"

namespace pf {

// costruttore con parametri iniziali
Simulation::Simulation(double newA, double newB, double newC, double newD,
                       double newx_0, double newy_0, double newdt)
    : A{newA}, B{newB}, C{newC}, D{newD}, x_0{newx_0}, y_0{newy_0}, dt{newdt} {}

// Getter degli elementi dei vector
const std::vector<double>& Simulation::gett() const { return t; }
const std::vector<double>& Simulation::getValuex() const { return data.x; }
const std::vector<double>& Simulation::getValuey() const { return data.y; }
const std::vector<double>& Simulation::getValueH() const { return data.H; }

double Simulation::e2_x() { return D / C; }  // coordinata x p.to equilibrio 2
double Simulation::e2_y() { return A / B; }  // coordinata y p.to equilibrio 2

// inizializzare i vettori x, y, H
void Simulation::startingVectors() {
  data.x.push_back(x_0);
  data.y.push_back(y_0);
  data.H.push_back(std::abs(-D * log(x_0) + C * x_0 + B * y_0 - A * log(y_0)));

  datarel.x_rel.push_back(x_0 * C / D);
  datarel.y_rel.push_back(y_0 * B / A);
  //  datarel.h_wsgn.push_back(-D * log(x_0) + C * x_0 + B * y_0 - A *
  //  log(y_0));
}

// metodo evolve
void Simulation::evolve() {
  double x_0_rel = x_0 * C / D;
  double y_0_rel = y_0 * B / A;

  double x_i_rel = x_0_rel + (A - B * y_0_rel) * x_0_rel * dt;
  double y_i_rel = y_0_rel + (C * x_0_rel - D) * y_0_rel * dt;

  double x_i = x_i_rel * D / C;
  double y_i = y_i_rel * A / B;

  assert(x_i_rel > 0 && "x_i  must be > 0, otherwise prey  extincted");
  assert(y_i_rel > 0 && " y_i must be > 0, otherwise predator extincted");
  /* if (x_i <= 0.0 || y_i <= 0.0) {//provsre con assert e throw
     std::cerr << "Errore: una delle due specie si è estinta. Simulazione "
                  "interrotta.\n";
     return;
   }*/

  // double h_wsgn_i = -D * std::log(x_i) + C * x_i + B * y_i - A *
  // std::log(y_i);
  double H_i = std::abs(-D * log(x_0) + C * x_0 + B * y_0 - A * log(y_0));

  data.x.push_back(x_i);
  data.y.push_back(y_i);
  data.H.push_back(H_i);

  datarel.x_rel.push_back(x_i_rel);
  datarel.y_rel.push_back(y_i_rel);
  //  datarel.h_wsgn.push_back(h_wsgn_i);

  x_0 = x_i;
  y_0 = y_i;

  assert(!data.x.empty());
  assert(!data.y.empty());
  assert(!data.H.empty());
  assert(!datarel.x_rel.empty());
  assert(!datarel.y_rel.empty());
  // assert(!datarel.h_wsgn.empty());
}

void Simulation::nextstep(int n) {  // provare a usare algoritmo
  for (int i = 1; i <= n; ++i) {
    evolve();
    t.push_back(dt * i);
  }
};
void Simulation::statistical_data() {
  if (data.x.empty() || data.y.empty() || data.H.empty()) {
    std::cout << "No data available.\n";
    return;
  }

  stats.mean_x = std::accumulate(data.x.begin(), data.x.end(), 0.0) /
                 static_cast<double>(data.x.size());
  stats.mean_y = std::accumulate(data.y.begin(), data.y.end(), 0.0) /
                 static_cast<double>(data.y.size());
  stats.mean_H = (std::accumulate(data.H.begin(), data.H.end(), 0.0) /
                  static_cast<double>(data.H.size()));

  auto [min_x, max_x] = std::minmax_element(
      data.x.begin(),
      data.x.end());  // structured binding auto [min_y, max_y] =
  auto [min_y, max_y] = std::minmax_element(data.y.begin(), data.y.end());
  auto [min_H, max_H] = std::minmax_element(data.H.begin(), data.H.end());

  stats.min_x = *min_x;
  stats.max_x = *max_x;
  stats.min_y = *min_y;
  stats.max_y = *max_y;
  stats.min_H = *min_H;
  stats.max_H = *max_H;
};

void Simulation::parameters_txt() {
  std::ofstream file("parameters_and_Coordinates.txt");
  if (file.is_open()) {
    file << "(" << e2_x() << "; " << e2_y() << ");" << "\n";
    file << "\n";
    file << "Parametri: " << "\n";
    file << "A " << A << "\n";
    file << "B " << B << "\n";
    file << "C " << C << "\n";
    file << "D " << D << "\n";
    file << "\n";
    file << "prede: " << data.x[0] << "; predatori: " << data.y[0] << "\n";

  } else {
    std::cerr << "Impossibile aprere il file" << "\n";
  }
  file.close();
}

void Simulation::simulation_txt() const {
  std::ofstream file("value_iteration.txt");
  if (file.is_open()) {
    std::size_t i = 0;
    std::for_each(t.begin(), t.end(), [&](const auto& t_val) {
      file << t_val << " " << data.x[i] << " " << data.y[i] << " " << data.H[i]
           << "\n";
      ++i;
    });
  } else {
    std::cerr << "Impossibile aprere il file" << "\n";
  }
  file.close();
}

void Simulation::stats_txt() const {
  if (stats.mean_x == 0.0 && stats.mean_y == 0.0 && stats.mean_H == 0.0) {
    std::cout << "Stats not calculated. Call statistical_data() first.\n";
    return;
  }
  std::ofstream file("statistic.txt");
  if (file.is_open()) {
    file << "Prede: " << "\n";
    file << "x_mean: " << stats.mean_x << "; x_min: " << stats.min_x
         << "; x_max: " << stats.max_x << "\n";
    file << "\n";
    file << "Predatori: " << "\n";
    file << "y_mean: " << stats.mean_y << "; y_min: " << stats.min_y
         << "; y_max: " << stats.max_y << "\n";
    file << "\n";
    file << "Integrale primo: " << "\n";
    file << "H_mean: " << stats.mean_H << "; H_min: " << stats.min_H
         << "; H_max: " << stats.max_H << "\n";
  } else {
    std::cerr << "Impossibile aprere il file" << "\n";
  }
  file.close();
};

void Simulation::createOrbit() {
  unsigned int width = 900.0f;
  unsigned int height = 775.0f;
  float edge = 25.0f;
  // float y_edge = 21.5f;
  float edge_ax = 20.0f;
  // float y_edge_ax = 17.2f;
  //  float x_factor = 360.0f;
  // float y_factor = 310.f;
  sf::RenderWindow window(sf::VideoMode(width, height), "Orbit graph",
                          sf::Style::Default);
  window.setFramerateLimit(60);

  sf::Vertex x_axis[2] = {
      sf::Vertex(sf::Vector2f(0.0f, static_cast<float>(height) - edge_ax),
                 sf::Color::White),
      sf::Vertex(sf::Vector2f(static_cast<float>(width),
                              static_cast<float>(height) - edge_ax),
                 sf::Color::White)};

  sf::Vertex y_axis[2] = {
      sf::Vertex(sf::Vector2f(edge_ax, 0.0f), sf::Color::White),
      sf::Vertex(sf::Vector2f(edge_ax, static_cast<float>(height)),
                 sf::Color::White)};
  // da Chat
  double minx = std::min(stats.min_x, e2_x());
  double maxx = std::max(stats.max_x, e2_x());
  double miny = std::min(stats.min_y, e2_y());
  double maxy = std::max(stats.max_y, e2_y());

  float x_factor =
      (static_cast<float>(width) - 2 * edge) / static_cast<float>(maxx - minx);
  float y_factor =
      (static_cast<float>(height) - 2 * edge) / static_cast<float>(maxy - miny);
  /*
    sf::VertexArray orbit(
        sf::LineStrip,
        data.x.size());  // LinesStrip    List of connected lines, a point uses
                         // the previous point to form a line.
    for (size_t i = 0; i < data.x.size(); ++i) {
      float px = static_cast<float>(x_margin + (data.x[i] - stats.min_x) /
                                                   (stats.max_x - stats.min_x) *
                                                   x_factor);
      //  static_cast<float>(width));

      float py = static_cast<float>(height) - y_margin -
                 static_cast<float>(
                     ((data.y[i] - stats.min_y) / (stats.max_y - stats.min_y)) *
                     y_factor);  /
                                 // static_cast<float>(height));
      orbit[i].position = sf::Vector2f(px, py);
      orbit[i].color = sf::Color::White;
    }*/
  // da ChatGPT

  sf::VertexArray orbit(sf::LineStrip, data.x.size());
  for (size_t i = 0; i < data.x.size(); ++i) {
    float px = edge + static_cast<float>((data.x[i] - minx) * x_factor);
    float py = static_cast<float>(height) - edge -
               static_cast<float>((data.y[i] - miny) * y_factor);
    orbit[i].position = sf::Vector2f(px, py);
    orbit[i].color = sf::Color::White;
  }
  sf::CircleShape equilibrium_point(5.0f);
  equilibrium_point.setOrigin(5.0f, 5.0f);
  equilibrium_point.setFillColor(sf::Color::Red);

  float eq_x = edge + static_cast<float>((e2_x() - minx) * x_factor);
  float eq_y = static_cast<float>(height) - edge -
               static_cast<float>((e2_y() - miny) * y_factor);
  equilibrium_point.setPosition(eq_x, eq_y);

  /*
    sf::CircleShape equilibrium_point(3.0f);
    equilibrium_point.setOrigin(3.0f, 3.0f);  // centra l'origine
    equilibrium_point.setFillColor(sf::Color::Red);
    equilibrium_point.setFillColor(sf::Color::Red);
    equilibrium_point.setPosition(
        static_cast<float>(x_margin) +
            static_cast<float>(
                std::abs(((e2_x() - stats.min_x)) / (stats.max_x - stats.min_x))
    * x_factor), static_cast<float>(height) - y_margin - static_cast<float>(
               std::abs(((e2_y() - stats.min_y)) / (stats.max_y - stats.min_y))
    * y_factor));*/

  sf::Font font;
  if (!font.loadFromFile("arial.ttf")) {
    std::cerr << ".Error while loading font" << '\n';
    return;
  }

  sf::Text predatori("asse y(predatori)", font, 17);
  predatori.setFillColor(sf::Color::White);
  predatori.setRotation(-90.0f);
  predatori.setPosition(1.0f, 160.0f);

  sf::Text prede("assex x(prede)", font, 17);
  prede.setFillColor(sf::Color::White);
  prede.setPosition(static_cast<float>(width) - 150.0f,
                    static_cast<float>(height) - 18.5f);

  sf::Text punto_eq("punto_eq", font, 17.0f);
  punto_eq.setFillColor(sf::Color::Red);
  punto_eq.setPosition(static_cast<float>(width) - 100.0f, 10.0f);

  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) window.close();
    }

    // nel ciclo di rendering:
    window.clear();
    window.draw(equilibrium_point);
    window.draw(x_axis, 2, sf::Lines);
    window.draw(y_axis, 2, sf::Lines);
    window.draw(orbit);
    window.draw(punto_eq);
    window.draw(prede);
    window.draw(predatori);
    window.display();
  }
}

void Simulation::sinusoidal() {
  float width = 900.0f;
  float height = 775.0f;
  float edge = 20.0f;

  sf::RenderWindow window(sf::VideoMode(static_cast<unsigned int>(width),
                                        static_cast<unsigned int>(height)),
                          "population during simulation", sf::Style::Default);
  window.setFramerateLimit(60);

  sf::Font font;
  if (!font.loadFromFile("arial.ttf")) {
    std::cerr << ".Error while loading font" << '\n';
    return;
  }

  sf::Vertex x_axis[2]{
      sf::Vertex(sf::Vector2f(0.0f, static_cast<float>(height) - edge),
                 sf::Color::White),
      sf::Vertex(sf::Vector2f(width, static_cast<float>(height) - edge),
                 sf::Color::White),

  };
  sf::Vertex y_axis[2]{
      sf::Vertex(sf::Vector2f(edge, 0.0f), sf::Color::White),
      sf::Vertex(sf::Vector2f(edge, height), sf::Color::White)};

  double t_min = t[0];
  double t_max = static_cast<double>(t.size() - 1) * dt;
  /*assert(t.size() == data.x.size());
  assert(t.size() == data.y.size());*/

  sf::VertexArray prede(sf::LineStrip, data.x.size());
  sf::VertexArray predatori(sf::LineStrip, data.y.size());
  for (size_t i = 0; i < data.x.size(); ++i) {
    float px = edge + ((static_cast<float>(i) * static_cast<float>(dt) -
                        static_cast<float>(t_min)) /
                       static_cast<float>(t_max - t_min) * (width - 2 * edge));
    float qx = height - edge -
               (static_cast<float>((data.x[i] - stats.min_x) /
                                   (stats.max_x - stats.min_x)) *
                (height - 2 * edge));

    float qy = height - edge -
               (static_cast<float>((data.y[i] - stats.min_y) /
                                   (stats.max_y - stats.min_y)) *
                (height - edge));
    predatori[i].position = sf::Vector2f(px, qy);
    predatori[i].color = sf::Color::Magenta;
    prede[i].position = sf::Vector2f(px, qx);
    prede[i].color = (sf::Color::Yellow);
  };

  sf::Text time("time", font, 14);
  time.setFillColor(sf::Color::White);
  time.setPosition(width - 50.0f, height - 20.0f);

  sf::Text population("population", font, 14);
  population.setRotation(-90.0f);
  population.setFillColor(sf::Color::White);
  population.setPosition(5.0f, 100.0f);

  sf::Text prey("prede", font, 14);
  prey.setFillColor(sf::Color::Yellow);
  prey.setPosition(50.0f, 20.0f);

  sf::Text predator("predator", font, 14);
  predator.setFillColor(sf::Color::Magenta);
  predator.setPosition(50.0f, 40.0f);

  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) window.close();
    }
    window.clear();
    window.draw(x_axis, 2, sf::Lines);
    window.draw(y_axis, 2, sf::Lines);
    window.draw(prede);
    window.draw(predatori);
    window.draw(time);
    window.draw(population);
    window.draw(prey);
    window.draw(predator);
    window.display();
  }
}

}  // namespace pf
