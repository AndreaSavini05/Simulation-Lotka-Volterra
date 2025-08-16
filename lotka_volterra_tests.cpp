#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "lotka_volterra.hpp"

TEST_CASE("First set of parameters and initial value") {  // foglio 1 excell
  // setting initial values and paramteteres
  pf::Simulation simulation1(0.7, 0.35, 0.6, 1.2, 50., 40., 0.001);

  simulation1.startingVectors();

  // subcases
  SUBCASE("testing values of equilibrium") {
    CHECK(simulation1.e2_x() == 2);
    CHECK(simulation1.e2_y() == 2);
  }
  SUBCASE("testing initial values of x,y,H") {
    CHECK(simulation1.getValuex()[0] == 50);
    CHECK(simulation1.getValuey()[0] == 40);
    CHECK(simulation1.getValueH()[0] == doctest::Approx(36.72336));

    CHECK(simulation1.getValuex().size() == 1);
    CHECK(simulation1.getValuey().size() == 1);
    CHECK(simulation1.getValueH().size() == 1);
  }
  SUBCASE("testing valus after one step") {
    simulation1.nextstep(1);
    CHECK(simulation1.getValuex()[1] == doctest::Approx(49.335));
    CHECK(simulation1.getValuey()[1] == doctest::Approx(41.152));
    CHECK(simulation1.getValueH()[1] == doctest::Approx(36.72375));

    CHECK(simulation1.getValuex().size() == 2);
    CHECK(simulation1.getValuey().size() == 2);
    CHECK(simulation1.getValueH().size() == 2);
  }
  SUBCASE("testing values after 0.005") {
    simulation1.nextstep(5);
    CHECK(simulation1.getValuex()[5] == doctest::Approx(46.56875));
    CHECK(simulation1.getValuey()[5] == doctest::Approx(45.92001));
    CHECK(simulation1.getValueH()[5] == doctest::Approx(36.72531));

    CHECK(simulation1.getValuex().size() == 6);
    CHECK(simulation1.getValuey().size() == 6);
    CHECK(simulation1.getValueH().size() == 6);
  }
  SUBCASE("testing dimension of vectors after 0.010") {
    simulation1.nextstep(10);
    CHECK(simulation1.getValuex()[10] == doctest::Approx(42.90816));
    CHECK(simulation1.getValuey()[10] == doctest::Approx(52.17563));
    CHECK(simulation1.getValueH()[10] == doctest::Approx(36.72726));

    CHECK(simulation1.getValuex().size() == 11);
    CHECK(simulation1.getValuey().size() == 11);
    CHECK(simulation1.getValueH().size() == 11);
  }

  SUBCASE("testing dimension of vectors after 7") {
    simulation1.nextstep(7000);
    CHECK(simulation1.getValuex()[7000] == doctest::Approx(2.71e-13));
    CHECK(simulation1.getValuey()[7000] == doctest::Approx(0.029687));
    CHECK(simulation1.getValueH()[7000] == doctest::Approx(37.19687));

    CHECK(simulation1.getValuex().size() == 7001);
    CHECK(simulation1.getValuey().size() == 7001);
    CHECK(simulation1.getValueH().size() == 7001);
  }

  SUBCASE("testing dimension of vectors after 21.275") {
    simulation1.nextstep(7000);
    CHECK(simulation1.getValuex()[21275] == doctest::Approx(5.85e-9));
    CHECK(simulation1.getValuey()[21275] == doctest::Approx(1.07e-9));
    CHECK(simulation1.getValueH()[21275] == doctest::Approx(37.20826));

    CHECK(simulation1.getValuex().size() == 21276);
    CHECK(simulation1.getValuey().size() == 21276);
    CHECK(simulation1.getValueH().size() == 21276);
  }
}

TEST_CASE("Second set of parameters an initial value") {  // foglio 5
  pf::Simulation simulation2(1.1, 0.4, 0.1, 0.4, 10, 5, 0.001);

  simulation2.startingVectors();

  SUBCASE("testing values of equilibrium") {
    CHECK(simulation2.e2_x() == 2.75);
    CHECK(simulation2.e2_y() == 4);
  }
  SUBCASE("testing initial values of x,y,H") {
    CHECK(simulation2.getValuex()[0] == 10);
    CHECK(simulation2.getValuey()[0] == 5);
    CHECK(simulation2.getValueH()[0] == doctest::Approx(-23.91416));

    CHECK(simulation2.getValuex().size() == 1);
    CHECK(simulation2.getValuey().size() == 1);
    CHECK(simulation2.getValueH().size() == 1);
  }
  SUBCASE("testing valus after one step") {
    simulation2.nextstep(1);
    CHECK(simulation2.getValuex()[1] == doctest::Approx(9.991));
    CHECK(simulation2.getValuey()[1] == doctest::Approx(5.003));
    CHECK(simulation2.getValueH()[1] == doctest::Approx(-23.91685));

    CHECK(simulation2.getValuex().size() == 2);
    CHECK(simulation2.getValuey().size() == 2);
    CHECK(simulation2.getValueH().size() == 2);
  }

  SUBCASE("testing values after 0.004") {
    simulation2.nextstep(4);
    CHECK(simulation2.getValuex()[4] == doctest::Approx(9.963977));
    CHECK(simulation2.getValuey()[4] == doctest::Approx(5.011984));
    CHECK(simulation2.getValueH()[4] == doctest::Approx(-23.92486));

    CHECK(simulation2.getValuex().size() == 5);
    CHECK(simulation2.getValuey().size() == 5);
    CHECK(simulation2.getValueH().size() == 5);
  }

  SUBCASE("testing dimension of vectors after 0.020") {
    simulation2.nextstep(20);
    CHECK(simulation2.getValuex()[20] == doctest::Approx(9.8193));
    CHECK(simulation2.getValuey()[20] == doctest::Approx(5.05948));
    CHECK(simulation2.getValueH()[20] == doctest::Approx(-23.966));

    CHECK(simulation2.getValuex().size() == 21);
    CHECK(simulation2.getValuey().size() == 21);
    CHECK(simulation2.getValueH().size() == 21);
  }


  SUBCASE("testing values after 1") {
    simulation2.nextstep(1000);
    CHECK(simulation2.getValuex()[1000] == doctest::Approx(2.8785));
    CHECK(simulation2.getValuey()[1000] == doctest::Approx(6.06861));
    CHECK(simulation2.getValueH()[1000] == doctest::Approx(-21.3482));

    CHECK(simulation2.getValuex().size() == 1001);
    CHECK(simulation2.getValuey().size() == 1001);
    CHECK(simulation2.getValueH().size() == 1001);
  }
}

TEST_CASE(
    "testing the case in which preys goes extinct and after a long time "
    "predators too") {  // foglio 3
  pf::Simulation simulation3(0.4, 2.6, 0.3, 0.9, 4, 7, 0.001);

  simulation3.startingVectors();

  simulation3.nextstep(18000);

  CHECK(simulation3.getValuex()[18000] == doctest::Approx(0.0).epsilon(0.001));
  CHECK(simulation3.getValuey()[18000] == doctest::Approx(0.0).epsilon(0.001));
}

TEST_CASE(
    "testing the case in which predators goes extinct and prey increse "
    "indefinitetely") {
  pf::Simulation simulation4(1.8, 0.1, 0.001, 4, 10, 3, 0.001);

  simulation4.startingVectors();

  CHECK(simulation4.getValuey()[3000] == doctest::Approx(0).epsilon(0.001));
}

TEST_CASE("Testing equilibrium") {
  pf::Simulation simulation5(2.0, 0.01, 0.02, 1.5, 0, 0, 0.001);
  simulation5.startingVectors();
  simulation5.nextstep(57000);
  CHECK(simulation5.getValuex()[568911] == doctest::Approx(0).epsilon(0.001));
  CHECK(simulation5.getValuey()[6530] == doctest::Approx(0).epsilon(0.001));
  CHECK(simulation5.getValuex()[1] == doctest::Approx(0).epsilon(0.001));
  CHECK(simulation5.getValuey()[3556] == doctest::Approx(0).epsilon(0.001));
  CHECK(simulation5.getValuex()[12368] == doctest::Approx(0).epsilon(0.001));
  CHECK(simulation5.getValuey()[55278] == doctest::Approx(0).epsilon(0.001));
}
