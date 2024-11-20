#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <exception>
#include <memory>
#include <algorithm>
#include <functional>
#include <numeric>

#include <RungeKuttaMethods/equation/equation.h>
#include <RungeKuttaMethods/methods_fabric.inl>


int main()
{
    using value_type       = double;
    using index_type       = int;
    using vector_of_values = std::vector<value_type>;
    using equation_type    = RKequation::equation<index_type, value_type>;

    // enter datas for solvers start
    value_type start_t       = 0.;
    vector_of_values start_y = {0.0, 0.0};
    value_type end           = 1.;
    value_type step          = 0.01;

    // enter solver type
    // possible solvers
    // RK2Heuns
    // RK2Butcher
    // RK3Kutta
    // RK3Heuns
    // RK4Classic
    std::string solver_type = "RK4Classic";

    // init solver start datas
    equation_type eq(start_t, start_y, end, step);


    // call solvers fabric
    auto solver = Solvers_Fabric::Solvers_Fabric::apply<
        index_type, value_type, equation_type
    >(solver_type, eq);

    solver->solve();

    // if print this -> OK
    std::cout << "Programm Finished" << std::endl;

    return 0;
}
