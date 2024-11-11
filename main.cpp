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
    using value_type = double;
    using index_type = int;

    using equation_type = RKequation::equation<index_type, value_type>;

    // enter datas for solvers start
    value_type start = 0.;
    value_type end = 5.;
    value_type step = 0.1;

    // enter solver type
    std::string solver_type = "RK2";

    // init solver start datas
    equation_type eq(start, end, step);


    // call solvers fabric
    auto solver = Solvers_Fabric::Solvers_Fabric::apply<
        index_type, value_type, equation_type
    >(solver_type, eq);

    solver->solve();

    // if print this -> OK
    std::cout << "Programm Finished" << std::endl;

    return 0;
}
