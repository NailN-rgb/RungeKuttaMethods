#pragma once
#include <RungeKuttaMethods/methods/BaseRK.h>
#include <RungeKuttaMethods/methods/RK2/RK2.h>

namespace Solvers_Fabric
{

struct Solvers_Fabric
{
    template<
        typename index_type,
        typename value_type,
        typename equation_type
    > static auto apply(std::string solver_id, equation_type equation)
    {
        using base_solver_type = solvers::BaseRK<index_type, value_type>;
        using rk2_solver_type  = solvers::RK2::RK2<base_solver_type>;

        if(solver_id == "RK2")
        {
            return std::make_unique<rk2_solver_type>(equation);
        }
        else
        {
            throw std::runtime_error("unknown type of solver");
        }
    }
};

}