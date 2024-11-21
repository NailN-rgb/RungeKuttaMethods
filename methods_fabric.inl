#pragma once
#include <RungeKuttaMethods/methods/BaseRK.h>

// RK2
#include <RungeKuttaMethods/methods/RK2/HeunsMethod.h>
#include <RungeKuttaMethods/methods/RK2/ButcherMethod.h>

//RK3
#include <RungeKuttaMethods/methods/RK3/KuttaMethod.h>
#include <RungeKuttaMethods/methods/RK3/Heun3dMethod.h>

//RK4
#include <RungeKuttaMethods/methods/RK4/ClassicMethod.h>

// implicit
#include <RungeKuttaMethods/methods/ImplicitMethods/bakward_euler.h>
#include <RungeKuttaMethods/methods/ImplicitMethods/crank_nicolcon.h>


namespace Solvers_Fabric
{

struct Solvers_Fabric
{
    template<
        typename index_type,
        typename value_type,
        typename equation_type
    > static std::unique_ptr<solvers::BaseRK<index_type, value_type>> apply(std::string solver_id, equation_type equation) 
    {
        using base_solver_type      = solvers::BaseRK<index_type, value_type>;
        //RK2
        using Heuns_solver_type     = solvers::HeunsMethod::HeunsMethod<base_solver_type>;
        using Butcher_solver_type   = solvers::ButcherMethod::ButcherMethod<base_solver_type>;
        //RK3 
        using Kutta_solver_type     = solvers::KuttaMethod::KuttaMethod<base_solver_type>;
        using Heun3d_solver_type    = solvers::Heun3dMethod::Heun3dMethod<base_solver_type>;
        //RK4
        using ClassicRK_solver_type = solvers::ClassicRungeKutta::ClassicRungeKutta<base_solver_type>;

        // implicit methods
        using backward_euler_solver_type = solvers::BackwardEuler::BackwardEuler<base_solver_type>;
        using crank_nicolson_solver_type = solvers::CrankNicolson::CrankNicolson<base_solver_type>;


        if(solver_id == "RK2Heuns")
        {
            return std::make_unique<Heuns_solver_type>(equation);
        }
        else if(solver_id == "RK2Butcher")
        {
            return std::make_unique<Butcher_solver_type>(equation);
        }
        else if(solver_id == "RK3Kutta")
        {
            return std::make_unique<Kutta_solver_type>(equation);
        }
        else if(solver_id == "RK3Heuns")
        {
            return std::make_unique<Heun3d_solver_type>(equation);
        }
        else if(solver_id == "RK4Classic")
        {
            return std::make_unique<ClassicRK_solver_type>(equation);
        }
        else if(solver_id == "BacwardEuler")
        {
            return std::make_unique<backward_euler_solver_type>(equation);
        }
        else if(solver_id == "BackwardCrankNicolson")
        {
            return std::make_unique<crank_nicolson_solver_type>(equation);
        }
        else
        {
            throw std::runtime_error("unknown type of solver");
        }
    }
};

}