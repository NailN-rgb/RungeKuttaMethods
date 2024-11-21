#pragma once

#include<RungeKuttaMethods/methods/utils/vector_operations.h>


/* Solver with tableau
*  0 |   0    0
*  1 | 1/2  1/2
* ___|_________
*    | 1/2  1/2
*/

namespace solvers::CrankNicolson
{

template<
    typename BaseSolver
>
class CrankNicolson : public BaseSolver
{
private:
    using value_type       = typename BaseSolver::value_type;
    using index_type       = typename BaseSolver::index_type;
    using equation_type    = typename BaseSolver::equation_type;
    using vector_of_values = typename BaseSolver::vector_of_values;

private:
    equation_type m_equation;
    value_type m_step;

    index_type m_max_iter = 10;

public:
    explicit CrankNicolson(equation_type eq) 
    : BaseSolver(eq),
      m_equation{eq},
      m_step{eq.get_step()}
    {}

public:
    vector_of_values get_next_solution(
        value_type old_point,
        vector_of_values old_solution
    )
    {
        vector_of_values k_start(old_solution.size());
        // initial approximation
        std::fill(k_start.begin(), k_start.end(), 0);

        auto k1 = m_equation.f(old_point, old_solution);
        auto k2 = this->simple_iteration_nicolson(old_point, old_solution, k_start, k1);


        return old_solution + 1./2. * m_step * (k1 + k2); 
    }

private:
    auto simple_iteration_nicolson(
        value_type t,
        const vector_of_values& old_sol,
        const vector_of_values& k_start,
        const vector_of_values& k_1
    ) -> vector_of_values
    {
        try
        {
            vector_of_values k_new(old_sol.size());
            vector_of_values k_old(old_sol.size());

            k_old = k_start;

            for(auto i = 0; i < m_max_iter; i++)
            {
                k_new = m_equation.f(t + m_step, old_sol + 1./2. * m_step * k_1 + 1./2. * m_step * k_old);

                if(norm(k_old, k_new) < std::pow(m_step, 2)) { break;}
                else                                         { k_old = k_new; }

                if(i == m_max_iter - 1)                 { std::cout << "Iteration Limit" << std::endl; }
            }

            return k_new;
        }
        catch(const std::exception& e)
        {
            throw std::runtime_error("Simple iteration error: " +  std::string(e.what()));
        }
    }
};

} // 