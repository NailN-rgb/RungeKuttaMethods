#pragma once


/* Solver with tableau
*  1 |   1
* ___|________
*    |   1
*/

namespace solvers::BackwardEuler
{

template<
    typename BaseSolver
>
class BackwardEuler : public BaseSolver
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
    explicit BackwardEuler(equation_type eq) 
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

        vector_of_values k1 = this->simple_iteration(old_point, old_solution, k_start);
        
        return old_solution + m_step * k1; 
    }

private:
    auto simple_iteration(
        value_type t,
        const vector_of_values& old_sol,
        const vector_of_values& k_start
    ) -> vector_of_values
    {
        try
        {
            vector_of_values k_new(old_sol.size());
            vector_of_values k_old(old_sol.size());

            k_old = k_start;

            for(auto i = 0; i < m_max_iter; i++)
            {
                k_new = m_equation.f(t + m_step, old_sol + m_step * k_old);

                if(norm(k_old, k_new) < std::pow(m_step, 2)) { break;}
                else                                         { k_old = k_new; }

                if(i == m_max_iter - 1)                      { std::cout << "Iteration Limit" << std::endl; }
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