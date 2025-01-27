#pragma once


/* Solver with tableau
*   1 |   0   0   0   0   0
* 1/3 | 1/3   0   0   0   0
* 1/3 | 1/3 1/6   0   0   0
* 1/2 | 1/8   0 3/8   0   0
*   1 | 1/2   0-3/2   2   0
* ____|____________________
*     | 1/6   0 4/6   0 1/6
*/

namespace solvers::MersonMethod
{

template<
    typename BaseSolver
>
class MersonMethod : public BaseSolver
{
private:
    using value_type       = typename BaseSolver::value_type;
    using index_type       = typename BaseSolver::index_type;
    using equation_type    = typename BaseSolver::equation_type;
    using vector_of_values = typename BaseSolver::vector_of_values;

private:
    equation_type m_equation;
    value_type m_step;
    value_type m_estimate_epsilon;
    value_type m_max_step;
    value_type m_min_step;

public:
    explicit MersonMethod(equation_type eq) 
    : BaseSolver(eq),
      m_equation{eq},
      m_step{eq.get_step()},
      m_estimate_epsilon{eq.get_epsilon()},
      m_min_step{eq.get_min_step()},
      m_max_step{eq.get_max_step()}
    {}

public:
    vector_of_values get_next_solution(
        value_type old_point,
        vector_of_values old_solution
    )
    {
        vector_of_values next_sol;
        
        auto k1 = m_equation.f(old_point, old_solution);
        auto k2 = m_equation.f(old_point + 1./3. * m_step, old_solution + 1./3. * m_step * k1);
        auto k3 = m_equation.f(old_point + 1./3. * m_step, old_solution + 1./6. * m_step * k1 + 1./6. * m_step * k2);
        auto k4 = m_equation.f(old_point + 1./2. * m_step, old_solution + 1./8. * m_step * k1 + 3./8. * m_step * k3);
        auto k5 = m_equation.f(old_point + m_step, old_solution + 1./2. * m_step * k1 - 
                                                    3./2. * m_step * k3 + 2 * m_step * k4);

        // auto E = 1./30. * (2 * k1 - 9 * k3 + 8 * k4 - k5);

        next_sol = 1./6. * m_step * (k1 + 4 * k4 + k5);

        return old_solution + next_sol;
    }
};

} // 