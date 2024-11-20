#pragma once


/* Solver with tableau
*   0 |   0   0   0   0
* 1/2 | 1/2   0   0   0
* 1/2 |   0 1/2   0   0
*   1 |   0   0   1   0
* ____|_________________
*     | 1/6 1/3 1/3  1/6
*/

namespace solvers::ClassicRungeKutta
{

template<
    typename BaseSolver
>
class ClassicRungeKutta : public BaseSolver
{
private:
    using value_type       = typename BaseSolver::value_type;
    using index_type       = typename BaseSolver::index_type;
    using equation_type    = typename BaseSolver::equation_type;
    using vector_of_values = typename BaseSolver::vector_of_values;

private:
    equation_type m_equation;
    value_type m_step;

public:
    explicit ClassicRungeKutta(equation_type eq) 
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
        vector_of_values next_sol;

        auto k1 = m_equation.f(old_point, old_solution);
        auto k2 = m_equation.f(old_point + 1./2. * m_step, old_solution + 1./2. * m_step * k1);
        auto k3 = m_equation.f(old_point + 1./2. * m_step, old_solution + 1./2. * m_step * k2);
        auto k4 = m_equation.f(old_point + m_step, old_solution + m_step * k3);
        
        return old_solution + m_step * (1./6. * k1 + 1./3. * k2 + 1./3. * k3 +  1./6. * k4); 
    }
};

} // 