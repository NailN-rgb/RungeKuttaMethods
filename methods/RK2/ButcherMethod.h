#pragma once


/* Solver with tableau
*  0 |   0   0
* 2/3| 2/3   0
* ___|________
*    | 1/4 3/4
*/

namespace solvers::ButcherMethod
{

template<
    typename BaseSolver
>
class ButcherMethod : public BaseSolver
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
    explicit ButcherMethod(equation_type eq) 
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
        auto k2 = m_equation.f(old_point + 2./3. * m_step, old_solution + 2./3. * m_step * k1);
        
        return old_solution + m_step * (1./4. * k1 + 3./4. * k2); 
    }
};

} // 