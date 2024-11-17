#pragma once


namespace solvers::RK2
{

template<
    typename BaseSolver
>
class RK2 : public BaseSolver
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
    explicit RK2(equation_type eq) 
    : BaseSolver(eq),
      m_equation{eq},
      m_step{eq.get_step()}
    {}

public:
    vector_of_values get_next_solution(
        vector_of_values old_point,
        vector_of_values old_solution
    )
    {
        vector_of_values next_sol;

        auto k1 = m_equation.f(old_point);
        auto k2 = m_equation.f(old_point + (2/3 * m_step));

        return old_solution + m_step * (0.25 * k1 + 0.75 * k2);
    }
};

} // 