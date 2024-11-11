#pragma once


namespace solvers::RK2
{

template<
    typename BaseSolver
>
class RK2 : public BaseSolver
{
private:
    using value_type    = typename BaseSolver::value_type;
    using index_type    = typename BaseSolver::index_type;
    using equation_type = typename BaseSolver::equation_type;

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
    value_type get_next_solution(
        value_type old_point,
        value_type old_solution
    )
    {
        auto k1 = m_equation.f(old_point);
        auto k2 = m_equation.f(old_point + 2/3 * m_step);

        return old_solution + m_step * (0.25 * k1 + 0.75 * k2);
    }
};

} // 