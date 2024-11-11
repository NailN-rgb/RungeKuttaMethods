#pragma once

#include <RungeKuttaMethods/plot/matplotlibcpp.h>

namespace solvers
{

template<
    typename IndexType,
    typename ValueType
>
class BaseRK
{
public:
    using index_type = IndexType;
    using value_type = ValueType;

    using vector_of_values = std::vector<value_type>;
    using equation_type    = RKequation::equation<index_type, value_type>;

protected:
    equation_type m_equation;
    vector_of_values m_calc_solution;
    vector_of_values m_expl_solution;
    vector_of_values m_x_points;
    value_type m_step;
    index_type m_steps_limit = 10000;

public:
    explicit BaseRK(equation_type eq) 
    : m_equation{eq},
      m_step{m_equation.get_step()} 
      {};
    
    BaseRK(const BaseRK& RK)          = default;
    ~BaseRK()                         = default;

public:
    bool solve();

private:
    bool set_initial_value();

private:
    virtual value_type get_next_solution(
        value_type old_point,
        value_type old_solution
    ) = 0;

private:
    value_type get_next_step_size() {return m_step;}

private:
    bool get_max_error();

public:
    bool visualize();
};


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::solve()
{
    try
    {
        set_initial_value();

        auto x_old = m_equation.get_start_point();
        auto x_max = m_equation.get_end_point();

        for(std::size_t i = 1; i < m_steps_limit; i++)
        {
            // calculate solution
            m_calc_solution.push_back(get_next_solution(
                x_old,
                m_calc_solution.back() 
            ));

            x_old += m_step;

            // true solution
            m_expl_solution.push_back(
                m_equation.expl_sol(x_old)
            );

            m_x_points.push_back(x_old);

            if(x_old > x_max - m_step)
            {
                // logic for last point
                break;
            }
        }

        get_max_error();

        visualize();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    
    return true;
}


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::set_initial_value()
{
    m_calc_solution.push_back(
        m_equation.expl_sol(
            m_equation.get_start_point()
        )
    );

    m_expl_solution.push_back(
        m_equation.expl_sol(
            m_equation.get_start_point()
        )
    );

    m_x_points.push_back(m_equation.get_start_point());

    return true;
}


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::get_max_error()
{
    if(m_calc_solution.size() != m_expl_solution.size())
    {
        throw std::runtime_error("Calculated and explixit solution sizes are dirrerent\n");
    }

    vector_of_values err(m_calc_solution.size());

    std::transform(
        m_calc_solution.begin(),
        m_calc_solution.end(),
        m_expl_solution.begin(),
        err.begin(),
        [&](auto calc, auto expl)
        {
            return std::fabs(calc) - std::fabs(expl);
        }
    );

    std::sort(err.begin(), err.end());

    std::cout << "Max error: " << err.front() << std::endl;

    return true;
}


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::visualize()
{
    namespace plt = matplotlibcpp;

    plt::figure_size(1200, 780);

    plt::named_plot("calc", m_x_points, m_calc_solution);

    plt::named_plot("Expl", m_x_points, m_expl_solution);

    plt::title("Correlation of calc & expl solutions");

    plt::legend();

    const char* filename = "./basic.png";
    plt::save(filename);

    return true;
}

}//