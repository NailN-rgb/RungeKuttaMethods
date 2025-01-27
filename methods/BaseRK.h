#pragma once

#include <RungeKuttaMethods/methods/utils/vector_operations.h>

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

    using vector_of_values  = std::vector<value_type>;
    using vector_of_vectors = std::vector<vector_of_values>; 
    using equation_type     = RKequation::equation<index_type, value_type>;

protected:
    equation_type m_equation;
    vector_of_vectors m_calc_solution;
    vector_of_vectors m_expl_solution;
    vector_of_values m_t_points;
    value_type m_step;
    index_type m_steps_limit = 1e7;

public:
    explicit BaseRK(equation_type eq) 
    : m_equation{eq},
      m_step{m_equation.get_step()} 
    {
        if (m_step <= 0) 
        {
            throw std::invalid_argument("Step size must be positive");
        }
    };
    
    BaseRK(const BaseRK& RK)          = default;
    ~BaseRK()                         = default;

public:
    bool solve();

protected:
    bool set_initial_value(bool is_real_problem);

protected:
    virtual vector_of_values get_next_solution(
        value_type old_point,
        vector_of_values old_solution
    ) = 0;

protected:
    value_type get_next_step_size() {return m_step;}

protected:
    bool get_max_error();

protected:
    bool visualize();  

protected:
    bool visualize_trajectory();
};


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::solve()
{
    try
    {
        if (!set_initial_value(true)) {
            throw std::runtime_error("Failed to set initial value.");
        }

        auto t_old   = m_equation.get_start_t();
        auto y_start = m_equation.get_start_y();
        auto t_max   = m_equation.get_last_t();

        for (std::size_t i = 1; i < m_steps_limit; i++)
        {
            if (m_calc_solution.empty()) 
            {
                throw std::runtime_error("No previous solution for calculating next solution.");
            }

            m_calc_solution.push_back(get_next_solution(t_old, m_calc_solution.back()));

            t_old += m_step;

            m_expl_solution.push_back(m_equation.expl_sol(t_old));
            m_t_points.push_back(t_old);

            if (t_old > t_max - m_step) 
            {
                break;
            }
        }

        if (!get_max_error()) 
        {
            throw std::runtime_error("Failed to compute max error.");
        }

        // if (!visualize()) 
        // {
        //     throw std::runtime_error("Failed to visualize results.");
        // }

        if (!visualize_trajectory()) 
        {
            throw std::runtime_error("Failed to visualize results.");
        }
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error("Error in solve(): " + std::string(e.what()));
    }
    
    return true;
}


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::set_initial_value(bool is_real_problem)
{
    try
    {
        if(is_real_problem)
        {
            m_t_points.push_back(m_equation.get_start_t());
            m_calc_solution.push_back(m_equation.get_start_y());
            m_expl_solution.push_back(m_equation.get_start_y());
        }
        else
        {
            auto start_point = m_equation.get_start_t();

            m_calc_solution.push_back(m_equation.expl_sol(start_point));
            m_expl_solution.push_back(m_equation.expl_sol(start_point));
            m_t_points.push_back(start_point);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error in set_initial_value(): " << e.what() << '\n';
        return false;
    }

    return true;
}


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::get_max_error()
{
    try
    {
        if (m_calc_solution.size() != m_expl_solution.size()) 
        {
            throw std::runtime_error("Calculated and explicit solution sizes differ.");
        }

        vector_of_vectors err(m_calc_solution.size());

        std::transform(
            m_calc_solution.begin(),
            m_calc_solution.end(),
            m_expl_solution.begin(),
            err.begin(),
            [](const vector_of_values& calc, const vector_of_values& expl)
            {
                if (calc.size() != expl.size()) {
                    throw std::runtime_error("Solution vector sizes differ.");
                }
                vector_of_values err_res(calc.size());
                std::transform(calc.begin(), calc.end(), expl.begin(), err_res.begin(),
                               [](value_type calculated_val, value_type explicit_val) {
                                   return std::fabs(calculated_val - explicit_val);
                               });
                return err_res;
            }
        );

        vector_of_values error_by_time_steps(err.size());
        std::transform(
            err.begin(),
            err.end(),
            error_by_time_steps.begin(),
            [](const vector_of_values& errors_at_point)
            {
                return std::accumulate(errors_at_point.begin(), errors_at_point.end(), 0.0);
            }
        );

        auto max_error = *std::max_element(error_by_time_steps.begin(), error_by_time_steps.end());
        std::cout << "Max error: " << max_error << std::endl;
        std::cout << "Max error divided by stepsize " << max_error / std::pow(m_step, 4) << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error in get_max_error(): " << e.what() << '\n';
        return false;
    }

    return true;
}


template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::visualize()
{
    namespace plt = matplotlibcpp;

    try
    {
        if (m_calc_solution.empty() || m_expl_solution.empty() || m_t_points.empty()) {
            throw std::runtime_error("No data to visualize.");
        }

        index_type equations_count = m_calc_solution[0].size();
        for (index_type i = 0; i < equations_count; i++)
        {
            vector_of_values calculated_vector_solution(m_calc_solution.size());
            vector_of_values explicit_vector_solution(m_expl_solution.size());
            vector_of_values points_by_t(m_t_points.size());

            std::transform(
                m_calc_solution.begin(),
                m_calc_solution.end(),
                calculated_vector_solution.begin(),
                [i](const vector_of_values& calculated) { return calculated[i]; }
            );

            std::transform(
                m_expl_solution.begin(),
                m_expl_solution.end(),
                explicit_vector_solution.begin(),
                [i](const vector_of_values& explicit_sol) { return explicit_sol[i]; }
            );

            plt::figure_size(1200, 780);
            plt::named_plot("calc", m_t_points, calculated_vector_solution);
            plt::named_plot("Expl", m_t_points, explicit_vector_solution);
            plt::title("Correlation of calc & expl solutions");
            plt::legend();

            std::string filename = "./basic" + std::to_string(i) + ".png";
            plt::save(filename);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error when visualize(): " << e.what() << '\n';
        return false;
    }

    return true;
}

template<
    typename IndexType,
    typename ValueType
> bool BaseRK<IndexType, ValueType>::visualize_trajectory()
{
    namespace plt = matplotlibcpp;

    try
    {
        if (m_calc_solution.empty() || m_expl_solution.empty() || m_t_points.empty()) {
            throw std::runtime_error("No data to visualize.");
        }

        vector_of_values calculated_x_vector_solution(m_calc_solution.size());
        vector_of_values calculated_y_vector_solution(m_expl_solution.size());
        // vector_of_values points_by_t(m_t_points.size());

        std::transform(
            m_calc_solution.begin(),
            m_calc_solution.end(),
            calculated_x_vector_solution.begin(),
            [](const vector_of_values& calculated) { return calculated[0]; }
        );

        std::transform(
            m_calc_solution.begin(),
            m_calc_solution.end(),
            calculated_y_vector_solution.begin(),
            [](const vector_of_values& calculated) { return calculated[1]; }
        );

        plt::figure_size(1200, 780);
        plt::named_plot("body trajectory", calculated_x_vector_solution, calculated_y_vector_solution);
        plt::plot({1.}, {0.}, {{"color", "red"}, {"marker", "o"}, {"linestyle", "--"}});
        plt::plot({0.}, {0.}, {{"color", "green"}, {"marker", "o"}, {"linestyle", "--"}});
        // plt::xlim(-0.1, 1.1);
        // plt::ylim(-0.6, 0.6);
        plt::title("Resticted three body system");
        plt::legend();
        plt::xlabel("X");
        plt::ylabel("Y");

        std::string filename = "./traj.png";
        plt::save(filename);
        
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error in visualize(): " << e.what() << '\n';
        return false;
    }

    return true;
}

}//