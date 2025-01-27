#pragma once

namespace RKequation
{
// for 1 equation
template<
    typename IndexType,
    typename ValueType
>
class equation
{
private:
    using index_type       = IndexType;
    using value_type       = ValueType;
    using vector_of_values = std::vector<value_type>;


private:
    value_type start_t;
    vector_of_values start_y;
    value_type end_point;
    value_type step;
    value_type m_m = 0.012277471;
    value_type m_M = 1 - m_m;

    // for adaptibe methods
    value_type estimate_eps = 0.001;
    value_type min_step = 0.0001;
    value_type max_step = 0.1;
    
public:
    equation(
        value_type start_t, 
        vector_of_values start_y,  
        value_type end, 
        value_type step,
        value_type estimate_eps = 1.0e-5,
        value_type min_step = 0.0001,
        value_type max_step = 0.1
    )
    : start_t(start_t),
    start_y(start_y),
    end_point(end),
    step(step),
    estimate_eps(estimate_eps),
    min_step(min_step),
    max_step(max_step)
    {}
    equation() = default;
    equation(const equation& eq) = default;
    ~equation() = default;

private:
    value_type r1(value_type x, value_type y)
    {
        return std::pow(std::pow(x + m_m, 2) + std::pow(y, 2), 1.5);
    }

private: 
    value_type r2(value_type x, value_type y)
    {
        return std::pow(std::pow(x - m_M, 2) + std::pow(y, 2), 1.5);
    }

public:
    vector_of_values f(value_type t, const vector_of_values& points) // t, {x1, y1, x2, y2}
    {   
        auto x = points[0], y = points[1], u = points[2], v = points[3];
        return vector_of_values{
            u,
            v,
            x + 2 * v - m_M * (x + m_m) / r1(x, y) - m_m * (x - m_M) / r2(x, y),
            y - 2 * u - m_M *  y        / r1(x, y) - m_m * y         / r2(x, y) 
        };
    }

public:
    vector_of_values expl_sol(value_type t)
    {
        vector_of_values res = {
            0,
            0,
            0,
            0
        };

        return res;
    }

public:
    value_type get_start_t()       { return start_t; }
    vector_of_values get_start_y() { return start_y; }
    value_type get_last_t()        { return end_point; }
    value_type get_step()          { return step; }
    // for adaptibe methods
    value_type get_epsilon()       { return estimate_eps; }
    value_type get_min_step()      { return min_step; }
    value_type get_max_step()      { return max_step; }
};

} //


