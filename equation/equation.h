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
    value_type start_x;
    value_type start_y;
    value_type end_point;
    value_type step;
    
public:
    equation(
        value_type start_t, 
        value_type start_x, 
        value_type start_y,  
        value_type end, 
        value_type step
    )
    : start_t(start_t),
    start_x(start_x),
    start_y(start_y),
    end_point(end),
    step(step)
    {}
    equation() = default;
    equation(const equation& eq) = default;
    ~equation() = default;

public:
    vector_of_values f(vector_of_values points) // t,x,y
    {    
        return vector_of_values{
            2 * points[1] + 3 * points[2],
            -   points[1] +     points[2]
        };
    }

public:
    vector_of_values expl_sol(value_type t)
    {
        return vector_of_values{
            std::exp(2 * t) * std::cos(t),
            std::exp(2 * t) * std::sin(t)
        };
    }

public:
    vector_of_values get_start_point() { return vector_of_values{start_t, start_x, start_y}; }
    value_type get_end_point()   { return end_point; }
    value_type get_step()        { return step; }
};

} //


