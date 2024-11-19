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
    
public:
    equation(
        value_type start_t, 
        vector_of_values start_y,  
        value_type end, 
        value_type step
    )
    : start_t(start_t),
    start_y(start_y),
    end_point(end),
    step(step)
    {}
    equation() = default;
    equation(const equation& eq) = default;
    ~equation() = default;

public:
    vector_of_values f(value_type t, const vector_of_values& points) // t,x,y
    {    
        return vector_of_values{
            - 2. * points[0] + 4. * points[1],
            -      points[0] + 3. * points[1]
        };
    }

public:
    vector_of_values expl_sol(value_type t)
    {
        vector_of_values res = {
            4. * std::exp(-t) - std::exp(2. * t),
                 std::exp(-t) - std::exp(2. * t)
        };

        return res;
    }

public:
    value_type get_start_t()       { return start_t; }
    vector_of_values get_start_y() { return start_y; }
    value_type get_last_t()        { return end_point; }
    value_type get_step()          { return step; }
};

} //


