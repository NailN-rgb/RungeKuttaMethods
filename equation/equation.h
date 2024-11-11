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
    using index_type = IndexType;
    using value_type = ValueType;

private:
    value_type start_point;
    value_type end_point;
    value_type step;
    
public:
    equation(value_type start, value_type end, value_type step)
    : start_point(start),
    end_point(end),
    step(step)
    {}
    equation() = default;
    equation(const equation& eq) = default;
    ~equation() = default;

public:
    value_type f(value_type x)
    {
        return std::sin(x);
    }

public:
    value_type expl_sol(value_type x)
    {
        return -std::cos(x) + 1;
    }

public:
    value_type get_start_point() {return start_point;}
    value_type get_end_point() {return end_point;}
    value_type get_step() {return step;}
};

} //


