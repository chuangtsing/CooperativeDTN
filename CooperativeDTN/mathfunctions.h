#pragma once

#include <random>

int factorial(int);


template<typename type=double>
class pareto_distribution : public exponential_distribution<type>
{
public:

	pareto_distribution(type l_scale=1, type l_shape = 1) : exponential_distribution(l_shape), m_scale(l_scale), m_shape(l_shape)
    {
    }

	~pareto_distribution()
	{
	}

	type scale() const
	{
		return m_scale;
	}

	type shape() const
	{
		return m_shape;
	}

	type operator()(default_random_engine& g) 
	{
		type tmp = exponential_distribution::operator()(g);
		return exp(tmp)*m_scale;
	}

private:
	type m_shape;
	type m_scale;
};