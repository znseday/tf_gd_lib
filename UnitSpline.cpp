//          Copyright Sergey Tsynikin 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

//#include <cassert>

#include "UnitSpline.h"

using namespace std;
using namespace tf_gd_lib;

bool CubicSpline::BuildSpline(const std::vector<SinglePoint> &Points)
{
    size_t n = Points.size();

	if (n < 3) // Или просто if и игнор (return), если значения плохие?  или вообще возвращать bool?????
		return false;

	Splines.clear();
	//Splines.reserve(n); // is it necessary before resize?
	Splines.resize(n);
	//Splines.shrink_to_fit(); // is it necessary after resize?

	for (size_t i = 0; i < n; ++i)
	{
		Splines[i].x = Points[i].x;
		Splines[i].a = Points[i].y;
	}
	Splines[0].c = 0.0;

	vector<double> alpha(n-1);
	vector<double> beta(n-1);

	double A, B, C, F, h_i, h_i1, z;
	alpha[0] = beta[0] = 0.0;

	for (size_t i = 1; i < n-1; ++i)
	{
		h_i = Points[i].x - Points[i-1].x, h_i1 = Points[i+1].x - Points[i].x;
		A = h_i;
		C = 2.0 * (h_i + h_i1);
		B = h_i1;
		F = 6.0 * ((Points[i+1].y - Points[i].y) / h_i1 - (Points[i].y - Points[i-1].y) / h_i);
		z = (A * alpha[i-1] + C);
		alpha[i] = -B / z;
		beta[i] = (F - A * beta[i-1]) / z;
	}

	Splines[n-1].c = (F - A * beta[n-2]) / (C + A * alpha[n-2]);

	for (long long i = n - 2; i > 0; --i)
		Splines[i].c = alpha[i] * Splines[i+1].c + beta[i];

	for (long long i = n - 1; i > 0; --i)
	{
		h_i = Points[i].x - Points[i-1].x;
		Splines[i].d = (Splines[i].c - Splines[i-1].c) / h_i;
		Splines[i].b = h_i * (2.0 * Splines[i].c + Splines[i-1].c) / 6.0 + (Points[i].y - Points[i-1].y) / h_i;
	}
}
//---------------------------------------------------------------------------

double CubicSpline::operator()(double x) const
{
	if (Splines.empty())  // If splines don't exist - return NaN
		return std::numeric_limits<double>::quiet_NaN();

	SplinePart s;
	if (x <= Splines[0].x)
		s = Splines[1];
	else if (x >= Splines.back().x) 
		s = Splines.back();
	else 
	{
		size_t i = 0, j = Splines.size() - 1;
		while (i + 1 < j)
		{
			size_t k = i + (j - i) / 2;
			if (x <= Splines[k].x)
				j = k;
			else
				i = k;
		}
		s = Splines[j];
	}

	double dx = (x - s.x);
	return s.a + (s.b + (s.c / 2. + s.d * dx / 6.0) * dx) * dx;
}
//---------------------------------------------------------------------------
