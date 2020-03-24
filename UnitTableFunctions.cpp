﻿//          Copyright Sergey Tsynikin 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <tuple>
#include <fstream>

#include "UnitTableFunctions.h"

using namespace std;
using namespace tf_gd_lib;

inline double tf_gd_lib::LineInterpol(double x, double x1, double x2, double y1, double y2)
{
	return y1 + (x-x1)*(y2-y1)/(x2-x1);
}
//---------------------------------------------------------------------------

inline double tf_gd_lib::LineInterpolSafeMiddleVal(double x, double x1, double x2, double y1, double y2)
{
	double dx = x2 - x1;
	if (dx == 0)
	{
		if (x < x1)
			return y1;
		else if (x > x2)
			return y2;
		else
			return (y1+y2)/2.0;
	}
	else
		return y1 + (x-x1)*(y2-y1)/(x2-x1);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

std::tuple<double,double> TableFunction::operator[](size_t i)
{
	iCash = i;
	// return a copy
	// could we return a reference?
	return make_tuple(Points[i].x, Points[i].y);
}
//---------------------------------------------------------------------------

double TableFunction::GetValByBSearchFromX(double x) const
{
	if (Points.size() < 2)
	{
		return 0;
	}

	auto it = lower_bound(Points.begin(), Points.end(), SinglePoint(x, 0), [](const SinglePoint &a, const SinglePoint &b)
		{
			return a.x < b.x;
		});

	if (it == Points.begin())
	{
		iCash = 0;

		return LineInterpol(x, it->x, (it+1)->x, it->y, (it+1)->y);
	}
	else if (it == Points.end())
	{
		iCash = Points.size() - 1;

		return LineInterpol(x, (it-2)->x, (it-1)->x, (it-2)->y, (it-1)->y);
	}
	else
	{
		iCash = distance(Points.begin(), it) - 1;

		return LineInterpol(x, (it-1)->x, it->x, (it-1)->y, it->y);
	}

}
//---------------------------------------------------------------------------

double TableFunction::GetValFromRightX(double x) const
{
	if (Points.size() <2)
	{
		return 0;
	}

	if (x < Points[iCash].x)
	{
		return Points[iCash].y;  // Or, as an option, to do left extrapolation
	}

	for (size_t i = iCash; i < Points.size(); ++i)
	{
		if (x < Points[i].x)   // Interpolation
		{
			iCash = i-1;
			return LineInterpol(x, Points[i-1].x, Points[i].x, Points[i-1].y, Points[i].y);
		}
	}

	// Right Extrapolation 
	iCash = Points.size()-1;
	return LineInterpol(x, Points[Points.size()-2].x, Points[Points.size()-1].x,
						   Points[Points.size()-2].y, Points[Points.size()-1].y);

}
//---------------------------------------------------------------------------

double TableFunction::GetValFromLeftX(double x) const
{
	if (Points.size() <2)
	{
		return 0;
	}

	if (x < Points[0].x)       // Left Extrapolation
	{
		iCash = 0;
		return LineInterpol(x, Points[0].x, Points[1].x, Points[0].y, Points[1].y);
	}

	int n = min(iCash+2, Points.size());

	for (size_t i = 1; i < n; ++i)
	{
		if (x < Points[i].x)   // Interpolation
		{
			iCash = i-1;
			return LineInterpol(x, Points[i-1].x, Points[i].x, Points[i-1].y, Points[i].y);
		}
	}

	return Points[iCash].y;   // Or, as an option, to do right extrapolation
}
//---------------------------------------------------------------------------

void TableFunction::ClearAll()
{
	Points.clear();
    Points.shrink_to_fit();

	MinX = MaxX = 0;
	MinY = MaxY = 0;
	x_ForMinY = x_ForMaxY = 0;
	i_ForMinY = i_ForMaxY = 0;

    iCash = 0;

	Name.clear();

    Spline.Clear();
}
//---------------------------------------------------------------------------

void TableFunction::SetValAtPoint(size_t i, double y)
{
	Points[i].y = y;
}
//---------------------------------------------------------------------------

void TableFunction::SetPointByNumber(size_t i, const std::tuple<double,double> &point)
{
	Points[i].x = get<0>(point);
	Points[i].y = get<1>(point);
}
//---------------------------------------------------------------------------

void TableFunction::CreateNewFunction(size_t n, const string &_name)
{
	ClearAll();
	Points.reserve(n);
	for (size_t i = 0; i < n; ++i)
		Points.emplace_back(0.0, 0.0);

    //Points.shrink_to_fit();
	Name = _name;
}
//---------------------------------------------------------------------------

void TableFunction::CreateDemoFunction(size_t n, double a, double dx, std::function<double(double)> f, const std::string &_name)
{
	ClearAll();
	Points.reserve(n);
	for (size_t i = 0; i < n; ++i)
	{
		double x = a + i*dx;
		Points.emplace_back(x, f(x));
	}
	//Points.shrink_to_fit();

    CalcStat();

    Name = _name;
}
//---------------------------------------------------------------------------

void TableFunction::KillDuplicates()
{
	unique(Points.begin(), Points.end(), [](const SinglePoint &a, const SinglePoint &b)
		{
			return a.x == b.x;
		});
}
//---------------------------------------------------------------------------

void TableFunction::Sort()
{                                  
	sort(Points.begin(), Points.end(), [](const SinglePoint &a, const SinglePoint &b)
		{
			return a.x < b.x;
		});
}
//---------------------------------------------------------------------------

void TableFunction::CalcStat()
{
	if (Points.empty())
		return;

	MinX = MaxX = Points[0].x;
	MinY = MaxY = Points[0].y;

	x_ForMinY = x_ForMaxY = Points[0].x;

	for (size_t i = 1; i < Points.size(); ++i)
	{
		if (Points[i].x > MaxX)
			MaxX = Points[i].x;

		if (Points[i].x < MinX)
			MinX = Points[i].x;

		if (Points[i].y > MaxY)
		{
			MaxY = Points[i].y;
			x_ForMaxY = Points[i].x;
			i_ForMaxY = i;
		}

		if (Points[i].y < MinY)
		{
			MinY = Points[i].y;
			x_ForMinY = Points[i].x;
			i_ForMinY = i;
		}
	}
}
//---------------------------------------------------------------------------

bool TableFunction::LoadFromFile(const string &FileName) // vs. wstring
{
	ifstream f(FileName); // vs. wifstream
	if (!f)
        return false;

    LoadFromStream(f);  // TO DO: test on travis

	Name = FileName;
	size_t pos = Name.find_last_of("\\");  // TO DO: Check separator in different OS

	if (pos != string::npos)
	{
		Name.erase(0, pos+1);
	}

	return true;
}
//---------------------------------------------------------------------------

void TableFunction::LoadFromStream(istream &Stream) // vs. wistream
{
	ClearAll();

	string line;
	while (getline(Stream, line))
	{
		double x, y;
		sscanf(line.c_str(), "%lf%lf", &x, &y);
		Points.emplace_back(x, y);
	}
	Points.shrink_to_fit();

	Sort();
	CalcStat();
}
//---------------------------------------------------------------------------

bool TableFunction::BuildSpline()
{
    return Spline.BuildSpline(Points);
}
//---------------------------------------------------------------------------