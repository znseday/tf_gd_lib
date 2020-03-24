//          Copyright Sergey Tsynikin 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE tf_gd_lib_test_module

#include "UnitSpline.h"
#include "UnitTableFunctions.h"
#include "UnitGradDescent.h"

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <tuple>

#include <iostream>


using namespace std;
using namespace tf_gd_lib;

bool CmpFunc(double a, double b, double eps)
{
	return (fabs(a - b) < eps);
}

BOOST_AUTO_TEST_SUITE(tf_gd_lib_test_suite)

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_boost_test)
{
    BOOST_CHECK(1 > 0);
}

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_tabulated_function_test)
{
	TableFunction tf;
	tf.CreateDemoFunction(401, -10, 0.05);

	cout << tf(-10.03) << endl;
	BOOST_CHECK( CmpFunc(tf(-10.03),  0.22891692244, 0.0000001) ); // left extrapolation
	BOOST_CHECK( CmpFunc(tf( 10.04), -0.13059085494, 0.0000001) ); // right extrapolation
	BOOST_CHECK( CmpFunc(tf(  4.01),  0.67480799289, 0.0000001) ); // interpolation
	BOOST_CHECK( CmpFunc(tf(0), 0, 0.00000000001) );               // center point
	
	tf.ClearAll();
	tf.CreateNewFunction(5, "Four_points_function");

	BOOST_CHECK( tf.Size() == 5 );

	// create some unordered data with duplicates
	tf.SetPoint(0, -1, 0.5);
	tf.SetPointByNumber(1, make_tuple(1, -0.125));
	tf.SetPoint(2,  0, 1.5);
	tf.SetPoint(3, -1, 0.75); // it's duplicate of Point[0] by x
	tf.SetPointByNumber(4, make_tuple(2, 0.25));

	tf.Sort();
	tf.KillDuplicates();
	tf.CalcStat();

	cout << "size = " << tf.Size();
	BOOST_CHECK( tf.Size() == 4 );

	for (size_t i = 0; i < tf.Size(); ++i)
	{
		auto point = tf[i];
		cout << get<0>(point) << "\t" << get<1>(point) << endl;
	}

	BOOST_CHECK( tf.GetName() == "Four_points_function" );

	// check stats
	BOOST_CHECK( tf.GetMinX() == -1 );
	BOOST_CHECK( tf.GetMaxX() == 2 );
	BOOST_CHECK( tf.GetMinY() == -0.125 );
	BOOST_CHECK( tf.GetMaxY() == 1.5);

	BOOST_CHECK( tf.Get_x_ForMinY() == 1);
	BOOST_CHECK( tf.Get_x_ForMaxY() == 0);
	BOOST_CHECK( tf.Get_i_ForMinY() == 1);
	BOOST_CHECK( tf.Get_i_ForMaxY() == 2);

}


BOOST_AUTO_TEST_CASE(tf_gd_lib_test_tf_load_and_spline_test)
{
	TableFunction tf;
	BOOST_CHECK( tf.LoadFromFile("test_data.txt") );

	BOOST_CHECK( tf.BuildSpline() );

	const double e = 1e-3;
	BOOST_CHECK( CmpFunc(tf.Spline(-0.185), 1.930, e) ); // at random point
	BOOST_CHECK( CmpFunc(tf.Spline(0.707), -0.335, e) ); // at random point
	BOOST_CHECK( CmpFunc(tf.Spline(0.5), -0.250, e) );   // at existing point
}

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_gradient_descent_test)
{
	GradDescent gd;
	BOOST_CHECK( 1 > 0);

}


BOOST_AUTO_TEST_SUITE_END()