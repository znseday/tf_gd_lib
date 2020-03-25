//          Copyright Sergey Tsynikin 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE tf_gd_lib_test_module

#include "UnitSpline.h"
#include "UnitTableFunctions.h"
#include "UnitGradDescent.h"

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <cmath>
#include <tuple>

#include <iostream>

#include <fstream>

using namespace std;
using namespace tf_gd_lib;

bool CmpFunc(double a, double b, double eps)
{
	return (fabs(a - b) < eps);
}
//---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(tf_gd_lib_test_suite)

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_boost_test)
{
	cerr << "------------------> TEST TEST TEST <------------------" << endl;
    BOOST_CHECK(1 > 0);
}
//---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_tabulated_function_test)
{
	TableFunction tf;
	tf.CreateDemoFunction(201, -1, 0.01);

	const double e = 0.005; // accuracy for linear interpolation/extrapolation testing

	cout << "tf(-1.003) = " << tf(-1.003) << endl;
	BOOST_CHECK( CmpFunc(tf(-1.003),  0.5689446, e) ); // left extrapolation

	cout << "tf(1.004) = " << tf(1.004) << endl;
	BOOST_CHECK( CmpFunc(tf (1.004), -0.5771398, e) ); // right extrapolation

	cout << "tf(0.117) = " << tf(0.117) << endl;
	BOOST_CHECK( CmpFunc(tf( 0.117),  0.9207505, e) ); // interpolation

	cout << "tf(0) = " << tf(0) << endl;
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

	
	cout << "size = " << tf.Size() << endl;
	BOOST_CHECK( tf.Size() == 4 ); // one duplicate was removed

	cout << "x\ty" << endl;
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
	BOOST_CHECK( tf.Get_i_ForMinY() == 2); // considering that sort was done
	BOOST_CHECK( tf.Get_i_ForMaxY() == 1); // considering that sort was done

}
//---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_tf_load_and_spline_test)
{
	TableFunction tf;
	BOOST_CHECK( tf.LoadFromFile("test_data.txt") );

	BOOST_CHECK( tf.BuildSpline() );

	const double e = 0.001; // 1e-3
	BOOST_CHECK( CmpFunc(tf.Spline(-0.185), 1.930, e) ); // at some random point
	BOOST_CHECK( CmpFunc(tf.Spline(0.707), -0.335, e) ); // at some random point
	BOOST_CHECK( CmpFunc(tf.Spline(0.5), -0.250, e) );   // at existing point

	BOOST_CHECK( tf.Spline(-0.185) != tf(-0.185) ); // at some random point
	BOOST_CHECK( tf.Spline(0.707)  != tf( 0.707) ); // at some random point

	BOOST_CHECK( tf.Spline(0.5) == tf(0.5) );   // at existing point
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double linear_experimental(double x)
{
	return 1.23 * x - 0.123;
}

double linear_predict(double x, const std::vector<double> &p)
{
	return p[0] * x + p[1];
}

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_gd_linear_test)
{
	GradDescent gd;	
	
	TableFunction experimental;   // Generate some linear experimental data
	experimental.CreateDemoFunction(101, -10, 0.25, linear_experimental, "Linear_experimental");

	gd.SetSrcFunction(experimental);
	gd.SetDstFunction(linear_predict);

	gd.SetAlpha(0.5);   // value for momentum
	gd.SetEps(0.00001); // value to derivative 
	gd.SetEta_FirstJump(10);  
	gd.SetEta_k_inc(1.08);
	gd.SetEta_k_dec(2.0);
	gd.SetMin_Eta(1e-8);       // min descent rate (will be multiplied by FirstJump before get started)
	gd.SetFinDifMethod(false); // using a plain central derivative
	gd.SetMaxIters(10000);     // iteration limit
	gd.SetMaxTime(3);          // time limit (seconds)

	const int param_count = 2;
	vector<double> params(param_count);
	vector<double> min_constrains(param_count);
	vector<double> max_constrains(param_count);
	vector<double> rel_constrains(param_count); // relative constrains in %
	vector<bool> type_constrains(param_count);  // set 'false' to use absolute constrains, 'true' - for relative

	params[0] = params[1] = 0; // start parameters
	min_constrains[0] = min_constrains[1] = -1000;
	max_constrains[0] = max_constrains[1] =  1000;
	rel_constrains[0] = rel_constrains[1] = 0;
	type_constrains[0] = type_constrains[1] = false; // use absolute constrains

	gd.SetParams(params);
	gd.SetMinConstrains(min_constrains);
	gd.SetMaxConstrains(max_constrains);
	gd.SetRelConstrains(rel_constrains);
	gd.SetTypeConstrains(type_constrains);

	GradErrorType res = gd.Go();

	cout << "gd_linear result: ";
	switch (res)
	{
	case GradErrorType::Success:
		cout << "Success" << endl;
		break;
	case GradErrorType::VectorSizesNotTheSame:
		cout << "VectorSizesNotTheSame" << endl;
		break;
	case GradErrorType::CanceledByUser:
		cout << "CanceledByUser" << endl;
		break;
	case GradErrorType::TimeOut:
		cout << "TimeOut" << endl;
		break;
	case GradErrorType::ItersOverflow:
		cout << "ItersOverflow" << endl;
		break;
	}

	cout << "gd.GetLastCost() = " << gd.GetLastCost() << endl;

	params = gd.GetParams();

	BOOST_CHECK(CmpFunc(gd.GetLastCost(), 0, 1e-16)); 

	cout << "Found Params: " << endl;
	for (auto p : params)
		cout << p << endl;

	BOOST_CHECK(CmpFunc(params[0],  1.23,  1e-12)); // almost the same as in SrcFunction
	BOOST_CHECK(CmpFunc(params[1], -0.123, 1e-12)); // almost the same as in SrcFunction
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double polynominal_experimental(double x)
{
	double noise = rand() / (double)RAND_MAX / 25.0; // add some noise 
	return 0.5*x*x*x - 1.1*x*x + 0.75*x - 1.5 + noise;
}

double polynominal_predict(double x, const std::vector<double>& p)
{
	return p[0] * x * x * x + p[1] * x * x + p[2] * x + p[3];	              
}

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_gd_polynominal_test)
{
	GradDescent gd;
	
	TableFunction experimental; // Generate some polynominal experimental data
	experimental.CreateDemoFunction(101, -10, 0.25, polynominal_experimental, "polynominal_experimental");

	gd.SetSrcFunction(experimental);
	gd.SetDstFunction(polynominal_predict);

	gd.SetAlpha(0.4);   // value for momentum
	gd.SetEps(0.000001); // value to derivative 
	gd.SetEta_FirstJump(10);
	gd.SetEta_k_inc(1.1);
	gd.SetEta_k_dec(2.0);
	gd.SetMin_Eta(1e-10);      // min descent rate (will be multiplied by FirstJump before get started)
	gd.SetFinDifMethod(true);  // using a method of finite differences
	gd.SetMaxIters(5000);      // iteration limit
	gd.SetMaxTime(3);          // time limit (seconds)

	const int param_count = 4;
	vector<double> params(param_count, 0);             // all start parameters are 0
	vector<double> min_constrains(param_count, -10);   // all min_constrains are -10
	vector<double> max_constrains(param_count,  10);   // all max constrains are  10
	vector<double> rel_constrains(param_count, 0);     // relative constrains in %, aren't used
	vector<bool> type_constrains(param_count, false);  // 'false' - to use absolute constrains

	gd.SetParams(params);
	gd.SetMinConstrains(min_constrains);
	gd.SetMaxConstrains(max_constrains);
	gd.SetRelConstrains(rel_constrains);
	gd.SetTypeConstrains(type_constrains);

	GradErrorType res = gd.Go();

	cout << "gd_polynominal result: ";
	switch (res)
	{
	case GradErrorType::Success:
		cout << "Success" << endl;
		break;
	case GradErrorType::VectorSizesNotTheSame:
		cout << "VectorSizesNotTheSame" << endl;
		break;
	case GradErrorType::CanceledByUser:
		cout << "CanceledByUser" << endl;
		break;
	case GradErrorType::TimeOut:
		cout << "TimeOut" << endl;
		break;
	case GradErrorType::ItersOverflow:
		cout << "ItersOverflow" << endl;
		break;
	}

	cout << "gd.GetLastCost() = " << gd.GetLastCost() << endl;

	params = gd.GetParams();

	cout << "Found Params: " << endl;
	for (auto p : params)
		cout << p << endl;
	
	BOOST_CHECK(CmpFunc(params[0],  0.5,  0.001)); // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[1], -1.1,  0.001)); // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[2],  0.75, 0.001)); // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[3], -1.5,  0.05));  // close to parameters of SrcFunction
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


double damped_oscillations_experimental(double x)
{
	double noise = rand() / (double)RAND_MAX / 250.0; // add some noise 
	return 3.0 * sin(0.25 * x + 0.5) * exp(-0.02 * x) + 10.0 + noise;
}

double damped_oscillations_predict(double x, const std::vector<double>& p)
{
	return p[0] * sin(p[1] * x + p[2]) * exp(-p[3] * x) + p[4];
}

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_gd_damped_oscillations_test)
{
	GradDescent gd;

	TableFunction experimental; // Generate some damped_oscillations experimental data
	experimental.CreateDemoFunction(101, -20, 1.0, damped_oscillations_experimental, "damped_oscillations_experimental");

	gd.SetSrcFunction(experimental);
	gd.SetDstFunction(damped_oscillations_predict);

	gd.SetAlpha(0.45);    // value for momentum
	gd.SetEps(0.000001);  // value to derivative 
	gd.SetEta_FirstJump(10);
	gd.SetEta_k_inc(1.09);
	gd.SetEta_k_dec(2.0);
	gd.SetMin_Eta(1e-11);       // min descent rate (will be multiplied by FirstJump before get started)
	gd.SetFinDifMethod(false);  // using a plain central derivative
	gd.SetMaxIters(1000);       // iteration limit
	gd.SetMaxTime(3);           // time limit (seconds)

	const int param_count = 5;
	vector<double> params(param_count);             
	vector<double> min_constrains(param_count);   
	vector<double> max_constrains(param_count);  
	vector<double> rel_constrains(param_count, 0);     
	vector<bool> type_constrains(param_count, false);

	params[0] = 2.5; params[1] = 0.27; params[2] = 0; params[3] = 0.01; params[4] = 15;

	min_constrains[0] = 1;       max_constrains[0] = 3;
	min_constrains[1] = 0;       max_constrains[1] = 0;    // won't be use
	min_constrains[2] = -3.15;   max_constrains[2] = 3.15;
	min_constrains[3] = 0.001;   max_constrains[3] = 0.05;
	min_constrains[4] = 0;       max_constrains[4] = 30;

	rel_constrains[1] = 15;      // 15 % by params[1]
	type_constrains[1] = true;

	gd.SetParams(params);
	gd.SetMinConstrains(min_constrains);
	gd.SetMaxConstrains(max_constrains);
	gd.SetRelConstrains(rel_constrains);
	gd.SetTypeConstrains(type_constrains);

	gd.SetCallBackFreq(50);  // call callback function every 50 iterations

	gd.SetCallback([&gd]() {
		cout << "Last Cost = " << gd.GetLastCost() << endl;
		});

	GradErrorType res = gd.Go();

	cout << "gd_damped_oscillations result: ";
	switch (res)
	{
	case GradErrorType::Success:
		cout << "Success" << endl;
		break;
	case GradErrorType::VectorSizesNotTheSame:
		cout << "VectorSizesNotTheSame" << endl;
		break;
	case GradErrorType::CanceledByUser:
		cout << "CanceledByUser" << endl;
		break;
	case GradErrorType::TimeOut:
		cout << "TimeOut" << endl;
		break;
	case GradErrorType::ItersOverflow:
		cout << "ItersOverflow" << endl;
		break;
	}

	cout << "gd.GetLastCost() = " << gd.GetLastCost() << endl;

	params = gd.GetParams();

	cout << "Found Params: " << endl;
	for (auto p : params)
		cout << p << endl;

	BOOST_CHECK(CmpFunc(params[0], 3.0,  0.00001)); // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[1], 0.25, 0.0001));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[2], 0.5,  0.0001));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[3], 0.02, 0.0001));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[4], 10.0, 0.005));   // close to parameters of SrcFunction
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


double two_gaussian_distribution_experimental(double x)
{
	double noise = rand() / (double)RAND_MAX / 1.5; // add some noise 
	return 7500 * exp(-(x - 8600) * (x - 8600) / (2.0 * 80 * 80)) + // the first gaussian distribution
		   2250 * exp(-(x - 8900) * (x - 8900) / (2.0 * 85 * 85)) + // the second gaussian distribution
		   -0.1 * x + 1200 +                                        // linear background
		   noise;
}

double two_gaussian_distribution_predict(double x, const std::vector<double>& p)
{
	return p[0] * exp(-(x - p[1]) * (x - p[1]) / (2.0 * p[2] * p[2])) + // the first gaussian distribution
		   p[3] * exp(-(x - p[4]) * (x - p[4]) / (2.0 * p[5] * p[5])) + // the second gaussian distribution
		   p[6] * x + p[7];                                             // linear background
}

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_gd_two_gaussian_distribution_test)
{
	GradDescent gd;

	TableFunction experimental; // Generate some two_gaussian_distribution experimental data
	experimental.CreateDemoFunction(101, 8200, 12, two_gaussian_distribution_experimental, "two_gaussian_distribution_experimental");

	gd.SetSrcFunction(experimental);
	gd.SetDstFunction(two_gaussian_distribution_predict);

	gd.SetAlpha(0.55);    // value for momentum
	gd.SetEps(0.00001);   // value to derivative 
	gd.SetEta_FirstJump(10);
	gd.SetEta_k_inc(1.2);
	gd.SetEta_k_dec(2.0);
	gd.SetMin_Eta(1e-10);       // min descent rate (will be multiplied by FirstJump before get started)
	gd.SetFinDifMethod(false);  // using a plain central derivative
	gd.SetMaxIters(10000);      // iteration limit
	gd.SetMaxTime(10);          // time limit (seconds)

	const int param_count = 8;
	vector<double> params(param_count);
	vector<double> min_constrains(param_count);
	vector<double> max_constrains(param_count);
	vector<double> rel_constrains(param_count, 0);
	vector<bool> type_constrains(param_count, false);

	params[0] = 7000; params[1] = 8700; params[2] = 75; // the first gaussian distribution
	params[3] = 2400; params[4] = 8850; params[5] = 90; // the second gaussian distribution
	params[6] = 0; params[7] = 1000;                    // linear background

	min_constrains[0] = 5000;       max_constrains[0] = 10000;
	min_constrains[1] = 8500;       max_constrains[1] = 8800;   
	min_constrains[2] = 60;         max_constrains[2] = 100;
	min_constrains[3] = 1500;       max_constrains[3] = 4200;
	min_constrains[4] = 8800;       max_constrains[4] = 9000;
	min_constrains[5] = 30;         max_constrains[5] = 100;
	min_constrains[6] = -1;         max_constrains[6] = 1;
	min_constrains[7] = -3000;      max_constrains[7] = 3000;

	gd.SetParams(params);
	gd.SetMinConstrains(min_constrains);
	gd.SetMaxConstrains(max_constrains);
	gd.SetRelConstrains(rel_constrains);
	gd.SetTypeConstrains(type_constrains);

	gd.SetCallBackFreq(200);  // call callback function every 200 iterations

	gd.SetCallback([&gd]() {
		cout << "Last Cost = " << gd.GetLastCost() << endl;
		});

	GradErrorType res = gd.Go();

	cout << "gd_two_gaussian_distribution result: ";
	switch (res)
	{
	case GradErrorType::Success:
		cout << "Success" << endl;
		break;
	case GradErrorType::VectorSizesNotTheSame:
		cout << "VectorSizesNotTheSame" << endl;
		break;
	case GradErrorType::CanceledByUser:
		cout << "CanceledByUser" << endl;
		break;
	case GradErrorType::TimeOut:
		cout << "TimeOut" << endl;
		break;
	case GradErrorType::ItersOverflow:
		cout << "ItersOverflow" << endl;
		break;
	}

	cout << "gd.GetLastCost() = " << gd.GetLastCost() << endl;

	params = gd.GetParams();

	cout << "Found Params: " << endl;
	for (auto p : params)
		cout << p << endl;


	BOOST_CHECK(CmpFunc(params[0], 7500, 0.3));   // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[1], 8600, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[2], 80,   0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[3], 2250, 0.3));   // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[4], 8900, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[5], 85,   0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[6], -0.1, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[7], 1200, 1));     // close to parameters of SrcFunction
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double mix_experimental(double x)
{
	double noise = rand() / (double)RAND_MAX / 100.0; // add some noise 
	return -0.1 * x * x * x -0.6 * x * x + 4.0 * x -200.0 +            // polynominal
		   10.0 * sin(1.5 * x + 1.0) * exp(-0.1 * x) +                 // damped oscillations
		   40 * exp(-(x - -2.0) * (x - -2.0) / (2.0 * 2.0 * 2.0)) + // gaussian distribution
		   noise; 
}


double mix_predict(double x, const std::vector<double>& p)
{
	return p[0] * x * x * x + p[1] * x * x + p[2] * x + p[3] +           // polynominal
		   p[4] * sin(p[5] * x + p[6]) * exp(-p[7] * x) +                // damped oscillations
		   p[8] * exp(-(x - p[9]) * (x - p[9]) / (2.0 * p[10] * p[10])); // gaussian distribution
}

BOOST_AUTO_TEST_CASE(tf_gd_lib_test_gd_mix_test)
{
	GradDescent gd;

	TableFunction experimental; // Generate some mix_distribution experimental data
	experimental.CreateDemoFunction(101, -10, 0.18, mix_experimental, "mix_experimental");

	gd.SetSrcFunction(experimental);
	gd.SetDstFunction(mix_predict);

	gd.SetAlpha(0.55);    // value for momentum
	gd.SetEps(0.00001);   // value to derivative 
	gd.SetEta_FirstJump(10);
	gd.SetEta_k_inc(1.15);
	gd.SetEta_k_dec(2.0);
	gd.SetMin_Eta(1e-12);       // min descent rate (will be multiplied by FirstJump before get started)
	gd.SetFinDifMethod(false);  // using a plain central derivative
	gd.SetMaxIters(10000);      // iteration limit
	gd.SetMaxTime(10);          // time limit (seconds)

	const int param_count = 11;
	vector<double> params(param_count);
	vector<double> min_constrains(param_count);
	vector<double> max_constrains(param_count);
	vector<double> rel_constrains(param_count, 0);
	vector<bool> type_constrains(param_count, false);

	// -0.1 * x * x * x - 0.6 * x * x + 4.0 * x - 200.0 +        // polynominal
	//	10.0 * sin(1.5 * x + 1.0) * exp(-0.1 * x) +              // damped oscillations
	//	40 * exp(-(x - -2.0) * (x - -2.0) / (2.0 * 2.0 * 2.0))   // gaussian distribution

	params[0] = 0;   params[1] = 0;     params[2] = 0;    params[3] = -235;  // polynominal
	params[4] = 7;   params[5] = 1.4;   params[6] = 1.7;  params[7] = 0.3;   // damped oscillations
	params[8] = 60;  params[9] = -3.5;  params[10] = 5;                      // gaussian distribution

	min_constrains[0] = -10;       max_constrains[0] = 10;
	min_constrains[1] = -10;       max_constrains[1] = 10;
	min_constrains[2] = -10;       max_constrains[2] = 10;
	min_constrains[3] = -300;      max_constrains[3] = -100;
	min_constrains[4] = 5;         max_constrains[4] = 20;
	min_constrains[5] = 1.3;       max_constrains[5] = 1.7;
	min_constrains[6] = -6;        max_constrains[6] = 6;
	min_constrains[7] = 0.001;     max_constrains[7] = 0.5;
	min_constrains[8] = 20;        max_constrains[8] = 100;
	min_constrains[9] = -5;        max_constrains[9] = 2;
	min_constrains[10] = 1;        max_constrains[10] = 10;

	gd.SetParams(params);
	gd.SetMinConstrains(min_constrains);
	gd.SetMaxConstrains(max_constrains);
	gd.SetRelConstrains(rel_constrains);
	gd.SetTypeConstrains(type_constrains);

	gd.SetCallBackFreq(200);  // call callback function every 200 iterations

	gd.SetCallback([&gd]() {
		cout << "Last Cost = " << gd.GetLastCost() << endl;
		});

	GradErrorType res = gd.Go();

	cout << "mix result: ";
	switch (res)
	{
	case GradErrorType::Success:
		cout << "Success" << endl;
		break;
	case GradErrorType::VectorSizesNotTheSame:
		cout << "VectorSizesNotTheSame" << endl;
		break;
	case GradErrorType::CanceledByUser:
		cout << "CanceledByUser" << endl;
		break;
	case GradErrorType::TimeOut:
		cout << "TimeOut" << endl;
		break;
	case GradErrorType::ItersOverflow:
		cout << "ItersOverflow" << endl;
		break;
	}

	cout << "gd.GetLastCost() = " << gd.GetLastCost() << endl;

	params = gd.GetParams();

	cout << "Found Params: " << endl;
	for (auto p : params)
		cout << p << endl;

// -0.1 * x * x * x - 0.6 * x * x + 4.0 * x - 200.0 +          // polynominal
//	10.0 * sin(1.5 * x + 1.0) * exp(-0.1 * x) +                // damped oscillations
//	40 * exp(-(x - -2.0) * (x - -2.0) / (2.0 * 2.0 * 2.0))     // gaussian distribution

	BOOST_CHECK(CmpFunc(params[0], -0.1, 0.03));   // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[1], -0.6, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[2],  4.0, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[3], -200, 0.03));   // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[4], 10.0, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[5],  1.5, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[6],  1.0, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[7],  0.1, 0.05));     // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[8],  40,  0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[9], -2.0, 0.05));  // close to parameters of SrcFunction
	BOOST_CHECK(CmpFunc(params[10], 2.0, 0.05));     // close to parameters of SrcFunction
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
