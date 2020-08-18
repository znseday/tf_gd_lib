//          Copyright Sergey Tsynikin 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <array>

#include "UnitGradDescent.h"

using namespace std;
using namespace tf_gd_lib;

void GradDescent::CalcCost()
{
    LastCost = 0;

	if (IsUseUserTargetFunction)
	{
		LastCost = UserTargetFunction(Params);
	}
	else
	{
		double dfC;
		for (size_t i = 0; i < SrcFunction.Size(); ++i)
		{
			dfC = SrcFunction.GetY(i) - DstFunction(SrcFunction.GetX(i), Params);
			LastCost += dfC * dfC;
		}
	}
}
//---------------------------------------------------------------------------

GradErrorType GradDescent::Go()
{
	TimeStart = ClockType::now();

    size_t ParamsCount = Params.size();

    array<size_t,4> sizes = {
    	MinConstrains.size(),
    	MaxConstrains.size(),
    	RelConstrains.size(),
    	TypeConstrains.size() }; // is ok to initialize in rinetime?

    if ( !all_of(sizes.begin(), sizes.end(), [ParamsCount](size_t v){return v == ParamsCount;}) )
    {
        return GradErrorType::VectorSizesNotTheSame;
    }

    for (size_t i = 0; i < ParamsCount; ++i)
    {
        if (TypeConstrains[i])
        {
            //MinConstrains[i] = Params[i] - Params[i]*RelConstrains[i]/100.0;
            //MaxConstrains[i] = Params[i] + Params[i]*RelConstrains[i]/100.0;
			MinConstrains[i] = Params[i] * (100 - RelConstrains[i]) / 100.0;
			MaxConstrains[i] = Params[i] * (100 + RelConstrains[i]) / 100.0;
        }
    }

    Cur_Eta.resize(ParamsCount);
    fill(Cur_Eta.begin(), Cur_Eta.end(), Min_Eta * Eta_FirstJump);

    IsCalculating = true;

    vector<double> dCost_dp(ParamsCount);

    vector<double> dp(ParamsCount, 0); // fill(dp.begin(), dp.end(), 0.0);

    vector<double> old_p(ParamsCount);

    //size_t PointsCount = SrcFunction.Size();

    double OldCost;
    double dCost;

    CalcCost(); // calc Cost in the first time
    OldCost = LastCost;

    LastIters = 0;

    do
    {
        if (!IsCalculating)
        {
            return GradErrorType::CanceledByUser;
        }

        fill(dCost_dp.begin(), dCost_dp.end(), 0.0);

    	for (size_t j = 0; j < ParamsCount; ++j)
        {
            if (FinDifMethod)
            {
                Params[j] -= 2*Eps;

                CalcCost();
                double yL2 = LastCost;

                Params[j] += Eps;
                CalcCost();
                double yL = LastCost;

                Params[j] += 2*Eps;
                CalcCost();
                double yR = LastCost;

                Params[j] += Eps;
                CalcCost();
                double yR2 = LastCost;

                Params[j] -= 2*Eps;
                CalcCost();

                dCost_dp[j] = (yL2 - 8*yL + 8*yR - yR2)/(12.0*Eps);
            }
            else
            {
                Params[j] -= Eps;
                CalcCost();

                double yL = LastCost;
                Params[j] += 2*Eps;
                CalcCost();

                double yR = LastCost;
                Params[j] -= Eps; 
                CalcCost();

                dCost_dp[j] = (-yL+yR)/(2.0*Eps);
            }
        }

        old_p = Params;

        for (size_t j = 0; j < ParamsCount; ++j)
        {
            Params[j] -= Cur_Eta[j]*dCost_dp[j];
            Params[j] += Alpha*dp[j];

            if (Params[j] > MaxConstrains[j])
                Params[j] = MaxConstrains[j];
            if (Params[j] < MinConstrains[j])
                Params[j] = MinConstrains[j];

            dp[j] = Params[j] - old_p[j];

            OldCost = LastCost;
            CalcCost();
            dCost = OldCost - LastCost;

			if (dCost > 0)
			{
				Cur_Eta[j] *= Eta_k_inc;
			}
            else
            {
                if (Cur_Eta[j] > Min_Eta)
                {
                    Params[j] = old_p[j];
                    dp[j] = 0;

                    Cur_Eta[j] /= Eta_k_dec;          	
                }
            }

        }

        ++LastIters;

        if (LastIters > MaxIters)
        {
            IsCalculating = false;
            return GradErrorType::ItersOverflow;
        }

        if (LastIters % CallBackFreq == 0)
        {
            TimeEnd = ClockType::now();
            LastTime = (double)chrono::duration_cast<chrono::milliseconds>(TimeEnd - TimeStart).count();
            LastTime /= 1.0e3;

            if (LastTime > MaxTime)
            {
                IsCalculating = false;
            	return GradErrorType::TimeOut;
            }

            if (Callback)
            {
                Callback();
            }
        }

    }
    while ( any_of(Cur_Eta.begin(), Cur_Eta.end(), [this](double v){return v > Min_Eta;}) );

    TimeEnd = ClockType::now();
    LastTime = (double)chrono::duration_cast<chrono::milliseconds>(TimeEnd - TimeStart).count();
    LastTime /= 1.0e3;

    IsCalculating = false;
    return GradErrorType::Success;
}
//---------------------------------------------------------------------------

