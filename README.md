# Bisection
Bisection is written by C++ with package gsl.
```CPP
vector<double> find_zeros(gsl_function F, double init_x[2], double step, int pole_num) {
    double x0 = init_x[0];
    double xf = init_x[1];
    vector<double> sol_vec;

    for (double x = x0 + step; x <= xf; x += step) {
        double f_x0 = GSL_FN_EVAL(&F, x0);
        double f_x = GSL_FN_EVAL(&F, x);

        if (f_x0 * f_x < 0) {
            // Root exists in [x0, x], so use bisection
            double range[2] = {x0, x};
            double sol = bisection(F, range);

            if ( abs(GSL_FN_EVAL(&F, sol)) < 1) {
                sol_vec.push_back(sol);

                if (sol_vec.size() == static_cast<size_t>(pole_num)) {
                    break;
                }

                x0 = x; // Move to the next range
            } else {
                x0 = x; // Move to the next range even if no root is found
            }
        } else {
            // If the function does not change sign, continue to the next step
            continue;
        }
    }
    return sol_vec;
}
```
# Example
```CPP
#include <iostream>
#include <fstream>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <bisection.cpp>
using namespace std;

double energy_level(double ka, double n)
{
    return ka * sin(ka * n / 2.) - 5. * cos(ka) * cos(ka * (n / 2. - 1.));
}

struct my_f_params
{
    double n;
};
double find_energy_levels(double ka, void *p)
{
    struct my_f_params *params = (struct my_f_params *)p;
    double n = (params->n);
    return energy_level(ka, n);
}

int main()
{
    // struct my_f_params params = {5e0};
    vector<double> rv;
    int pole_num = 10;
    double step = 1e-2;
    gsl_function F;
    F.function = &find_energy_levels;
    // F.params = &params;
    double init_x[2] = {0e0, 10e0};

    ofstream ff("./energy_level.csv");
    ff.precision(17);
    for (int i = 0; i < pole_num; i++)
    {
        ff << "x" << i+1;
        if (i != pole_num - 1)
        {
            ff << ",";
        }
        else
        {
            ff << "\n";
        }
    }

    double n0 = 2e0;
    while (n0 <= 10e0)
    {
        struct my_f_params params = {n0};
        F.params = &params;
        rv = find_zeros(F, init_x, step, pole_num);

        for (int i = 0; i < pole_num; ++i)
        {   
            if (i < rv.size())
            {
                ff << rv[i];
                if (i != pole_num - 1)
                {
                    ff << ",";
                }
                else
                {
                    ff << "\n";
                }
            }
            else
            {
                ff << "NaN";
                if (i != pole_num - 1)
                {
                    ff << ",";
                }
                else
                {
                    ff << "\n";
                }
            }
            
        }
        n0 += 0.01;
    }
    ff.close();
}
```
# Run your file
```sh
clang++ file_name.cpp -o a -lgsl -lgslcblas -std=c++20
```