/**
 * spline_interpolator.cc - Utility classes for 1D cubic spline interpolation
 *
 * Copyright (C) 2017 University of Oxford
 */

/*  CCOPYRIGHT */

#include "spline_interpolator.h"

#include "armawrap/newmat.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using namespace NEWMAT;

double Spline::operator()(double x) const
{
    double dx = x - x0;
    return a + b * dx + c * dx * dx + d * dx * dx * dx;
}

std::string Spline::eqn() const
{
    std::stringstream s;
    s << d << "(x-" << x0 << ")^3 + " << c << "(x-" << x0 << ")^2 + " << b << "(x-" << x0 << ") + " << a;
    return s.str();
}

double SplineInterpolator::operator()(double x) const
{
    int j;
    for (j = 0; j < (int)m_splines.size(); j++)
    {
        if (m_splines[j].x0 > x)
        {
            if (j == 0)
                j++;
            break;
        }
    }
    j--;

    return m_splines[j](x);
}

/**
 * Based on code found on StackOverflow, which was in turn derived from an algorithm on Wikipedia
 * How trustworthy is that? Not very as it turned out, since the original code had a memory bug (fixed
 * in below)
 *
 * https://stackoverflow.com/questions/1204553/are-there-any-good-libraries-for-solving-cubic-splines-in-c
 */
NaturalSplineInterpolator::NaturalSplineInterpolator(vector<double> &x, vector<double> &y)
{
    int n = x.size() - 1;
    if (n < 1)
    {
        throw std::length_error("Need at least two control points");
    }
    if (x.size() != y.size())
    {
        throw std::length_error("Number of x values must match number of y values");
    }

    vector<double> a;
    a.insert(a.begin(), y.begin(), y.end());
    vector<double> b(n);
    vector<double> d(n);
    vector<double> h;

    for (int i = 0; i < n; ++i)
        h.push_back(x[i + 1] - x[i]);

    vector<double> alpha;
    alpha.push_back(0); // This element is never used but alpha must have size n
    for (int i = 1; i < n; ++i)
    {
        alpha.push_back(3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1]);
    }

    vector<double> c(n + 1);
    vector<double> l(n + 1);
    vector<double> mu(n + 1);
    vector<double> z(n + 1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n; ++i)
    {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for (int j = n - 1; j >= 0; --j)
    {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / 3 / h[j];
    }

    for (int i = 0; i < n; ++i)
    {
        m_splines.push_back(Spline(a[i], b[i], c[i], d[i], x[i]));
    }
}

/**
 * Based on code and documentation for the SCIPY Python implementation of
 * PCHIP interpolation
 */
PchipInterpolator::PchipInterpolator(vector<double> &x, vector<double> &y)
{
    if (x.size() < 2)
    {
        throw std::length_error("Need at least two control points");
    }
    if (x.size() != y.size())
    {
        throw std::length_error("Number of x values must match number of y values");
    }

    // Gradients at internal points
    vector<double> dk = get_derivatives(x, y);

    for (unsigned int n = 0; n < x.size() - 1; n++)
    {
        double x1 = x[n];
        double x2 = x[n + 1];
        double y1 = y[n];
        double y2 = y[n + 1];
        double m1 = dk[n];
        double m2 = dk[n + 1];
        m_splines.push_back(get_spline(x1, y1, m1, x2, y2, m2));
    }
}

/**
 * Get the equation of the cubic curve through (x1, y1) with gradient m1 and (x2, y2) with gradient m2
 */
Spline PchipInterpolator::get_spline(double x1, double y1, double m1, double x2, double y2, double m2)
{
    Matrix A(4, 4);

    A << 0 << 0 << 0 << 1 << (x2 - x1) * (x2 - x1) * (x2 - x1) << (x2 - x1) * (x2 - x1) << (x2 - x1) << 1 << 0 << 0 << 1
      << 0 << 3 * (x2 - x1) * (x2 - x1) << 2 * (x2 - x1) << 1 << 0;

    ColumnVector B(4);
    B << y1 << y2 << m1 << m2;
    ColumnVector C = A.i() * B;
    Spline spl(C(4), C(3), C(2), C(1), x1);
    return spl;
}

/**
 * Determine the gradient at the edge points of the interval
 *
 * h0, h1 are the x-sizes of the first two intervals next to the edge
 * m0, m1 are the linear interpolation gradients of these intervals
 */
double PchipInterpolator::edge_case(double h0, double h1, double m0, double m1)
{
    // one-sided three-point estimate for the derivative
    double d = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);

    // If gradient has different sign to gradient in edge interval
    // set to zero to prevent overshoot and preserve general shape
    if (d * m0 < 0)
    {
        d = 0;
    }
    // Don't let the gradient get too extreme compared to the
    // linear gradient of the edge interval?
    else if ((m0 * m1 < 0) && (abs(d) > 3.0 * abs(m0)))
    {
        d = 3.0 * m0;
    }

    return d;
}

/**
 * Get the value of the derivative to apply at each control point.
 *
 * This is a harmonic mean of the linear interpolation gradients either side
 * of the point, except at edge points where a special case is used.
 */
vector<double> PchipInterpolator::get_derivatives(vector<double> &x, vector<double> &y)
{
    vector<double> dk;

    // Gradients at internal points
    vector<double> mk, hk;
    for (unsigned int k = 0; k < x.size() - 1; k++)
    {
        hk.push_back((x[k + 1] - x[k]));
        mk.push_back((y[k + 1] - y[k]) / (x[k + 1] - x[k]));
    }

    if (y.size() == 2)
    {
        // only have two points, use linear interpolation
        dk.push_back(mk[0]);
        dk.push_back(mk[0]);
        return dk;
    }

    for (unsigned int k = 0; k < mk.size() - 1; k++)
    {
        if ((mk[k] * mk[k + 1] < 0) || (mk[k] == 0) || (mk[k + 1] == 0))
        {
            // If gradients either side have different signs, or either is zero,
            // set gradient at point to zero to prevent overshoot and preserve
            // shape of curve
            dk.push_back(0);
        }
        else
        {
            double w1 = 2 * hk[k + 1] + hk[k];
            double w2 = hk[k + 1] + 2 * hk[k];

            double whmean = (w1 / mk[k] + w2 / mk[k + 1]) / (w1 + w2);
            dk.push_back(1 / whmean);
        }
    }

    // special case endpoints, as suggested in
    // Cleve Moler, Numerical Computing with MATLAB, Chap 3.4
    dk.insert(dk.begin(), edge_case(hk[0], hk[1], mk[0], mk[1]));
    dk.push_back(edge_case(hk[hk.size() - 1], hk[hk.size() - 2], mk[mk.size() - 1], mk[mk.size() - 2]));

    return dk;
}
