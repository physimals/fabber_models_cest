#pragma once
/** 
 * spline_interpolator.h - Utility classes for 1D cubic spline interpolation
 *
 * Copyright (C) 2017 University of Oxford 
 */

/*  CCOPYRIGHT */

#include <vector>
#include <string>

/** 
 * Parameters for cubic spline starting at x=x0
 * 
 * Equation is d(x-x0)^3 + c(x-x0)^2 + b(x-x0) + 1
 */
struct Spline {
    Spline(double a=0, double b=0, double c=0, double d=0, double x0=0)
      : a(a), b(b), c(c), d(d), x0(x0) {}

    /** Evaluate the equation */
    double operator()(double x) const;

    /** Return human-readable equation to plug into FooPlot etc. */
    std::string eqn() const;

    double a, b, c, d, x0;
};

/**
 * Base class for anything which does spline interpolation
 */
class SplineInterpolator
{
public:
    /** Evaluate the interpolated function at an arbitrary x value */
    double operator()(double x) const;
protected:
    std::vector<Spline> m_splines;
};

/**
 * PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) spline interpolator, 
 * as used in Matlab version of the model
 *
 * The PCHIP method is designed to be shape-preserving at the expense of smoothness.
 * The second derivative is not necessarily continuous, however 'overshoots' are
 * avoided.
 */
class PchipInterpolator : public SplineInterpolator
{
public:
    /**
     * Generate a PCHIP interpolation with control points given by vectors of x 
     * and y points
     */
    PchipInterpolator(std::vector<double> &x, std::vector<double> &y);
private:
    Spline get_spline(double x1, double y1, double m1, double x2, double y2, double m2);
    double edge_case(double h0, double h1, double m0, double m1);
    std::vector<double> get_derivatives(std::vector<double> &x, std::vector<double> &y);
};

/**
 * 'Natural' spline interpolator for comparison
 *
 * Natural splines ensure continuous second derivateives, which are set to zero at the
 * endpoints of the interval of interpolation.
 */
class NaturalSplineInterpolator : public SplineInterpolator
{
public:
    /**
     * Generate a natural spline interpolation with control points given by vectors of x 
     * and y points
     */
    NaturalSplineInterpolator(std::vector<double> &x, std::vector<double> &y);
};
