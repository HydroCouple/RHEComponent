/*!
 *  \file    odesolver.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 *  Runge Kutta solves from
 *  This file and its associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2018
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */


#ifndef ODESOLVER_H
#define ODESOLVER_H

#include "rhecomponent_global.h"

#ifdef USE_CVODE
#include <cvode/cvode.h>
#endif

#include <vector>

class ODESolver;

/*!
 *
 */
typedef void (*ComputeDerivatives)(double t, double y[], double dydt[], void* userData);

/*!
 *
 */
typedef int (ODESolver::*Solve)(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData);


#ifdef USE_CVODE

struct RHECOMPONENT_EXPORT RedirectionData
{
    ComputeDerivatives deriv;
    void *userData;
};

#endif

class ODESolver
{

  public:

    enum SolverType
    {
      RK4,
      RKQS,

#ifdef USE_CVODE
      CVODE_ADAMS,
      CVODE_BDF
#endif
    };

    /*!
     * \brief ODESolver
     * \param maxSize
     * \param solverType
     */
    ODESolver(int size, SolverType solverType);

    ~ODESolver();

    void initialize();

    int size() const;

    void setSize(int size);

    SolverType solverType() const;

    void setSolverType(SolverType solverType);

    int maxIterations() const;

    void setMaxIterations(int iterations);

    int getIterations() const;

    int order() const;

    void setOrder(int order);

    double relativeTolerance() const;

    void setRelativeTolerance(double tolerance);

    double absoluteTolerance() const;

    void setAbsoluteTolerance(double tolerance);

    /*!
     * \brief solve
     * \param y
     * \param n
     * \param t
     * \param dt
     * \param yout
     * \param derivs
     * \param userData
     * \return
     */
    int solve(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData);

  private:

    /*!
     * \brief rk4 Given values for the variables y[1..n] and their derivatives dydx[1..n] known at t, use the
     fourth-order Runge-Kutta method to advance the solution over an interval dt and return the
     incremented variables as yout[1..n], which need not be a distinct array from y. The user
     supplies the routine derivs(t,y,dydt), which returns derivatives dydt at t.
     * \param y
     * \param n
     * \param t
     * \param dt
     * \param yout
     * \param derivs
     * \param userData
     * \return
     */
    int rk4(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData);

    /*!
     * \brief rkqsDriver Driver function for Runge-Kutta integration with adaptive
      stepsize control. Integrates starting n values in y[]
      from t to t+dt with accuracy eps. h1 is the initial stepsize
      guess and derivs is a user-supplied function that computes
      derivatives dy/dx of y. On completion, yout contains the
      new values of y at the end of the integration interval.
     * \param y
     * \param n
     * \param t
     * \param dt
     * \param yout
     * \param derivs
     * \param userData
     * \return
     */
    int rkqsDriver(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData);

    /*!
     * \brief rkqs Fifth-order Runge-Kutta integration step with monitoring of
      local truncation error to assure accuracy and adjust stepsize.
      Inputs are current value of t, trial step size (dtTry), and
      accuracy (eps). Outputs are stepsize taken (dtDid) and estimated
      next stepsize (dtNext). Also updated are the values y[].
     * \param t
     * \param dydt
     * \param yout
     * \param n
     * \param dtTry
     * \param dtDid
     * \param dtNext
     * \param derivs
     * \param userData
     * \return
     */
    int rkqs(double *t, double dydt[], double yout[], int n, double dtTry, double *dtDid, double *dtNext, ComputeDerivatives derivs, void* userData);

    /*!
     * \brief rkck Uses the Runge-Kutta-Cash-Karp method to advance y[] at x
       over stepsize h.
     * \param t
     * \param dydt
     * \param yout
     * \param n
     * \param dt
     * \param deriv
     * \param userData
     */
    void rkck(double t, double dydt[],  double yout[], int n, double dt, ComputeDerivatives deriv, void* userData);


#ifdef USE_CVODE
    /*!
     * \brief solveCVODE
     * \param y
     * \param n
     * \param t
     * \param dt
     * \param yout
     * \param derivs
     * \param userData
     */
    int solveCVODE(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData);

    static int ComputeDerivatives_CVODE(realtype t, N_Vector y, N_Vector dydt, void *user_data);
#endif

    void clearMemory();

  public:

    bool m_print = false;

  private:


    int m_size,
    m_maxSteps,
    m_order,
    m_currentIterations;

    //RK4 Parameters
    double m_safety,
    m_pgrow,
    m_pshrnk,
    m_errcon,
    m_relTol,
    m_absTol,
    *m_yscal,
    *m_yerr,
    *m_ytemp,
    *m_ak;

    SolverType m_solverType;
    Solve m_solver;

#ifdef USE_CVODE
    void* m_cvodeSolver;
    N_Vector m_cvy;
#endif

};







#endif // ODESOLVER_H
