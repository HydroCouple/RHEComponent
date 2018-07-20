#include "rhemodel.h"
#include "element.h"
#include "elementjunction.h"
#include "iboundarycondition.h"
#include "rhecomponent.h"

using namespace std;

void RHEModel::update()
{
  if(m_currentDateTime <= m_endDateTime)
  {

    applyBoundaryConditions(m_currentDateTime);

    //Retrieve external data from other coupled models
    if(m_retrieveCouplingDataFunction)
    {
      (*m_retrieveCouplingDataFunction)(this, m_currentDateTime);
    }

    if(m_component)
      m_component->applyInputValues();

    m_timeStep = computeTimeStep();

    //Solve heat transport first
    solveHeatTransport(m_timeStep);

    m_prevDateTime = m_currentDateTime;

    prepareForNextTimeStep();

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime = std::min(m_nextOutputTime + m_outputInterval / 86400.0 , m_endDateTime);
    }

    if(m_verbose)
    {
      printStatus();
    }

    m_currentDateTime = m_currentDateTime + m_timeStep / 86400.0;

  }
}

void RHEModel::prepareForNextTimeStep()
{

  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computeRadiationBalance(m_timeStep);
    m_totalIncomingSolarRadiation += element->totalIncomingSolarRadiation;
    m_totalNetSWSolarRadiation += element->totalNetSWSolarRadiation;
    m_totalBackLWRadiation += element->totalBackLWRadiation;
    m_totalSedNetSWSolarRadiation += element->totalSedNetSWSolarRadiation;
    m_totalAtmosphericLWRadiation +=  element->totalAtmosphericLWRadiation;
    m_totalLandCoverLWRadiation += element->totalLandCoverLWRadiation;
  }

}

void RHEModel::applyInitialConditions()
{

  //Initialize heat and solute balance trackers
  m_totalIncomingSolarRadiation = 0;
  m_totalNetSWSolarRadiation = 0;
  m_totalBackLWRadiation = 0;
  m_totalSedNetSWSolarRadiation = 0;
  m_totalAtmosphericLWRadiation =  0;
  m_totalLandCoverLWRadiation = 0;

  applyBoundaryConditions(m_currentDateTime);

  //Write initial output
  //writeOutput();

  //Set next output time
  m_nextOutputTime = m_currentDateTime; //m_outputInterval / 86400.0;
}

void RHEModel::applyBoundaryConditions(double dateTime)
{

//#ifdef USE_OPENMMP
//#pragma omp parallel for
//#endif
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->applyBoundaryConditions(dateTime);
  }
}

double RHEModel::computeTimeStep()
{
  double timeStep = m_maxTimeStep;
  double nextTime = m_currentDateTime + timeStep / 86400.0;

  if(nextTime > m_nextOutputTime)
  {
    timeStep = std::max(1.0,  (m_nextOutputTime - m_currentDateTime) *  86400.0);
  }

  timeStep = std::min(std::max(timeStep, 1.0), m_maxTimeStep);

  return timeStep;
}

void RHEModel::solveHeatTransport(double timeStep)
{
  //Set initial input and output values to current values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computeRadiation();
  }
}
