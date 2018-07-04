/*!
*  \file    stmproject.cpp
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
*  This file and its associated files and libraries are free software;
*  you can redistribute it and/or modify it under the terms of the
*  Lesser GNU General Public License as published by the Free Software Foundation;
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

#include "rhemodel.h"
#include "rhecomponent.h"
#include "spatial/point.h"
#include "spatial/network.h"
#include "element.cpp"
#include "elementjunction.h"
#include "spatial/edge.h"
#include "iboundarycondition.h"

using namespace std;

RHEModel::RHEModel(RHEComponent *component)
  : QObject(component),
    m_timeStep(0.0001), //seconds
    m_maxTimeStep(0.5), //seconds
    m_printFrequency(10),
    m_currentPrintCount(0),
    m_flushToDiskFrequency(10),
    m_currentflushToDiskCount(0),
    m_verbose(false),
    m_waterDensity(1000.0), //kg/m^3
    m_cp(4184.0), //4187.0 J/kg/C
    m_extinctionCoefficient(0.9675),
    m_albedo(0.0),
    m_atmLWReflection(0.03),
    m_emissWater(0.97),
    m_stefanBoltzmannConst(5.67e-8),
    m_atmEmissCoeff(0.5),
    #ifdef USE_NETCDF
    m_outputNetCDF(nullptr),
    #endif
    m_retrieveCouplingDataFunction(nullptr),
    m_component(component)
{
}

RHEModel::~RHEModel()
{

  for(Element *element : m_elements)
    delete element;

  m_elements.clear();
  m_elementsById.clear();


  for(ElementJunction *elementJunction : m_elementJunctions)
    delete elementJunction;

  m_elementJunctions.clear();
  m_elementJunctionsById.clear();


  closeOutputFiles();

  for(IBoundaryCondition *boundaryCondition : m_boundaryConditions)
    delete boundaryCondition;

  m_boundaryConditions.clear();
}

double RHEModel::maxTimeStep() const
{
  return m_maxTimeStep;
}

void RHEModel::setMaxTimeStep(double timeStep)
{
  m_maxTimeStep = timeStep;
}

double RHEModel::currentTimeStep() const
{
  return m_timeStep;
}

double RHEModel::startDateTime() const
{
  return m_startDateTime;
}

void RHEModel::setStartDateTime(double dateTime)
{
  m_startDateTime = dateTime;
}

double RHEModel::endDateTime() const
{
  return m_endDateTime;
}

void RHEModel::setEndDateTime(double dateTime)
{
  m_endDateTime = dateTime;
}

double RHEModel::outputInterval() const
{
  return m_outputInterval;
}

void RHEModel::setOutputInterval(double interval)
{
  m_outputInterval = interval;
}

double RHEModel::currentDateTime() const
{
  return m_currentDateTime;
}

double RHEModel::waterDensity() const
{
  return m_waterDensity;
}

void RHEModel::setWaterDensity(double value)
{
  m_waterDensity = value;
}

double RHEModel::specificHeatCapacityWater() const
{
  return m_cp;
}

void RHEModel::setSpecificHeatCapacityWater(double value)
{
  m_cp = value;
}

int RHEModel::numElementJunctions() const
{
  return m_elementJunctions.size();
}

ElementJunction *RHEModel::addElementJunction(const string &id, double x, double y, double z)
{
  if(m_elementJunctionsById.find(id) == m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = new ElementJunction(id, x, y, z, this);
    eJunction->index = m_elementJunctions.size();
    m_elementJunctions.push_back(eJunction);
    m_elementJunctionsById[id] = eJunction;
    return eJunction;
  }

  return nullptr;
}

void RHEModel::deleteElementJunction(const string &id)
{
  std::unordered_map<string,ElementJunction*>::iterator eJIter =  m_elementJunctionsById.find(id) ;

  if(eJIter != m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = eJIter->second;
    m_elementJunctionsById.erase(eJIter);

    std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
    if(it != m_elementJunctions.end())
    {
      m_elementJunctions.erase(it);
    }

    delete eJunction;
  }
}

void RHEModel::deleteElementJunction(int index)
{
  ElementJunction *eJunction = m_elementJunctions[index];

  m_elementJunctionsById.erase(eJunction->id);

  std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
  if(it != m_elementJunctions.end())
    m_elementJunctions.erase(it);

  delete eJunction;
}

ElementJunction *RHEModel::getElementJunction(const string &id)
{
  return m_elementJunctionsById[id];
}

ElementJunction *RHEModel::getElementJunction(int index)
{
  return m_elementJunctions[index];
}

int RHEModel::numElements() const
{
  return m_elements.size();
}

Element *RHEModel::addElement(const string &id, ElementJunction *upStream, ElementJunction *downStream)
{
  if(upStream && downStream)
  {
    Element *element = new Element(id, upStream, downStream, this);
    element->index = m_elements.size();
    m_elements.push_back(element);
    m_elementsById[id] = element;
    return element;
  }

  return nullptr;
}

void RHEModel::deleteElement(const string &id)
{
  unordered_map<string,Element*>::iterator eIter = m_elementsById.find(id);

  if(eIter != m_elementsById.end())
  {
    Element *element = eIter->second;
    m_elementsById.erase(eIter);

    vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);
    if(it != m_elements.end())
      m_elements.erase(it);

    delete element;
  }
}

void RHEModel::deleteElement(int index)
{
  Element *element = m_elements[index];
  m_elementJunctionsById.erase(element->id);

  vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);

  if(it != m_elements.end())
    m_elements.erase(it);

  delete element;
}

Element *RHEModel::getElement(const string &id)
{
  return m_elementsById[id];
}

Element *RHEModel::getElement(int index)
{
  return m_elements[index];
}

RetrieveCouplingData RHEModel::retrieveCouplingDataFunction() const
{
  return m_retrieveCouplingDataFunction;
}

void RHEModel::setRetrieveCouplingDataFunction(RetrieveCouplingData retrieveCouplingDataFunction)
{
  m_retrieveCouplingDataFunction = retrieveCouplingDataFunction;
}

bool RHEModel::initialize(list<string> &errors)
{
  bool initialized = initializeInputFiles(errors) &&
                     initializeTimeVariables(errors) &&
                     initializeElements(errors) &&
                     initializeSolver(errors) &&
                     initializeOutputFiles(errors) &&
                     initializeBoundaryConditions(errors);


  if(initialized)
  {
    applyInitialConditions();
  }

  return initialized;
}

bool RHEModel::finalize(std::list<string> &errors)
{
  closeOutputFiles();

  for(IBoundaryCondition *boundaryCondition : m_boundaryConditions)
    delete boundaryCondition;

  m_boundaryConditions.clear();

  return true;
}

bool RHEModel::initializeTimeVariables(std::list<string> &errors)
{
  if(m_startDateTime >= m_endDateTime)
  {
    errors.push_back("End datetime must be greater than startdatetime");
    return false;
  }

  if( (m_endDateTime - m_startDateTime) *  86400.0 < m_maxTimeStep )
  {
    errors.push_back("Make sure timestep is less than the simulation interval");
    return false;
  }

  if( m_maxTimeStep <= 0)
  {
    errors.push_back("Make sure time steps are greater 0");
    return false;
  }


  m_currentDateTime = m_startDateTime;
  m_nextOutputTime = m_currentDateTime;

  m_currentPrintCount = 0;
  m_currentflushToDiskCount = 0;

  return true;
}

bool RHEModel::initializeElements(std::list<string> &errors)
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
    size_t numElements = elementJunction->incomingElements.size() + elementJunction->outgoingElements.size();

    switch (numElements)
    {
      case 0:
        elementJunction->junctionType = ElementJunction::NoElement;
        break;
      case 1:
        elementJunction->junctionType = ElementJunction::SingleElement;
        break;
      case 2:
        elementJunction->junctionType = ElementJunction::DoubleElement;
        break;
      default:
        elementJunction->junctionType = ElementJunction::MultiElement;
        break;
    }
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->index = i;
    element->initialize();
  }


  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
  }

  return true;
}

bool RHEModel::initializeSolver(std::list<string> &errors)
{
//  m_heatSolver->setSize(m_elements.size());
//  m_heatSolver->initialize();

//  for(size_t i = 0; i < m_soluteSolvers.size(); i++)
//  {
//    ODESolver *solver =  m_soluteSolvers[i];
//    solver->setSize(m_elements.size());
//    solver->initialize();
//  }

  return true;
}

bool RHEModel::initializeBoundaryConditions(std::list<string> &errors)
{
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->clear();
    boundaryCondition->findAssociatedGeometries();
    boundaryCondition->prepare();
  }

  return true;
}

bool RHEModel::findProfile(Element *from, Element *to, std::list<Element *> &profile)
{
  for(Element *outgoing : from->downstreamJunction->outgoingElements)
  {
    if(outgoing == to)
    {
      profile.push_back(from);
      profile.push_back(outgoing);
      return true;
    }
    else if(findProfile(outgoing, to, profile))
    {
      profile.push_front(from);
      return true;
    }
  }

  return false;
}

