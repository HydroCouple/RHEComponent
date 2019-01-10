/*!
*  \file    radiativefluxtimeseriesbc.cpp
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
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
*  \todo Test transport on branching networks
*  \warning
*/

#include "stdafx.h"
#include "radiativefluxbc.h"
#include "rhemodel.h"
#include "element.h"
#include "core/datacursor.h"
#include "temporal/timeseries.h"

RadiativeFluxBC::RadiativeFluxBC(Element *startElement, Element *endElement, RHEModel *model)
  : QObject(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_model(model)
{
  m_dataCursor = new DataCursor();
}

RadiativeFluxBC::~RadiativeFluxBC()
{
  delete m_dataCursor;
}

void RadiativeFluxBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void RadiativeFluxBC::prepare()
{

}

void RadiativeFluxBC::applyBoundaryConditions(double dateTime)
{
  double value = 0;

  if(m_timeSeries->numColumns() == (int)m_profile.size())
  {
    for(size_t i = 0 ; i < m_profile.size(); i++)
    {
      if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
      {
        m_profile[i]->incomingSWSolarRadiation = value;
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t i = 0 ; i < m_profile.size(); i++)
      {
        m_profile[i]->incomingSWSolarRadiation = value;
      }
    }
  }
}

void RadiativeFluxBC::clear()
{

}

Element *RadiativeFluxBC::startElement() const
{
  return m_startElement;
}

void RadiativeFluxBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *RadiativeFluxBC::endElement() const
{
  return m_endElement;
}

void RadiativeFluxBC::setEndElement(Element *element)
{
  m_endElement = element;
}

QSharedPointer<TimeSeries> RadiativeFluxBC::timeSeries() const
{
  return m_timeSeries;
}

void RadiativeFluxBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}
