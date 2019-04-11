/*!
*  \file    meteorologytimeseries.cpp
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
*  \todo
*  \warning
*/
#include "stdafx.h"
#include "meteorologybc.h"
#include "element.h"
#include "rhemodel.h"
#include "core/datacursor.h"
#include "temporal/timeseries.h"


MeteorologyBC::MeteorologyBC(Element *startElement, Element *endElement, int variableIndex, RHEModel *model)
  : QObject(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_variableIndex(variableIndex),
    m_model(model)
{
  m_dataCursor = new DataCursor();
}

MeteorologyBC::~MeteorologyBC()
{
  delete m_dataCursor;
}

void MeteorologyBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void MeteorologyBC::prepare()
{

}

void MeteorologyBC::applyBoundaryConditions(double dateTime)
{
  double value = 0;

  if(m_timeSeries->numColumns() == (int)m_profile.size())
  {
    switch (m_variableIndex)
    {
      case 1:
        {
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
            {
              m_profile[i]->relativeHumidity = value;
            }
          }
        }
        break;
      case 2:
        {
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
            {
              m_profile[i]->airTemperature = value;
            }
          }
        }
        break;
      case 3:
        {
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
            {
              m_profile[i]->windSpeed = value;
            }
          }
        }
        break;
      case 4:
        {
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
            {
              m_profile[i]->landCoverTemperature = value;
            }
          }
        }
        break;
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      switch (m_variableIndex)
      {
        case 1:
          {
            for(size_t i = 0; i < m_profile.size(); i++)
            {
              m_profile[i]->relativeHumidity = value;
            }
          }
          break;
        case 2:
          {
            for(size_t i = 0; i < m_profile.size(); i++)
            {
              m_profile[i]->airTemperature = value;
            }
          }
          break;
        case 3:
          {
            for(size_t i = 0; i < m_profile.size(); i++)
            {
              m_profile[i]->windSpeed = value;
            }
          }
          break;
        case 4:
          {
            for(size_t i = 0; i < m_profile.size(); i++)
            {
              m_profile[i]->landCoverTemperature = value;
            }
          }
          break;
      }
    }
  }
}

void MeteorologyBC::clear()
{

}

Element *MeteorologyBC::startElement() const
{
  return m_startElement;
}

void MeteorologyBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *MeteorologyBC::endElement() const
{
  return m_endElement;
}

void MeteorologyBC::setEndElement(Element *element)
{
  m_endElement = element;
}

QSharedPointer<TimeSeries> MeteorologyBC::timeSeries() const
{
  return m_timeSeries;
}

void MeteorologyBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}
