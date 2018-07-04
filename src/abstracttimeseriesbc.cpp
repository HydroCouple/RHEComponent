/*!
*  \file    timeseriesbc.cpp
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
*  \todo Test transport on branching networks
*  \warning
*/


#include "abstracttimeseriesbc.h"
#include "element.h"
#include "elementjunction.h"
#include "rhemodel.h"

using namespace std;

AbstractTimeSeriesBC::AbstractTimeSeriesBC(RHEModel *model)
  : m_model(model),
    m_cursor(0)
{

}

AbstractTimeSeriesBC::~AbstractTimeSeriesBC()
{

}

void AbstractTimeSeriesBC::clear()
{
  m_cursor = 0;
}

bool AbstractTimeSeriesBC::addValue(double dateTime, double value)
{
  if(m_dateTimes.size() == 0 || dateTime > m_dateTimes[m_dateTimes.size() -1])
  {
    m_dateTimes.push_back(dateTime);
    m_values.push_back(value);
    return true;
  }

  return false;
}

bool AbstractTimeSeriesBC::remove(int index)
{
  if(index < (int)m_dateTimes.size())
  {
    m_dateTimes.erase(m_dateTimes.begin() + index);
    m_values.erase(m_values.begin() + index);
    return true;
  }

  return false;
}

double AbstractTimeSeriesBC::interpolate(double dateTime, bool &worked)
{
  int index = findDateTimeIndex(dateTime);

  if(index > -1)
  {
    double interpValue = 0.0;

    if(index == (int)m_dateTimes.size() - 1)
    {
      interpValue = m_values[index];
    }
    else
    {
      double downDate = m_dateTimes[index];
      double upDate = m_dateTimes[index + 1];
      double downValue = m_values[index];
      double upValue = m_values[index + 1];
      interpValue = downValue + ((upValue - downValue) * (dateTime -  downDate)) / (upDate - downDate);
    }

    worked = true;
    return interpValue;
  }

  worked = false;
  return -999999999999999;
}

int AbstractTimeSeriesBC::findDateTimeIndex(double dateTime)
{
  if(dateTime >= m_dateTimes[0] && dateTime <= m_dateTimes[m_dateTimes.size() - 1])
  {
    for(size_t i = m_cursor; i < m_dateTimes.size() - 1; i++)
    {
      if(dateTime >= m_dateTimes[i] && dateTime <= m_dateTimes[i+1])
      {
        m_cursor = i;
        return m_cursor;
      }
    }

    //resetart
    for(size_t i = 0; i < m_dateTimes.size() - 1; i++)
    {
      if(dateTime >= m_dateTimes[i] && dateTime <= m_dateTimes[i+1])
      {
        m_cursor = i;
        return m_cursor;
      }
    }
  }

  return -1;
}





