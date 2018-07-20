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

#include "meteorologytimeseriesbc.h"
#include "element.h"
#include "rhemodel.h"

MeteorologyTimeSeriesBC::MeteorologyTimeSeriesBC(Element *element, int variableIndex, RHEModel *model)
  : AbstractTimeSeriesBC(model),
    m_element(element),
    m_variableIndex(variableIndex)
{

}

MeteorologyTimeSeriesBC::~MeteorologyTimeSeriesBC()
{

}

void MeteorologyTimeSeriesBC::findAssociatedGeometries()
{

}

void MeteorologyTimeSeriesBC::prepare()
{

}

void MeteorologyTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableIndex)
    {
      case 1:
        {
          m_element->relativeHumidity = value;
        }
        break;
      case 2:
        {
          m_element->airTemperature = value;
        }
        break;
      case 3:
        {
          m_element->windSpeed = value;
        }
        break;
    }
  }
}

Element *MeteorologyTimeSeriesBC::element() const
{
  return m_element;
}

void MeteorologyTimeSeriesBC::setElement(Element *element)
{
  m_element = element;
}



UniformMeteorologyTimeSeriesBC::UniformMeteorologyTimeSeriesBC(Element *startElement, Element *endElement, int variableIndex, RHEModel *model)
  :AbstractTimeSeriesBC(model),
   m_startElement(startElement),
   m_endElement(endElement),
   m_variableIndex(variableIndex)
{

}

UniformMeteorologyTimeSeriesBC::~UniformMeteorologyTimeSeriesBC()
{

}

void UniformMeteorologyTimeSeriesBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void UniformMeteorologyTimeSeriesBC::prepare()
{

}

void UniformMeteorologyTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableIndex)
    {
      case 1:
        {
          for(Element *element : m_profile)
          {
            element->relativeHumidity = value;
          }
        }
        break;
      case 2:
        {
          for(Element *element : m_profile)
          {
            element->airTemperature = value;
          }
        }
        break;
      case 3:
        {
          for(Element *element : m_profile)
          {
            element->windSpeed = value;
          }
        }
        break;
    }
  }
}

Element *UniformMeteorologyTimeSeriesBC::startElement() const
{
  return m_startElement;
}

void UniformMeteorologyTimeSeriesBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *UniformMeteorologyTimeSeriesBC::endElement() const
{
  return m_endElement;
}

void UniformMeteorologyTimeSeriesBC::setEndElement(Element *element)
{
  m_endElement = element;
}

