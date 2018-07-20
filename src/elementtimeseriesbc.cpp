/*!
*  \file    junctiontimeseriesbc.cpp
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

#include "elementtimeseriesbc.h"
#include "element.h"

ElementTimeSeriesBC::ElementTimeSeriesBC(Element *element, int variableIndex, RHEModel *model)
  :AbstractTimeSeriesBC(model),
    m_element(element),
    m_variableIndex(variableIndex)
{

}

ElementTimeSeriesBC::~ElementTimeSeriesBC()
{

}

void ElementTimeSeriesBC::findAssociatedGeometries()
{

}

void ElementTimeSeriesBC::prepare()
{
//  switch (m_variableIndex)
//  {
//    case -1:
//      {
//        m_element->channelTemperature = true;
//      }
//      break;
//    default:
//      {
//        m_element->soluteConcs[m_variableIndex].isBC = true;
//      }
//      break;
//  }
}

void ElementTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableIndex)
    {
      case -1:
        {
          m_element->channelTemperature = value;
        }
        break;
      default:
        {
//          m_element->soluteConcs[m_variableIndex].value = value;
        }
        break;
    }
  }
}

Element *ElementTimeSeriesBC::element() const
{
  return m_element;
}

void ElementTimeSeriesBC::setElement(Element *element)
{
  m_element = element;
}
