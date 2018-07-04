/*!
*  \file    radiativefluxtimeseriesbc.h
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

#ifndef RADIATIVEFLUXTIMESERIESBC_H
#define RADIATIVEFLUXTIMESERIESBC_H

#include "abstracttimeseriesbc.h"

class RHECOMPONENT_EXPORT RadiativeFluxTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    RadiativeFluxTimeSeriesBC(Element *element, RHEModel *model);

    virtual ~RadiativeFluxTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *element() const;

    void setElement(Element *element);

  private:

    Element *m_element;

};


class RHECOMPONENT_EXPORT UniformRadiativeFluxTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    UniformRadiativeFluxTimeSeriesBC(Element *startElement, Element *endElement, RHEModel *model);

    virtual ~UniformRadiativeFluxTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

  private:

    std::list<Element*> m_profile;
    Element *m_startElement, *m_endElement;

};



#endif // RADIATIVEFLUXTIMESERIESBC_H
