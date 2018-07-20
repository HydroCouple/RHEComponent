/*!
*  \file    timeseriesbc.h
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

#ifndef ABSTRACTTIMESERIESBC_H
#define ABSTRACTTIMESERIESBC_H

#include "iboundarycondition.h"
#include "rhecomponent_global.h"

#include <vector>
#include <unordered_map>
#include <list>

class RHEModel;
struct Element;
struct ElementJunction;

class RHECOMPONENT_EXPORT AbstractTimeSeriesBC : public virtual IBoundaryCondition
{

  public:

    AbstractTimeSeriesBC(RHEModel *model);

    virtual ~AbstractTimeSeriesBC();

    void clear() override;

    bool addValue(double dateTime, double value);

    bool remove(int index);

    size_t size() const;

    double dateTime(size_t index) const;

    double value(size_t index) const;

  protected:

    double interpolate(double dateTime, bool &worked);

  private:

    int findDateTimeIndex(double dateTime);

  protected:

    RHEModel *m_model;
    static const std::unordered_map<std::string,int> m_variableIndexes;

  private:
    int m_cursor;
    std::vector<double> m_dateTimes;
    std::vector<double> m_values;
};

#endif // ABSTRACTTIMESERIESBC_H
