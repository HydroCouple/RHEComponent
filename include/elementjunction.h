/*!
*  \file    elementjunction.h
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

#ifndef ELEMENJUNCTION_H
#define ELEMENJUNCTION_H

#include "rhecomponent_global.h"

#include <string>
#include <set>

class RHEModel;
class HCVertex;
struct Element;

/*!
 * \brief The ElementJunction struct represents the
 */
struct  RHECOMPONENT_EXPORT ElementJunction
{

    enum JunctionType
    {
      NoElement = 0,
      SingleElement = 1,
      DoubleElement = 2,
      MultiElement = 3
    };

    /*!
     * \brief ElementJunction
     * \param numsolutes - Number of solutes
     * \param model -
     */
    ElementJunction(const std::string &id, double x, double y, double z, RHEModel *model);

    /*!
     * \brief ~ElementJunction - Deletes the ElementJunction and its associated data.
     */
    ~ElementJunction();

    /*!
     * \brief id unique identifier for this junction
     */
    std::string id;

    /*!
     * \brief index
     */
    int index;

    /*!
     * \brief x location
     */
    double x;

    /*!
     * \brief y location
     */
    double y;

    /*!
     * \brief z location
     */
    double z;

    /*!
     * \brief incomingElements
     */
    std::set<Element*> incomingElements;

    /*!
     * \brief outgoingElements
     */
    std::set<Element*> outgoingElements;


    /*!
     * \brief junctionType
     */
    JunctionType junctionType;

    /*!
     * \brief model
     */
    RHEModel *model;


};

#endif // ELEMENJUNCTION_H
