/*!
*  \file    variable.h
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

#ifndef VARIABLE_H
#define VARIABLE_H

#include "rhecomponent_global.h"

/*!
 * \brief The Variable struct represents a variable value
 */
struct RHECOMPONENT_EXPORT Variable
{

    /*!
     * \brief value - Magnitude of the variable
     */
    double value;

    /*!
     * \brief associatedValue - A related value that may be stored with this variable.
     */
    double associatedValue;

    /*!
     * \brief isBC - Boolean indicating whether this variable is a boundary condition or not.
     */
    bool isBC;

    Variable():
       value(0.0),
       associatedValue(0.0),
       isBC(false)
    {

    }

    void copy(const Variable & variable)
    {
      value = variable.value;
      associatedValue = variable.associatedValue;
      isBC = variable.isBC;
    }

};



#endif // VARIABLE_H
