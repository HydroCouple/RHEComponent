/*!
*  \file    element.h
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

#ifndef ELEMENT_H
#define ELEMENT_H

#include "rhecomponent_global.h"

#include <string>

struct Element;
struct ElementJunction;
class RHEModel;


/*!
 * \brief This struct represents the channel control volume
 */
struct RHECOMPONENT_EXPORT Element
{
    /*!
    * \brief Element - Creates an instance of the control volume element used to represent a computational
    * element in a reach.
    * \param numSolutes - Number of solutes that are going to be transported in addition to temperature.
    * \param from - The upstream junction of this element.
    * \param to - The downstream junction of this element.
    * \param project
    */
    Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  RHEModel *model);

    /*!
    * \brief ~Element - Destructor for this class.
    */
    ~Element();

    /*!
    * \brief index unique identifier for element
    */
    int index;

    /*!
    * \brief id
    */
    std::string id;

    /*!
    * \brief x
    */
    double x;

    /*!
    * \brief y
    */
    double y;

    /*!
    * \brief z
    */
    double z;

    /*!
    * \brief temperature (C)
    */
    double channelTemperature;

    /*!
    * \brief length (m) not needed
    */
    double length;

    /*!
    * \brief depth (m)
    */
    double channelDepth;

    /*!
    * \brief width (m) not needed
    */
    double channelWidth;

    /*!
    * \brief upstreamJunction
    */
    ElementJunction *upstreamJunction;

    /*!
    * \brief toJunction
    */
    ElementJunction *downstreamJunction;

    /*!
    * \brief relativeHumidity (%)
    */
    double relativeHumidity;

    /*!
    * \brief saturationVaporPressure (kPa)
    */
    double saturationVaporPressureAir;

    /*!
    * \brief vaporPressure  (kPa)
    */
    double vaporPressureAir;

    /*!
    * \brief windVelocity  (m/s)
    */
    double windSpeed;

    /*!
    * \brief airTemperature (C)
    */
    double airTemperature;


    /*!
     * \brief landCoverTemperature (C)
     */
    double landCoverTemperature;

    /*!
    * \brief inComingSolarRadiation (W/m^2)
    */
    double incomingSWSolarRadiation;

    /*!
    * \brief skyViewFactor
    */
    double skyViewFactor;

    /*!
     * \brief shadeFactor
     */
    double shadeFactor;

    /*!
     * \brief shadeFactorMultiplier
     */
    double shadeFactorMultiplier;

    /*!
     * \brief landCoverEmiss
     */
    double landCoverEmiss;

    /*!
    * \brief netSolarRadiation  (W/m^2)
    */
    double netSWSolarRadiation;

    /*!
    * \brief backRadiation  (W/m^2)
    */
    double backLWRadiation;

    /*!
    * \brief sedNetSolarRadiation  (W/m^2)
    */
    double sedNetSWSolarRadiation;

    /*!
    * \brief atmosphericLWRadiation  (W/m^2)
    */
    double atmosphericLWRadiation;

    /*!
    * \brief landCoverLWRadiation
    */
    double landCoverLWRadiation;

    /*!
     * \brief sumRadiation
     */
    double netMCRadiation;

    /*!
    * \brief totalIncomingSolarRadiation (kJ/m^2)
    */
    double totalIncomingSolarRadiation;

    /*!
    * \brief netSolarRadiation  (kJ/m^2)
    */
    double totalNetSWSolarRadiation;

    /*!
    * \brief backRadiation  (kJ/m^2)
    */
    double totalBackLWRadiation;

    /*!
    * \brief sedNetSolarRadiation  (kJ/m^2)
    */
    double totalSedNetSWSolarRadiation;

    /*!
    * \brief atmosphericLWRadiation  (kJ/m^2)
    */
    double totalAtmosphericLWRadiation;

    /*!
    * \brief landCoverLWRadiation (kJ/m^2)
    */
    double totalLandCoverLWRadiation;

    /*!
    * \brief upstreamElement
    */
    Element *upstreamElement;

    /*!
    * \brief upstreamElementDirection
    */
    double upstreamElementDirection;

    /*!
    * \brief downstreamElement
    */
    Element *downstreamElement;

    /*!
    * \brief downstreamElementDirection
    */
    double downstreamElementDirection;


    /*!
     * \brief distanceFromUpStreamJunction
     */
    double distanceFromUpStreamJunction;

    /*!
    * \brief model
    */
    RHEModel *model;

    /*!
    * \brief initializeSolutes
    * \param numSolutes
    */
    void initialize();

    /*!
    * \brief computeRadiation - Computes the time derivative of temperature based on data generated by the ODE solver.
    * \param dt - The timestep over which to compute the solute gradient.
    * \param T - The temperature array for all elements.
    * \return
    */
    double computeRadiation();

    void computeNetSWSolarRadiation();

    void computeBackLWRadiation();

    void computeAtmosphericLWRadiation();

    void computeLandCoverLWRadiation();

    /*!
    * \brief computeHeatBalance
    */
    void computeRadiationBalance(double timeStep);

  private:

    /*!
    * \brief setUpstreamElement
    */
    void setUpstreamElement();

    /*!
    * \brief setDownStreamElement
    */
    void setDownStreamElement();

};

#endif // ELEMENT_H
