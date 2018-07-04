/*!
*  \file    element.cpp
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

#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "rhemodel.h"

#include <math.h>

using namespace std;

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  RHEModel *model)
  : id(id),
    channelTemperature(0),
    length(0),
    channelDepth(0),
    channelWidth(0),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    relativeHumidity(0.0),
    saturationVaporPressureAir(0.0),
    vaporPressureAir(0.0),
    windSpeed(0.0),
    airTemperature(0.0),
    incomingSWSolarRadiation(0.0),
    skyViewFactor(1.0),
    landCoverEmiss(0.96),
    netSWSolarRadiation(0.0),
    backLWRadiation(0.0),
    sedNetSWSolarRadiation(0.0),
    atmosphericLWRadiation(0.0),
    landCoverLWRadiation(0.0),
    totalIncomingSolarRadiation(0.0),
    totalNetSWSolarRadiation(0.0),
    totalBackLWRadiation(0.0),
    totalSedNetSWSolarRadiation(0.0),
    totalAtmosphericLWRadiation(0.0),
    totalLandCoverLWRadiation(0.0),
    upstreamElement(nullptr),
    downstreamElement(nullptr),
    model(model)
{

  upstream->outgoingElements.insert(this);
  downstream->incomingElements.insert(this);

  x = (upstream->x +  downstream->x) / 2.0;
  y = (upstream->y +  downstream->y) / 2.0;
  z = (upstream->z +  downstream->z) / 2.0;

}

Element::~Element()
{
  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);
}

void Element::initialize()
{
  //set upstream and downstream elements
  setUpstreamElement();
  setDownStreamElement();

  netSWSolarRadiation = incomingSWSolarRadiation = backLWRadiation =
      sedNetSWSolarRadiation = atmosphericLWRadiation =
      landCoverLWRadiation = totalIncomingSolarRadiation =
      totalNetSWSolarRadiation = totalBackLWRadiation =
      totalSedNetSWSolarRadiation = totalAtmosphericLWRadiation =
      totalLandCoverLWRadiation = netMCRadiation = 0.0;

}

double Element::computeRadiation()
{
  computeNetSWSolarRadiation();
  computeBackLWRadiation();
  computeAtmosphericLWRadiation();
  computeLandCoverLWRadiation();

  netMCRadiation = netSWSolarRadiation + backLWRadiation + atmosphericLWRadiation + landCoverLWRadiation;

  return netMCRadiation;
}

void Element::computeNetSWSolarRadiation()
{
  netSWSolarRadiation = (1.0  - model->m_albedo) * incomingSWSolarRadiation * (1.0 - shadeFactor);
  sedNetSWSolarRadiation = netSWSolarRadiation * exp(-model->m_extinctionCoefficient * channelDepth);
  netSWSolarRadiation -= sedNetSWSolarRadiation;
}

void Element::computeBackLWRadiation()
{
  double channelTempK = channelTemperature  + 273.15;
  backLWRadiation = -model->m_emissWater * model->m_stefanBoltzmannConst * (channelTempK  * channelTempK * channelTempK * channelTempK);
}

void Element::computeAtmosphericLWRadiation()
{
  saturationVaporPressureAir = 0.61275 * exp(17.27 * airTemperature / (237.3 + airTemperature));
  vaporPressureAir = relativeHumidity * saturationVaporPressureAir / 100.0;

  double airTempK = airTemperature + 273.15;
  double emiss = model->m_atmEmissCoeff + 0.031 * sqrt(vaporPressureAir/ /*0.134580245*/ 0.133322387415);

  atmosphericLWRadiation = model->m_stefanBoltzmannConst * (airTempK * airTempK * airTempK * airTempK) *
                           emiss * (1.0 - model->m_atmLWReflection);

  //  double esatair = 4.596 * exp(17.27 * airTemperature / (237.3 + airTemperature));
  //  double eair = relativeHumidity * esatair / 100.0;
  //  atmosphericLWRadiation = 0.000000117 * (airTempK * airTempK * airTempK * airTempK) * (0.5 + 0.031 *  sqrt(eair)) * (1. - model->m_atmLWReflection);
  //  atmosphericLWRadiation *= 4.184 * 10000.0  / 86400.0;
}

void Element::computeLandCoverLWRadiation()
{
  double airTempK = airTemperature + 273.15;
  landCoverLWRadiation = landCoverEmiss * model->m_stefanBoltzmannConst * (1.0 - skyViewFactor) * (airTempK * airTempK * airTempK * airTempK);
}

void Element::computeRadiationBalance(double timeStep)
{
  totalIncomingSolarRadiation += incomingSWSolarRadiation * timeStep / 1000.0;
  totalNetSWSolarRadiation += netSWSolarRadiation * timeStep / 1000.0;
  totalSedNetSWSolarRadiation += sedNetSWSolarRadiation * timeStep / 1000.0;
  totalBackLWRadiation += backLWRadiation * timeStep / 1000.0;
  totalAtmosphericLWRadiation += atmosphericLWRadiation * timeStep / 1000.0;
  totalLandCoverLWRadiation += landCoverLWRadiation * timeStep / 1000.0;
}

void Element::setUpstreamElement()
{
  upstreamElement = nullptr;
  upstreamElementDirection = 1.0;

  if(upstreamJunction->junctionType == ElementJunction::DoubleElement)
  {
    for(Element *element : upstreamJunction->incomingElements)
    {
      if(element != this)
      {
        upstreamElement = element;
        upstreamElementDirection = 1.0;
        return;
      }
    }

    for(Element *element : upstreamJunction->outgoingElements)
    {
      if(element != this)
      {
        upstreamElement = element;
        upstreamElementDirection = -1.0;
        return;
      }
    }
  }
}

void Element::setDownStreamElement()
{
  downstreamElement = nullptr;
  downstreamElementDirection = 1.0;

  if(downstreamJunction->junctionType == ElementJunction::DoubleElement)
  {
    for(Element *element : downstreamJunction->outgoingElements)
    {
      if(element != this)
      {
        downstreamElement = element;
        downstreamElementDirection = 1.0;
        return;
      }
    }

    for(Element *element : upstreamJunction->incomingElements)
    {
      if(element != this)
      {
        downstreamElement = element;
        downstreamElementDirection = -1.0;
        return;
      }
    }
  }
}

