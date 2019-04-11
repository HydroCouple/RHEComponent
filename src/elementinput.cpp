/*!
   *  \file    elementinput.h
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
#include "rhecomponent.h"
#include "elementinput.h"
#include "spatial/point.h"
#include "spatial/linestring.h"
#include "rhemodel.h"
#include "element.h"
#include "temporal/timedata.h"
#include "hydrocoupletemporal.h"

#include <QDebug>
#include <omp.h>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace HydroCouple::SpatioTemporal;
using namespace HydroCouple::Temporal;

ElementInput::ElementInput(const QString &id,
                           Dimension *timeDimension,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           VariableType varType,
                           RHEComponent *modelComponent)
  : TimeGeometryMultiInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
                                 valueDefinition, modelComponent),
    m_component(modelComponent),
    m_varType(varType)
{

}

bool ElementInput::addProvider(HydroCouple::IOutput *provider)
{

  if(AbstractMultiInput::addProvider(provider))
  {
    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;
    ITimeIdBasedComponentDataItem *timeIdBaseDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)))
    {
      if(timeGeometryDataItem->geometryCount())
      {
        std::vector<bool> mapped(timeGeometryDataItem->geometryCount(), false);

        for(int i = 0; i < geometryCount() ; i++)
        {
          HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

          if(lineString->pointCount())
          {
            HCPoint *p1 = lineString->pointInternal(0);
            HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

            for(int j = 0; j < timeGeometryDataItem->geometryCount() ; j++)
            {
              if(!mapped[j])
              {
                ILineString *lineStringProvider = dynamic_cast<ILineString*>(timeGeometryDataItem->geometry(j));

                IPoint *pp1 = lineStringProvider->point(0);
                IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);

                double deltap1p1 = hypot(p1->x() - pp1->x() , p1->y() - pp1->y());
                double deltap2p2 = hypot(p2->x() - pp2->x() , p2->y() - pp2->y());

                double deltap1p2 = hypot(p1->x() - pp2->x() , p1->y() - pp2->y());
                double deltap2p1 = hypot(p2->x() - pp1->x() , p2->y() - pp1->y());

                if( deltap1p1 < 1e-3 && deltap2p2 < 1e-3)
                {
                  m_geometryMapping[provider][i] = j;
                  //                  m_geometryMappingOrientation[provider][i] = 1.0;
                  mapped[j] = true;
                  break;
                }
                else if(deltap1p2 < 1e-3 &&  deltap2p1 < 1e-3)
                {
                  m_geometryMapping[provider][i] = j;
                  //                  m_geometryMappingOrientation[provider][i] = -1.0;
                  mapped[j] = true;
                  break;
                }
              }
            }
          }
        }
      }

      return true;
    }
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)))
    {
      if(geometryDataItem->geometryCount())
      {
        std::vector<bool> mapped(geometryDataItem->geometryCount(), false);

        for(int i = 0; i < geometryCount() ; i++)
        {
          HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

          if(lineString->pointCount())
          {
            HCPoint *p1 = lineString->pointInternal(0);
            HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

            for(int j = 0; j < geometryDataItem->geometryCount() ; j++)
            {
              if(!mapped[j])
              {
                ILineString *lineStringProvider = dynamic_cast<ILineString*>(geometryDataItem->geometry(j));

                IPoint *pp1 = lineStringProvider->point(0);
                IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);

                double deltap1p1 = hypot(p1->x() - pp1->x() , p1->y() - pp1->y());
                double deltap2p2 = hypot(p2->x() - pp2->x() , p2->y() - pp2->y());

                double deltap1p2 = hypot(p1->x() - pp2->x() , p1->y() - pp2->y());
                double deltap2p1 = hypot(p2->x() - pp1->x() , p2->y() - pp1->y());

                if( deltap1p1 < 1e-3 && deltap2p2 < 1e-3)
                {
                  m_geometryMapping[provider][i] = j;
                  mapped[j] = true;
                  break;
                }
                else if(deltap1p2 < 1e-3 &&  deltap2p1 < 1e-3)
                {
                  m_geometryMapping[provider][i] = j;
                  mapped[j] = true;
                  break;
                }
              }
            }
          }
        }
      }

      return true;
    }
    else if((timeIdBaseDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)))
    {

      QStringList identifiers = timeIdBaseDataItem->identifiers();

      if(identifiers.size())
      {
        std::vector<bool> mapped(identifiers.size(), false);

        for(int i = 0; i < geometryCount() ; i++)
        {
          HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

          for(int j = 0; j < identifiers.size(); j++)
          {
            QString identifier = identifiers[j];

            if(!mapped[j] && !identifier.compare(lineString->id()))
            {
              m_geometryMapping[provider][i] = j;
              mapped[j] = true;
              break;
            }
          }
        }
      }

      return true;
    }
  }

  return false;
}

bool ElementInput::removeProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::removeProvider(provider))
  {
    m_geometryMapping.erase(provider);
    //    m_geometryMappingOrientation.erase(provider);
    return true;
  }

  return false;
}

bool ElementInput::canConsume(HydroCouple::IOutput *provider, QString &message) const
{
  ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
  IGeometryComponentDataItem *geometryDataItem = nullptr;
  ITimeIdBasedComponentDataItem *timeIdBaseDataItem = nullptr;

  if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
     (timeGeometryDataItem->geometryType() == IGeometry::LineString ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZ ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZM) &&
     (provider->valueDefinition()->type() == QVariant::Double ||
      provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)) &&
          (geometryDataItem->geometryType() == IGeometry::LineString ||
           geometryDataItem->geometryType() == IGeometry::LineStringZ ||
           geometryDataItem->geometryType() == IGeometry::LineStringZM) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((timeIdBaseDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)))
  {
    return true;
  }

  message = "Provider must be a LineString";

  return false;
}

void ElementInput::retrieveValuesFromProvider()
{
  moveDataToPrevTime();
  int currentTimeIndex = m_times.size() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_component->modelInstance()->currentDateTime());

  for(HydroCouple::IOutput *provider : providers())
  {
    provider->updateValues(this);
  }
}

void ElementInput::applyData()
{
  double currentTime = m_component->modelInstance()->currentDateTime();

  for(HydroCouple::IOutput *provider : providers())
  {

    std::unordered_map<int,int>& geomap = m_geometryMapping[provider];
    //    std::unordered_map<int,double>& geoorient = m_geometryMappingOrientation[provider];

    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;
    ITimeIdBasedComponentDataItem *timeIdBaseDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)))
    {
      int currentTimeIndex = timeGeometryDataItem->timeCount() - 1;
      int previousTimeIndex = std::max(0 , timeGeometryDataItem->timeCount() - 2);

      double providerCurrentTime = timeGeometryDataItem->time(currentTimeIndex)->julianDay();
      double providerPreviousTime = timeGeometryDataItem->time(previousTimeIndex)->julianDay();

      if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
      {
        double factor = 0.0;

        if(providerCurrentTime > providerPreviousTime)
        {
          double denom = providerCurrentTime - providerPreviousTime;
          double numer = currentTime - providerPreviousTime;
          factor = numer / denom;
        }

        switch (m_varType)
        {
          case Depth:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelDepth = value2 + factor *(value1 - value2);
              }
            }
            break;
          case TopWidth:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelWidth = value2 + factor *(value1 - value2);
              }
            }
            break;
          case Temperature:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelTemperature = value2 + factor *(value1 - value2);
              }

            }
            break;
          case LCTemperature:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->landCoverTemperature = value2 + factor *(value1 - value2);
              }

            }
            break;
          case SkyviewFactor:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->skyViewFactor = value2 + factor *(value1 - value2);
              }

            }
            break;
          case ShadeFactor:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactor = value2 + factor *(value1 - value2);
              }
            }
            break;
          case ShadeFactorMultiplier:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactorMultiplier = value2 + factor *(value1 - value2);
              }
            }
            break;
        }
      }
      else
      {
        switch (m_varType)
        {
          case Depth:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelDepth = value;
              }
            }
            break;
          case TopWidth:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelWidth = value;
              }
            }
            break;
          case Temperature:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelTemperature = value;
              }
            }
            break;
          case LCTemperature:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->landCoverTemperature = value;
              }
            }
            break;
          case SkyviewFactor:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->skyViewFactor = value;
              }
            }
            break;
          case ShadeFactor:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactor = value;
              }
            }
            break;
          case ShadeFactorMultiplier:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactorMultiplier = value;
              }
            }
            break;
        }
      }
    }
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)))
    {
      switch (m_varType)
      {
        case Depth:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->channelDepth = value;
            }
          }
          break;
        case TopWidth:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->channelWidth = value;
            }
          }
          break;
        case Temperature:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->channelTemperature = value;
            }
          }
          break;
        case SkyviewFactor:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->skyViewFactor = value;
            }
          }
          break;
        case ShadeFactor:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->shadeFactor = value;
            }
          }
          break;
        case ShadeFactorMultiplier:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->shadeFactorMultiplier = value;
            }
          }
          break;
      }
    }
    else if((timeIdBaseDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)))
    {
      int currentTimeIndex = timeIdBaseDataItem->timeCount() - 1;
      int previousTimeIndex = std::max(0 , timeIdBaseDataItem->timeCount() - 2);

      double providerCurrentTime = timeIdBaseDataItem->time(currentTimeIndex)->julianDay();
      double providerPreviousTime = timeIdBaseDataItem->time(previousTimeIndex)->julianDay();

      if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
      {
        double factor = 0.0;

        if(providerCurrentTime > providerPreviousTime)
        {
          double denom = providerCurrentTime - providerPreviousTime;
          double numer = currentTime - providerPreviousTime;
          factor = numer / denom;
        }

        switch (m_varType)
        {
          case Depth:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBaseDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelDepth = value2 + factor *(value1 - value2);
              }
            }
            break;
          case TopWidth:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBaseDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelWidth = value2 + factor *(value1 - value2);
              }
            }
            break;
          case Temperature:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBaseDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelTemperature = value2 + factor *(value1 - value2);
              }

            }
            break;
          case SkyviewFactor:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBaseDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->skyViewFactor = value2 + factor *(value1 - value2);
              }

            }
            break;
          case ShadeFactor:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBaseDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactor = value2 + factor *(value1 - value2);
              }
            }
            break;
          case ShadeFactorMultiplier:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBaseDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactorMultiplier = value2 + factor *(value1 - value2);
              }
            }
            break;
        }
      }
      else
      {
        switch (m_varType)
        {
          case Depth:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelDepth = value;
              }
            }
            break;
          case TopWidth:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelWidth = value;
              }
            }
            break;
          case Temperature:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->channelTemperature = value;
              }
            }
            break;
          case SkyviewFactor:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->skyViewFactor = value;
              }
            }
            break;
          case ShadeFactor:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactor = value;
              }
            }
            break;
          case ShadeFactorMultiplier:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBaseDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->shadeFactorMultiplier = value;
              }
            }
            break;
        }
      }
    }

  }
}

ElementInput::VariableType ElementInput::variableType() const
{
  return m_varType;
}

void ElementInput::setVariableType(VariableType variableType)
{
  m_varType = variableType;
}


//ElementHeatSourceInput::ElementHeatSourceInput(const QString &id,
//                                               Dimension *timeDimension,
//                                               Dimension *geometryDimension,
//                                               ValueDefinition *valueDefinition,
//                                               SourceType srcType,
//                                               RHEComponent *modelComponent)
//  : TimeGeometryMultiInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
//                                 valueDefinition, modelComponent),
//    m_component(modelComponent),
//    m_srcType(srcType)
//{

//}

//bool ElementHeatSourceInput::addProvider(HydroCouple::IOutput *provider)
//{
//  if(AbstractMultiInput::addProvider(provider))
//  {
//    ITimeGeometryComponentDataItem *timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider);

//    std::unordered_map<int, int> geometryMapping;

//    if(timeGeometryDataItem->geometryCount())
//    {
//      std::vector<bool> mapped(timeGeometryDataItem->geometryCount(), false);

//      for(int i = 0; i < geometryCount() ; i++)
//      {
//        HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

//        if(lineString->pointCount())
//        {
//          HCPoint *p1 = lineString->pointInternal(0);
//          HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

//          for(int j = 0; j < timeGeometryDataItem->geometryCount() ; j++)
//          {
//            if(!mapped[j])
//            {
//              ILineString *lineStringProvider = dynamic_cast<ILineString*>(timeGeometryDataItem->geometry(j));

//              IPoint *pp1 = lineStringProvider->point(0);
//              IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);


//              if(hypot(p1->x() - pp1->x() , p1->y() - pp1->y()) < 1e-4 && hypot(p2->x() - pp2->x() , p2->y() - pp2->y()) < 1e-4)
//              {
//                geometryMapping[i] = j;
//                mapped[j] = true;
//                break;
//              }
//              else if(hypot(p1->x() - pp2->x() , p1->y() - pp2->y()) < 1e-4 && hypot(p2->x() - pp1->x() , p2->y() - pp1->y()) < 1e-4)
//              {
//                geometryMapping[i] = j;
//                mapped[j] = true;
//                break;
//              }
//            }
//          }
//        }
//      }
//    }

//    m_geometryMapping[provider] = geometryMapping;

//    return true;
//  }

//  return false;
//}

//bool ElementHeatSourceInput::removeProvider(HydroCouple::IOutput *provider)
//{
//  if(AbstractMultiInput::removeProvider(provider))
//  {
//    m_geometryMapping.erase(provider);
//    return true;
//  }

//  return false;
//}

//void ElementHeatSourceInput::retrieveValuesFromProvider()
//{
//  moveDataToPrevTime();
//  int currentTimeIndex = m_times.size() - 1;
//  m_times[currentTimeIndex]->setJulianDay(m_component->modelInstance()->currentDateTime());

//  QList<IOutput*>::iterator it;

//#ifdef USE_OPENMP
////#pragma omp parallel for
//#endif
//  for(it = m_providers.begin(); it !=  m_providers.end() ; it++)
//  {
//    (*it)->updateValues(this);
//  }
//}

//void ElementHeatSourceInput::applyData()
//{
//  double currentTime = m_component->modelInstance()->currentDateTime();

//  for(IOutput *provider : m_providers)
//  {

//    std::unordered_map<int,int> &geometryMapping = m_geometryMapping[provider];

//    ITimeGeometryComponentDataItem *timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider);

//    int currentTimeIndex = timeGeometryDataItem->timeCount() - 1;
//    int previousTimeIndex = std::max(0 , timeGeometryDataItem->timeCount() - 2);

//    double providerCurrentTime = timeGeometryDataItem->time(currentTimeIndex)->julianDay();
//    double providerPreviousTime = timeGeometryDataItem->time(previousTimeIndex)->julianDay();

//    if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
//    {
//      double factor = 0.0;

//      if(providerCurrentTime > providerPreviousTime)
//      {
//        double denom = providerCurrentTime - providerPreviousTime;
//        double numer =currentTime - providerPreviousTime;
//        factor = numer / denom;
//      }

//      switch (m_srcType)
//      {
//        case RadiativeFlux:
//          {
//            for(auto it : geometryMapping)
//            {
//              double value1 = 0;
//              double value2 = 0;

//              timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
//              timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

//              Element *element =  m_component->modelInstance()->getElement(it.first);
//              element->radiationFluxes += value2 + factor *(value1 - value2);

//            }
//          }
//          break;
//        case HeatFlux:
//          {
//            for(auto it : geometryMapping)
//            {
//              double value1 = 0;
//              double value2 = 0;

//              timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
//              timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

//              Element *element =  m_component->modelInstance()->getElement(it.first);
//              element->externalHeatFluxes += value2 + factor *(value1 - value2);
//            }
//          }
//          break;
//      }
//    }
//    else
//    {
//      switch (m_srcType)
//      {
//        case RadiativeFlux:
//          {
//            for(auto it : geometryMapping)
//            {
//              double value = 0;
//              timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value);

//              Element *element =  m_component->modelInstance()->getElement(it.first);
//              element->radiationFluxes += value;

//            }
//          }
//          break;
//        case HeatFlux:
//          {
//            for(auto it : geometryMapping)
//            {
//              double value = 0;
//              timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value);

//              Element *element =  m_component->modelInstance()->getElement(it.first);
//              element->externalHeatFluxes += value;
//            }
//          }
//          break;
//      }
//    }
//  }
//}

//ElementHeatSourceInput::SourceType ElementHeatSourceInput::sourceType() const
//{
//  return m_srcType;
//}

//void ElementHeatSourceInput::setSourceType(SourceType srcType)
//{
//  m_srcType = srcType;
//}
