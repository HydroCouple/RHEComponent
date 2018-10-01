/*!
 *  \file    rhecomponentinfo.cpp
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
#include "rhecomponentinfo.h"
#include "rhecomponent.h"
#include "spatial/geometryfactory.h"

using namespace HydroCouple;

RHEComponentInfo::RHEComponentInfo(QObject *parent)
  :AbstractModelComponentInfo(parent)
{
  GeometryFactory::registerGDAL();

  setId("Radiative Heat Exchange Model 1.0.0");
  setCaption("RHE Component");
  setIconFilePath(":/RHEComponent/rhecomponenticon");
  setDescription("A one-dimensional channel heat and solute transport model.");
  setCategory("Hydrodyanmics\\Heat Transport");
  setCopyright("");
  setVendor("");
  setUrl("www.hydrocouple.org");
  setEmail("caleb.buahin@gmail.com");
  setVersion("1.0.0");

  QStringList documentation;
  documentation << "Several sources";
  setDocumentation(documentation);

}

RHEComponentInfo::~RHEComponentInfo()
{
}

IModelComponent *RHEComponentInfo::createComponentInstance()
{
  QString id =  QUuid::createUuid().toString();
  RHEComponent *component = new RHEComponent(id, this);
  component->setDescription("RHE Model Instance");
  return component;
}
