/*!
*  \file    rhemodelio.cpp
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
#include "rhemodel.h"
#include "element.h"
#include "elementjunction.h"
#include "radiativefluxbc.h"
#include "hydraulicsbc.h"
#include "meteorologybc.h"
#include "elementbc.h"
#include "temporal/timedata.h"
#include "threadsafenetcdf/threadsafencfile.h"
#include "threadsafenetcdf/threadsafencdim.h"
#include "threadsafenetcdf/threadsafencatt.h"
#include "temporal/timeseries.h"

#include <QDir>
#include <QDate>
#include <cstdlib>
#include <errno.h>

#ifdef USE_NETCDF

using namespace netCDF;
using namespace netCDF::exceptions;

#endif

using namespace std;

bool RHEModel::verbose() const
{
  return m_verbose;
}

void RHEModel::setVerbose(bool verbose)
{
  m_verbose = verbose;
}

int RHEModel::printFrequency() const
{
  return m_printFrequency;
}

void RHEModel::setPrintFrequency(int printFreq)
{
  m_printFrequency = printFreq;
}

int RHEModel::flushToDiskFrequency() const
{
  return m_flushToDiskFrequency;
}

void RHEModel::setFlushToDiskFrequency(int diskFlushFrequency)
{
  m_flushToDiskFrequency = diskFlushFrequency;
}

QFileInfo RHEModel::inputFile() const
{
  return m_inputFile;
}

void RHEModel::setInputFile(const QFileInfo &inputFile)
{
  m_inputFile = inputFile;
}

QFileInfo RHEModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void RHEModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

QFileInfo RHEModel::outputNetCDFFile() const
{
  return m_outputNetCDFFileInfo;
}

void RHEModel::setOutputNetCDFFile(const QFileInfo &outputNetCDFFile)
{
  m_outputNetCDFFileInfo = outputNetCDFFile;
}

void RHEModel::printStatus()
{
  m_currentPrintCount++;

  if (m_currentPrintCount >= m_printFrequency)
  {

    printf("RHE TimeStep (s): %f\tDateTime: %f\t {TotalIncomingSWSolarRadiation: %g (KJ/m^2)\tTotalNetSWSolarRadiation: %g (KJ/m^2)\t"
           "TotalBackLWRadiation: %g (KJ/m^2)\tTotalSedNetSWSolarRadiation: %g (KJ/m^2)\tTotalAtmosphericLWRadiation: %g (KJ/m^2)\t"
           "TotalLandCoverLWRadiation: %g (KJ/m^2)}\n", m_timeStep, m_currentDateTime,
           m_totalIncomingSolarRadiation, m_totalNetSWSolarRadiation, m_totalBackLWRadiation, m_totalSedNetSWSolarRadiation,
           m_totalAtmosphericLWRadiation, m_totalLandCoverLWRadiation);

    m_currentPrintCount = 0;
  }
}

void RHEModel::saveAs(const QFileInfo &filePath)
{
  QFileInfo fileInfo;

  if (filePath.isRelative())
  {
    fileInfo = relativePathToAbsolute(filePath);
  }
  else
  {
    fileInfo = filePath;
  }

  QString file = fileInfo.absoluteFilePath();

  if (!file.isEmpty() && !file.isNull() && filePath.absoluteDir().exists() && fileInfo.isFile())
  {


  }
}

bool RHEModel::initializeInputFiles(list<string> &errors)
{

  if (QFile::exists(m_inputFile.absoluteFilePath()))
  {
    QFile file(m_inputFile.absoluteFilePath());

    if (file.open(QIODevice::ReadOnly))
    {
      m_outNetCDFVariablesOnOff.clear();
      m_timeSeries.clear();

      m_delimiters = QRegExp("(\\,|\\t|\\;|\\s+)");
      int currentFlag = -1;

      QTextStream streamReader(&file);
      int lineCount = 0;
      while (!streamReader.atEnd())
      {
        QString line = streamReader.readLine().trimmed();
        lineCount++;

        if (!line.isEmpty() && !line.isNull())
        {
          bool readSuccess = true;
          QString error = "";

          auto it = m_inputFileFlags.find(line.toStdString());

          if (it != m_inputFileFlags.cend())
          {
            currentFlag = it->second;
          }
          else if (!QStringRef::compare(QStringRef(&line, 0, 2), ";;"))
          {
            //commment do nothing
          }
          else
          {
            switch (currentFlag)
            {
              case 1:
                readSuccess = readInputFileOptionTag(line, error);
                break;
              case 2:
                readSuccess = readInputFileOutputTag(line, error);
                break;
              case 3:
                readSuccess = readInputFileElementJunctionsTag(line, error);
                break;
              case 4:
                readSuccess = readInputFileElementsTag(line, error);
                break;
              case 5:
                readSuccess = readInputFileElementBCTag(line, error);
                break;
              case 6:
                readSuccess = readInputFileHydraulicsTag(line, error);
                break;
              case 7:
                readSuccess = readInputFileRadiativeFluxesTag(line, error);
                break;
              case 8:
                readSuccess = readInputFileMeteorologyTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileTimeSeriesTag(line, error);
                break;
            }
          }

          if (!readSuccess)
          {
            errors.push_back("Line " + std::to_string(lineCount) + " : " + error.toStdString());
            file.close();
            return false;
            break;
          }
        }
      }

      file.close();
    }
  }

  return true;
}

bool RHEModel::initializeOutputFiles(list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool RHEModel::initializeCSVOutputFile(list<string> &errors)
{
  if (m_outputCSVFileInfo.isRelative())
  {
    m_outputCSVFileInfo = relativePathToAbsolute(m_outputCSVFileInfo);
  }

  QString file = m_outputCSVFileInfo.absoluteFilePath();

  if (!file.isEmpty() && !file.isNull() && !m_outputCSVFileInfo.absoluteDir().exists())
  {
    errors.push_back("Output shapefile directory does not exist: " + file.toStdString());
    return false;
  }

  if (!file.isEmpty() && !file.isNull())
  {
    if(m_outputCSVFileInfo.isDir())
      return true;

    if (m_outputCSVStream.device() == nullptr)
    {
      QFile *device = new QFile(file, this);
      m_outputCSVStream.setDevice(device);
    }

    if (m_outputCSVStream.device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_outputCSVStream.setRealNumberPrecision(10);
      m_outputCSVStream.setRealNumberNotation(QTextStream::SmartNotation);
      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Depth, Width, ChannelTemperature, "
                           "IncomingSWSolarRadiation, NetSWSolarRadiation, BackLWRadiation, SedNetSWSolarRadiation, AtmosphericLWRadiation, LandCoverLWRadiation,"
                           "TotalIncomingSWSolarRadiation, TotalNetSWSolarRadiation, TotalBackLWRadiation, TotalSedNetSWSolarRadiation, "
                           "TotalAtmosphericLWRadiation, TotalLandCoverLWRadiation" << endl;

    }

    m_outputCSVStream.flush();

    return true;
  }

  return false;
}

bool RHEModel::initializeNetCDFOutputFile(list<string> &errors)
{

#ifdef USE_NETCDF

  if (m_outputNetCDFFileInfo.isRelative())
  {
    m_outputNetCDFFileInfo = relativePathToAbsolute(m_outputNetCDFFileInfo);
  }

  if (m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() || m_outputNetCDFFileInfo.isDir())
  {
    return true;
  }
  else if (!m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() &&
           !m_outputNetCDFFileInfo.absoluteFilePath().isNull() &&
           !m_outputNetCDFFileInfo.absoluteDir().exists())
  {
    std::string message = "NetCDF output file directory does not exist: " + m_outputNetCDFFileInfo.absoluteFilePath().toStdString();
    errors.push_back(message);
    return false;
  }

  bool returnValue = false;


  closeOutputNetCDFFile();

  try
  {

    m_outNetCDFVariables.clear();

    m_outputNetCDF = new ThreadSafeNcFile(m_outputNetCDFFileInfo.absoluteFilePath().toStdString(), NcFile::replace);


    //time variable
    ThreadSafeNcDim timeDim =  m_outputNetCDF->addDim("time");
    ThreadSafeNcVar timeVar =  m_outputNetCDF->addVar("time", NcType::nc_DOUBLE, timeDim);
    timeVar.putAtt("long_name", "Time");
    timeVar.putAtt("standard_name", "time");
    timeVar.putAtt("calendar", "julian");
    m_outNetCDFVariables["time"] = timeVar;


    //Add element junctions
    ThreadSafeNcDim junctionDim =  m_outputNetCDF->addDim("element_junctions", m_elementJunctions.size());

    ThreadSafeNcVar junctionIdentifiers =  m_outputNetCDF->addVar("element_junction_id", NcType::nc_STRING, junctionDim);
    junctionIdentifiers.putAtt("long_name", "Element Junction Identifier");
    m_outNetCDFVariables["element_junction_id"] = junctionIdentifiers;

    ThreadSafeNcVar junctionX =  m_outputNetCDF->addVar("x", NcType::nc_FLOAT, junctionDim);
    junctionX.putAtt("long_name", "Junction X-Coordinate");
    junctionX.putAtt("units", "m");
    m_outNetCDFVariables["x"] = junctionX;

    ThreadSafeNcVar junctionY =  m_outputNetCDF->addVar("y", NcType::nc_FLOAT, junctionDim);
    junctionY.putAtt("long_name", "Junction Y-Zoordinate");
    junctionY.putAtt("units", "m");
    m_outNetCDFVariables["y"] = junctionY;

    ThreadSafeNcVar junctionZ =  m_outputNetCDF->addVar("z", NcType::nc_FLOAT, junctionDim);
    junctionZ.putAtt("long_name", "junction z-coordinate");
    junctionZ.putAtt("units", "m");
    m_outNetCDFVariables["z"] = junctionZ;

    float *vertx = new float[m_elementJunctions.size()];
    float *verty = new float[m_elementJunctions.size()];
    float *vertz = new float[m_elementJunctions.size()];
    char **junctionIds = new char *[m_elementJunctions.size()];

    //write other relevant junction attributes here.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elementJunctions.size(); i++)
    {
      ElementJunction *junction = m_elementJunctions[i];

      junctionIds[i] = new char[junction->id.size() + 1];
      strcpy(junctionIds[i], junction->id.c_str());

      vertx[i] = junction->x;
      verty[i] = junction->y;
      vertz[i] = junction->z;
    }

    junctionX.putVar(vertx);
    junctionY.putVar(verty);
    junctionZ.putVar(vertz);
    junctionIdentifiers.putVar(junctionIds);

    delete[] vertx;
    delete[] verty;
    delete[] vertz;

    for (size_t i = 0; i < m_elementJunctions.size(); i++)
    {
      delete[] junctionIds[i];
    }

    delete[] junctionIds;

    //Add Elements
    ThreadSafeNcDim elementsDim =  m_outputNetCDF->addDim("elements", m_elements.size());

    ThreadSafeNcVar elementIdentifiers =  m_outputNetCDF->addVar("element_id", NcType::nc_STRING, elementsDim);
    elementIdentifiers.putAtt("long_name", "Element Identifier");
    m_outNetCDFVariables["element_id"] = elementIdentifiers;

    ThreadSafeNcVar elementFromJunction =  m_outputNetCDF->addVar("from_junction", NcType::nc_INT64, elementsDim);
    elementFromJunction.putAtt("long_name", "Upstream Junction");
    m_outNetCDFVariables["from_junction"] = elementFromJunction;

    ThreadSafeNcVar elementToJunction =  m_outputNetCDF->addVar("to_junction", NcType::nc_INT64, elementsDim);
    elementToJunction.putAtt("long_name", "Downstream Junction");
    m_outNetCDFVariables["to_junction"] = elementToJunction;

    //    ThreadSafeNcVar elementsVar =  m_outputNetCDF->addVar("elements", NcType::NC_FLOAT, elementsDim);
    //    elementsVar.putAtt("long_name", "Distance");
    //    elementsVar.putAtt("units", "m");
    //    m_outNetCDFVariables["elements"] = elementsVar;

    int *fromJunctions = new int[m_elements.size()];
    int *toJunctions = new int[m_elements.size()];
    char **elementIds = new char *[m_elements.size()];
    //    float *els = new float[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      elementIds[i] = new char[element->id.size() + 1];
      strcpy(elementIds[i], element->id.c_str());

      fromJunctions[i] = element->upstreamJunction->index;
      toJunctions[i] = element->downstreamJunction->index;
      //      els[i] = element->distanceFromUpStreamJunction;

    }

    elementIdentifiers.putVar(elementIds);
    elementFromJunction.putVar(fromJunctions);
    elementToJunction.putVar(toJunctions);
    //    elementsVar.putVar(els);

    delete[] fromJunctions;
    delete[] toJunctions;
    //    delete[] els;

    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      delete[] elementIds[i];
    }

    delete[] elementIds;

    auto varOnOff = [this](const std::string& name) -> bool
    {
      return m_outNetCDFVariablesOnOff.find(name) != m_outNetCDFVariablesOnOff.end() ? m_outNetCDFVariablesOnOff[name] : true;
    };

    //hydraulics variables
    if((m_outNetCDFVariablesOnOff["depth"] = varOnOff("depth")))
    {
      ThreadSafeNcVar depthVar =  m_outputNetCDF->addVar("depth", "float",
                                                         std::vector<std::string>({"time", "elements"}));
      depthVar.putAtt("long_name", "Flow Depth");
      depthVar.putAtt("units", "m");
      m_outNetCDFVariables["depth"] = depthVar;

      m_outNetCDFVariablesIOFunctions["depth"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *depth = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          depth[i] = static_cast<float>(elements[i]->channelDepth);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), depth);
        delete[] depth;
      };
    }

    if((m_outNetCDFVariablesOnOff["width"] = varOnOff("width")))
    {
      ThreadSafeNcVar widthVar =  m_outputNetCDF->addVar("width", "float",
                                                         std::vector<std::string>({"time", "elements"}));
      widthVar.putAtt("long_name", "Flow Top Width");
      widthVar.putAtt("units", "m");
      m_outNetCDFVariables["width"] = widthVar;
      m_outNetCDFVariablesIOFunctions["width"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *width = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          width[i] = static_cast<float>(elements[i]->channelWidth);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), width);
        delete[] width;
      };
    }

    if((m_outNetCDFVariablesOnOff["channel_temperature"] = varOnOff("channel_temperature")))
    {
      ThreadSafeNcVar temperatureVar =  m_outputNetCDF->addVar("channel_temperature", "float",
                                                               std::vector<std::string>({"time", "elements"}));
      temperatureVar.putAtt("long_name", "Channel Temperature");
      temperatureVar.putAtt("units", "°C");
      m_outNetCDFVariables["channel_temperature"] = temperatureVar;
      m_outNetCDFVariablesIOFunctions["channel_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *channel_temperature = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          channel_temperature[i] = static_cast<float>(elements[i]->channelTemperature);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), channel_temperature);
        delete[] channel_temperature;
      };
    }

    if((m_outNetCDFVariablesOnOff["air_temperature"] = varOnOff("air_temperature")))
    {
      ThreadSafeNcVar airTemperatureVar =  m_outputNetCDF->addVar("air_temperature", "float",
                                                                  std::vector<std::string>({"time", "elements"}));
      airTemperatureVar.putAtt("long_name", "Air Temperature");
      airTemperatureVar.putAtt("units", "°C");
      m_outNetCDFVariables["air_temperature"] = airTemperatureVar;
      m_outNetCDFVariablesIOFunctions["air_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *air_temperature = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          air_temperature[i] = static_cast<float>(elements[i]->airTemperature);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), air_temperature);
        delete[] air_temperature;
      };
    }


    if((m_outNetCDFVariablesOnOff["landcover_temperature"] = varOnOff("landcover_temperature")))
    {
      ThreadSafeNcVar landCoverTemperatureVar =  m_outputNetCDF->addVar("landcover_temperature", "float",
                                                                        std::vector<std::string>({"time", "elements"}));
      landCoverTemperatureVar.putAtt("long_name", "Landcover Temperature");
      landCoverTemperatureVar.putAtt("units", "°C");
      m_outNetCDFVariables["landcover_temperature"] = landCoverTemperatureVar;
      m_outNetCDFVariablesIOFunctions["landcover_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *landcover_temperature = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          landcover_temperature[i] = static_cast<float>(elements[i]->landCoverTemperature);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), landcover_temperature);
        delete[] landcover_temperature;
      };
    }

    if((m_outNetCDFVariablesOnOff["incoming_shortwave_solar_radiation"] = varOnOff("incoming_shortwave_solar_radiation")))
    {
      ThreadSafeNcVar incomingSWSolarRadiationVar =  m_outputNetCDF->addVar("incoming_shortwave_solar_radiation", "float",
                                                                            std::vector<std::string>({"time", "elements"}));
      incomingSWSolarRadiationVar.putAtt("long_name", "Incoming Shortwave Solar Radiation");
      incomingSWSolarRadiationVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["incoming_shortwave_solar_radiation"] = incomingSWSolarRadiationVar;
      m_outNetCDFVariablesIOFunctions["air_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *incoming_shortwave_solar_radiation = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          incoming_shortwave_solar_radiation[i] = static_cast<float>(elements[i]->incomingSWSolarRadiation);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), incoming_shortwave_solar_radiation);
        delete[] incoming_shortwave_solar_radiation;
      };
    }

    if((m_outNetCDFVariablesOnOff["net_shortwave_solar_radiation"] = varOnOff("net_shortwave_solar_radiation")))
    {
      ThreadSafeNcVar netSWSolarRadiationVar =  m_outputNetCDF->addVar("net_shortwave_solar_radiation", "float",
                                                                       std::vector<std::string>({"time", "elements"}));
      netSWSolarRadiationVar.putAtt("long_name", "Net Shortwave Solar Radiation");
      netSWSolarRadiationVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["net_shortwave_solar_radiation"] = netSWSolarRadiationVar;
      m_outNetCDFVariablesIOFunctions["net_shortwave_solar_radiation"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *net_shortwave_solar_radiation = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          net_shortwave_solar_radiation[i] = static_cast<float>(elements[i]->netSWSolarRadiation);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), net_shortwave_solar_radiation);
        delete[] net_shortwave_solar_radiation;
      };
    }

    if((m_outNetCDFVariablesOnOff["back_longwave_radiation"] = varOnOff("back_longwave_radiation")))
    {
      ThreadSafeNcVar backLWRadiationVar =  m_outputNetCDF->addVar("back_longwave_radiation", "float",
                                                                   std::vector<std::string>({"time", "elements"}));
      backLWRadiationVar.putAtt("long_name", "Back Longwave Radiation");
      backLWRadiationVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["back_longwave_radiation"] = backLWRadiationVar;
      m_outNetCDFVariablesIOFunctions["air_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *back_longwave_radiation = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          back_longwave_radiation[i] = static_cast<float>(elements[i]->backLWRadiation);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), back_longwave_radiation);
        delete[] back_longwave_radiation;
      };
    }

    if((m_outNetCDFVariablesOnOff["sediment_net_shortwave_solar_radiation"] = varOnOff("sediment_net_shortwave_solar_radiation")))
    {
      ThreadSafeNcVar sedNetSWSolarRadiationVar =  m_outputNetCDF->addVar("sediment_net_shortwave_solar_radiation", "float",
                                                                          std::vector<std::string>({"time", "elements"}));
      sedNetSWSolarRadiationVar.putAtt("long_name", "Net Shortwave Solar Radiation Reaching Sediment");
      sedNetSWSolarRadiationVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["sediment_net_shortwave_solar_radiation"] = sedNetSWSolarRadiationVar;
      m_outNetCDFVariablesIOFunctions["air_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *sediment_net_shortwave_solar_radiation = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          sediment_net_shortwave_solar_radiation[i] = static_cast<float>(elements[i]->sedNetSWSolarRadiation);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), sediment_net_shortwave_solar_radiation);
        delete[] sediment_net_shortwave_solar_radiation;
      };
    }

    if((m_outNetCDFVariablesOnOff["atmospheric_longwave_radiation"] = varOnOff("atmospheric_longwave_radiation")))
    {
      ThreadSafeNcVar atmosphericLWRadiationVar =  m_outputNetCDF->addVar("atmospheric_longwave_radiation", "float",
                                                                          std::vector<std::string>({"time", "elements"}));
      atmosphericLWRadiationVar.putAtt("long_name", "Atmospheric Longwave Radiation");
      atmosphericLWRadiationVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["atmospheric_longwave_radiation"] = atmosphericLWRadiationVar;
      m_outNetCDFVariablesIOFunctions["air_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *atmospheric_longwave_radiation = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          atmospheric_longwave_radiation[i] = static_cast<float>(elements[i]->atmosphericLWRadiation);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), atmospheric_longwave_radiation);
        delete[] atmospheric_longwave_radiation;
      };
    }

    if((m_outNetCDFVariablesOnOff["landcover_longwave_radiation"] = varOnOff("landcover_longwave_radiation")))
    {
      ThreadSafeNcVar landCoverLWRadiationVar =  m_outputNetCDF->addVar("landcover_longwave_radiation", "float",
                                                                        std::vector<std::string>({"time", "elements"}));
      landCoverLWRadiationVar.putAtt("long_name", "Landcover Longwave Radiation");
      landCoverLWRadiationVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["landcover_longwave_radiation"] = landCoverLWRadiationVar;
      m_outNetCDFVariablesIOFunctions["landcover_longwave_radiation"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *landcover_longwave_radiation = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          landcover_longwave_radiation[i] = static_cast<float>(elements[i]->landCoverLWRadiation);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), landcover_longwave_radiation);
        delete[] landcover_longwave_radiation;
      };
    }

    if((m_outNetCDFVariablesOnOff["net_main_channel_radiation"] = varOnOff("net_main_channel_radiation")))
    {
      ThreadSafeNcVar sumMCRadiationVar =  m_outputNetCDF->addVar("net_main_channel_radiation", "float",
                                                                  std::vector<std::string>({"time", "elements"}));
      sumMCRadiationVar.putAtt("long_name", "Net Radiation Reaching Main Channel");
      sumMCRadiationVar.putAtt("units", "W/m^2");
      m_outNetCDFVariables["net_main_channel_radiation"] = sumMCRadiationVar;
      m_outNetCDFVariablesIOFunctions["net_main_channel_radiation"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *net_main_channel_radiation = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          net_main_channel_radiation[i] = static_cast<float>(elements[i]->netMCRadiation);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), net_main_channel_radiation);
        delete[] net_main_channel_radiation;
      };
    }

    if((m_outNetCDFVariablesOnOff["shade_factor"] = varOnOff("shade_factor")))
    {
      ThreadSafeNcVar shadeVar =  m_outputNetCDF->addVar("shade_factor", "float",
                                                         std::vector<std::string>({"time", "elements"}));
      shadeVar.putAtt("long_name", "Shade Factor");
      shadeVar.putAtt("units", "");
      m_outNetCDFVariables["shade_factor"] = shadeVar;
      m_outNetCDFVariablesIOFunctions["shade_factor"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *shade_factor = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          shade_factor[i] = static_cast<float>(elements[i]->shadeFactor);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), shade_factor);
        delete[] shade_factor;
      };
    }

    if((m_outNetCDFVariablesOnOff["shade_factor_multiplier"] = varOnOff("shade_factor_multiplier")))
    {
      ThreadSafeNcVar shadeMultVar =  m_outputNetCDF->addVar("shade_factor_multiplier", "float",
                                                             std::vector<std::string>({"time", "elements"}));
      shadeMultVar.putAtt("long_name", "Shade Factor Multiplier");
      shadeMultVar.putAtt("units", "");
      m_outNetCDFVariables["shade_factor_multiplier"] = shadeMultVar;
      m_outNetCDFVariablesIOFunctions["shade_factor_multiplier"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *shade_factor_multiplier = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          shade_factor_multiplier[i] = static_cast<float>(elements[i]->shadeFactorMultiplier);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), shade_factor_multiplier);
        delete[] shade_factor_multiplier;
      };
    }

    m_optionalOutputVariables.clear();
    m_optionalOutputVariables.reserve(m_outNetCDFVariablesOnOff.size());

    for (const auto& pair : m_outNetCDFVariablesOnOff)
    {
      if(pair.second)
        m_optionalOutputVariables.push_back(pair.first);
    }

    m_outputNetCDF->sync();

    returnValue = true;

  }
  catch (NcException &e)
  {
    std::string message = std::string(e.what());
    printf("%s\n", e.what());
    errors.push_back(message);
    returnValue = false;
  }

#endif

  return returnValue;
}

bool RHEModel::readInputFileOptionTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  std::string optionsFlag = options[0].toStdString();
  auto it = m_optionsFlags.find(optionsFlag);

  if (it != m_optionsFlags.end())
  {
    int optionsIndex = it->second;

    switch (optionsIndex)
    {
      case 1:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (SDKTemporal::DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_startDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 2:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (SDKTemporal::DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_endDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 3:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_outputInterval = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Report interval error";
            return false;
          }
        }
        break;
      case 4:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_maxTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Max time step error";
            return false;
          }
        }
        break;
      case 5:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_albedo = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Albedo error";
            return false;
          }
        }
        break;
      case 6:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_emissWater = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Water emmisivity error";
            return false;
          }
        }
        break;
      case 7:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_atmEmissCoeff = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Water emmisivity error";
            return false;
          }
        }
        break;
      case 8:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_atmLWReflection = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Atmospheric longwater reflection error";
            return false;
          }
        }
        break;
      case 9:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_stefanBoltzmannConst = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Boltzmann constant error";
            return false;
          }
        }
        break;
      case 10:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_verbose = QString::compare(options[1], "No", Qt::CaseInsensitive) && QString::compare(options[1], "False", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Verbose error";
            return false;
          }
        }
        break;
      case 11:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_flushToDiskFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Flush to disk frequency error";
            return false;
          }
        }
        break;
      case 12:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_printFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Print frequency error";
            return false;
          }
        }
        break;
      case 13:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_extinctionCoefficient = options[1].toDouble(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Extinction coefficient error";
            return false;
          }
        }
        break;
    }
  }

  return true;
}

bool RHEModel::readInputFileOutputTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  QString optionsFlag = options[0];

  if (options.size() == 2)
  {
    if (!QString::compare(optionsFlag, "csv", Qt::CaseInsensitive))
    {
      m_outputCSVFileInfo = QFileInfo(options[1]);
    }
    else if (!QString::compare(optionsFlag, "netcdf", Qt::CaseInsensitive))
    {
      m_outputNetCDFFileInfo = QFileInfo(options[1]);
    }
  }
  else
  {
    errorMessage = "Output file error";
    return false;
  }

  return true;
}

bool RHEModel::readInputFileElementJunctionsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];

    bool workedX;
    bool workedY;
    bool workedZ;

    double x = columns[1].toDouble(&workedX);

    double y = columns[2].toDouble(&workedY);

    double z = columns[3].toDouble(&workedZ);

    if (workedX && workedY && workedZ)
    {
      addElementJunction(id.toStdString(), x, y, z);
    }
    else
    {
      errorMessage = "Junctions error";
      return false;
    }
  }
  else
  {
    errorMessage = "Junctions error";
    return false;
  }

  return true;
}

bool RHEModel::readInputFileElementsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() >= 10)
  {
    QString id = columns[0];
    QString fromId = columns[1];
    QString toId = columns[2];

    auto fromIt = m_elementJunctionsById.find(fromId.toStdString());
    auto toIt = m_elementJunctionsById.find(toId.toStdString());

    if (fromIt != m_elementJunctionsById.end() &&
        toIt != m_elementJunctionsById.end())
    {
      ElementJunction *ej1 = m_elementJunctionsById[fromId.toStdString()];
      ElementJunction *ej2 = m_elementJunctionsById[toId.toStdString()];

      bool lengthOk;
      double length = columns[3].toDouble(&lengthOk);

      bool depthOk ;
      double depth = columns[4].toDouble(&depthOk);

      bool widthOk ;
      double width = columns[5].toDouble(&widthOk);

      bool tempOk ;
      double temp = columns[6].toDouble(&tempOk);

      bool shadeOk;
      double shade = columns[7].toDouble(&shadeOk);

      bool skyViewOk;
      double skyView = columns[8].toDouble(&skyViewOk);

      bool lcEmmisOk;
      double lcEmmis = columns[9].toDouble(&lcEmmisOk);

      if (lengthOk && depthOk && skyViewOk &&
          widthOk && lcEmmisOk && tempOk && shadeOk)
      {
        Element *element = addElement(id.toStdString(), ej1, ej2);
        element->length = length;
        element->channelDepth = depth;
        element->channelWidth = width;
        element->channelTemperature = temp;
        element->skyViewFactor = skyView;
        element->landCoverEmiss = lcEmmis;
        element->shadeFactor = shade;

        if(columns.size() > 10)
        {
          shade = columns[10].toDouble(&shadeOk);
          element->shadeFactorMultiplier = shade;
        }
      }
      else
      {
        errorMessage = "";
        return false;
      }
    }
    else
    {
      errorMessage = "Wrong upstream or downstream junction";
      return false;
    }
  }

  return true;
}

bool RHEModel::readInputFileElementBCTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      bool found = false;
      QString variable = columns[2].trimmed();
      QString type = columns[3];

      if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
      {
        if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
        {
          int variableIndex = -2;

          if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
          {
            variableIndex = -1;
          }


          bool valueOk ;
          double value = columns[4].toDouble(&valueOk);

          if (valueOk)
          {

            ElementBC *elementBC = new ElementBC(fromElement, toElement,
                                                 variableIndex, this);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);
            m_timeSeries[ts->id().toStdString()] = ts;

            elementBC->setTimeSeries(ts);
            m_boundaryConditions.push_back(elementBC);

            found = true;
          }
          else
          {
            errorMessage = "Temperature BC value is invalid";
            return false;
          }
        }


        if (!found)
        {
          return false;
        }
      }
      else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
      {
        int variableIndex = -2;

        if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
        {
          variableIndex = -1;
        }

        if (variableIndex > -2)
        {
          std::string tsId = columns[4].toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            ElementBC *elementBC = new ElementBC(fromElement, toElement, variableIndex, this);
            elementBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(elementBC);
          }
          else
          {
            errorMessage = "Specified BC timeseries does not exist";
            return false;
          }
        }
        else
        {
          errorMessage = "Variable specified for BC is not valid";
          return false;
        }
      }
    }
    else
    {
      errorMessage = "Boundary condition junction not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool RHEModel::readInputFileHydraulicsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString variableType = columns[2];
      QString valueType = columns[3];
      QString varValue = columns[4].trimmed();

      auto it = m_hydraulicVariableFlags.find(variableType.toStdString());

      if (it != m_hydraulicVariableFlags.end())
      {
        int variableIndex = it->second;

        if (!QString::compare(valueType, "VALUE", Qt::CaseInsensitive))
        {

          bool valueOk;
          double value =  varValue.toDouble(&valueOk);

          if (valueOk)
          {
            HydraulicsBC *hydraulicsBC = new HydraulicsBC(fromElement, toElement,variableIndex, this);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);
            m_timeSeries[ts->id().toStdString()] = ts;

            hydraulicsBC->setTimeSeries(ts);
            m_boundaryConditions.push_back(hydraulicsBC);
          }
          else
          {
            errorMessage = "Uniform hydraulics value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = varValue.toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if(tsIt != m_timeSeries.end())
          {
            HydraulicsBC *hydraulicsBC = new HydraulicsBC(fromElement, toElement, variableIndex, this);
            hydraulicsBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(hydraulicsBC);
          }
          else
          {
            errorMessage = "Specified hydraulics timeseries does not exist";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Variable specified for uniform hydraulic is incorrect";
        return false;
      }
    }
    else
    {
      errorMessage = "Uniform hydraulic condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool RHEModel::readInputFileRadiativeFluxesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString type = columns[2];
      QString varValue = columns[3].trimmed();

      if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
      {

        bool valueOk;
        double value =  varValue.toDouble(&valueOk);

        if (valueOk)
        {
          RadiativeFluxBC *radiationFluxBC = new RadiativeFluxBC(fromElement, toElement, this);

          QUuid uid = QUuid::createUuid();
          QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
          ts->addRow(m_startDateTime, value);
          ts->addRow(m_endDateTime, value);
          m_timeSeries[ts->id().toStdString()] = ts;

          radiationFluxBC->setTimeSeries(ts);
          m_boundaryConditions.push_back(radiationFluxBC);
        }
        else
        {
          errorMessage = "Radiation BC value is invalid";
          return false;
        }
      }
      else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
      {
        std::string tsId = varValue.toStdString();
        auto tsIt = m_timeSeries.find(tsId);

        if(tsIt != m_timeSeries.end())
        {
          RadiativeFluxBC *radiationFluxBC = new RadiativeFluxBC(fromElement, toElement, this);
          radiationFluxBC->setTimeSeries(tsIt->second);
          m_boundaryConditions.push_back(radiationFluxBC);
        }
        else
        {
          errorMessage = "Specified radiation flux timeseries does not exist";
          return false;
        }
      }
    }
    else
    {
      errorMessage = "Boundary condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool RHEModel::readInputFileMeteorologyTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString variableType = columns[2];
      QString valueType = columns[3];
      QString varValue = columns[4].trimmed();

      auto it = m_meteorologicalVariableFlags.find(variableType.toStdString());

      if (it != m_meteorologicalVariableFlags.end())
      {
        int variableIndex = it->second;

        if (!QString::compare(valueType, "VALUE", Qt::CaseInsensitive))
        {

          bool valueOk;
          double value =  varValue.toDouble(&valueOk);

          if (valueOk)
          {
            MeteorologyBC *meteorologyBC = new MeteorologyBC(fromElement, toElement,
                                                             variableIndex, this);
            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);
            m_timeSeries[ts->id().toStdString()] = ts;

            meteorologyBC->setTimeSeries(ts);
            m_boundaryConditions.push_back(meteorologyBC);

          }
          else
          {
            errorMessage = "Uniform meteorology value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = varValue.toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if(tsIt != m_timeSeries.end())
          {
            MeteorologyBC *meteorologyBC = new MeteorologyBC(fromElement, toElement, variableIndex, this);
            meteorologyBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(meteorologyBC);
          }
          else
          {
            errorMessage = "Specified meteorology timeseries does not exist";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Variable specified for uniform meteorology is incorrect";
        return false;
      }
    }
    else
    {
      errorMessage = "Uniform meteorology boundary condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool RHEModel::readInputFileTimeSeriesTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2)
  {
    QFileInfo fileInfo(options[1].trimmed());

    if (fileInfo.isRelative())
      fileInfo = relativePathToAbsolute(fileInfo);

    if(QFile::exists(fileInfo.absoluteFilePath()))
    {
      QSharedPointer<TimeSeries> timeSeries(TimeSeries::createTimeSeries(options[0], fileInfo, this));

      if(!timeSeries.isNull())
      {
        m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
      }
      else
      {
        errorMessage = "Timeseries specified is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified filepath does not exist";
      return false;
    }
  }
  else
  {
    errorMessage = "TimeSeries must have two columns";
    return false;
  }

  return true;
}

void RHEModel::writeOutput()
{
  m_currentflushToDiskCount++;

  if (m_currentflushToDiskCount >= m_flushToDiskFrequency)
  {
    m_flushToDisk = true;
    m_currentflushToDiskCount = 0;
  }
  else
  {
    m_flushToDisk = false;
  }

  writeCSVOutput();
  writeNetCDFOutput();
}

void RHEModel::writeCSVOutput()
{
  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    for (size_t i = 0; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];


      //      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Depth, Width, ChannelTemperature, "
      //                           "IncomingSWSolarRadiation, NetSWSolarRadiation, BackLWRadiation, SedNetSWSolarRadiation, AtmosphericLWRadiation, LandCoverLWRadiation,";
      //      "TotalIncomingSWSolarRadiation, TotalNetSWSolarRadiation, TotalBackLWRadiation, TotalSedNetSWSolarRadiation, TotalAtmosphericLWRadiation, TotalLandCoverLWRadiation" << endl;

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->index
                        << ", " << element->x << ", " << element->y << ", " << element->z
                        << ", " << element->channelDepth
                        << ", " << element->channelWidth
                        << ", " << element->channelTemperature
                        << ", " << element->incomingSWSolarRadiation
                        << ", " << element->netSWSolarRadiation
                        << ", " << element->backLWRadiation
                        << ", " << element->sedNetSWSolarRadiation
                        << ", " << element->atmosphericLWRadiation
                        << ", " << element->landCoverLWRadiation
                        << ", " << element->totalIncomingSolarRadiation
                        << ", " << element->totalNetSWSolarRadiation
                        << ", " << element->totalBackLWRadiation
                        << ", " << element->totalSedNetSWSolarRadiation
                        << ", " << element->totalAtmosphericLWRadiation
                        << ", " << element->totalLandCoverLWRadiation
                        << endl;
    }

    if (m_flushToDisk)
    {
      m_outputCSVStream.flush();
    }
  }
}

void RHEModel::writeNetCDFOutput()
{
#ifdef USE_NETCDF

  if(m_outputNetCDF)
  {

    size_t currentTime = m_outNetCDFVariables["time"].getDim(0).getSize();
    m_outNetCDFVariables["time"].putVar(std::vector<size_t>({currentTime}), m_currentDateTime);

    int nVars = static_cast<int>(m_optionalOutputVariables.size());

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nVars; i++)
    {
      std::string varName = m_optionalOutputVariables[static_cast<size_t>(i)];
      (m_outNetCDFVariablesIOFunctions[varName])(currentTime, m_outNetCDFVariables[varName], m_elements);
    }


    //    float *depth = new float[m_elements.size()];
    //    float *width = new float[m_elements.size()];
    //    float *temperature = new float[m_elements.size()];
    //    float *airTemperature = new float[m_elements.size()];
    //    float *landCoverTemperature = new float[m_elements.size()];
    //    float *incomingSWSolarRadiation = new float[m_elements.size()];
    //    float *netSWSolarRadiation = new float[m_elements.size()];
    //    float *backwaterLWRadiation = new float[m_elements.size()];
    //    float *sedNetSWSolarRadiation = new float[m_elements.size()];
    //    float *atmosphericLWRadiation = new float[m_elements.size()];
    //    float *landCoverLWRadiation = new float[m_elements.size()];
    //    float *sumMCRadiation = new float[m_elements.size()];
    //    float *shadeFactor = new float[m_elements.size()];
    //    float *shadeFactorMult = new float[m_elements.size()];

    //#ifdef USE_OPENMP
    //#pragma omp parallel for
    //#endif
    //    for (int i = 0; i < (int)m_elements.size(); i++)
    //    {
    //      Element *element = m_elements[i];
    //      depth[i] = element->channelDepth;
    //      width[i] = element->channelWidth;
    //      temperature[i] = element->channelTemperature;
    //      airTemperature[i] = element->airTemperature;
    //      landCoverTemperature[i] = element->landCoverTemperature;
    //      incomingSWSolarRadiation[i] = element->incomingSWSolarRadiation;
    //      netSWSolarRadiation[i] = element->netSWSolarRadiation;
    //      backwaterLWRadiation[i] = element->backLWRadiation;
    //      sedNetSWSolarRadiation[i] = element->sedNetSWSolarRadiation;
    //      atmosphericLWRadiation[i] = element->atmosphericLWRadiation;
    //      landCoverLWRadiation[i] = element->landCoverLWRadiation;
    //      sumMCRadiation[i] = element->netMCRadiation;
    //      shadeFactor[i] = element->shadeFactor;
    //      shadeFactorMult[i] = element->shadeFactorMultiplier;
    //    }

    //    m_outNetCDFVariables["depth"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), depth);

    //    m_outNetCDFVariables["width"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), width);

    //    m_outNetCDFVariables["channel_temperature"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), temperature);

    //    m_outNetCDFVariables["air_temperature"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), airTemperature);

    //    m_outNetCDFVariables["landcover_temperature"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), landCoverTemperature);

    //    m_outNetCDFVariables["net_main_channel_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), sumMCRadiation);

    //    m_outNetCDFVariables["incoming_shortwave_solar_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), incomingSWSolarRadiation);

    //    m_outNetCDFVariables["net_shortwave_solar_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), netSWSolarRadiation);

    //    m_outNetCDFVariables["back_longwave_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), backwaterLWRadiation);

    //    m_outNetCDFVariables["sediment_net_shortwave_solar_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), sedNetSWSolarRadiation);

    //    m_outNetCDFVariables["atmospheric_longwave_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), atmosphericLWRadiation);

    //    m_outNetCDFVariables["landcover_longwave_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), landCoverLWRadiation);

    //    m_outNetCDFVariables["shade_factor"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), shadeFactor);

    //    m_outNetCDFVariables["shade_factor_multiplier"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), shadeFactorMult);

    //    delete[] depth;
    //    delete[] width;
    //    delete[] temperature;
    //    delete[] airTemperature;
    //    delete[] landCoverTemperature;
    //    delete[] incomingSWSolarRadiation;
    //    delete[] netSWSolarRadiation;
    //    delete[] backwaterLWRadiation;
    //    delete[] sedNetSWSolarRadiation;
    //    delete[] atmosphericLWRadiation;
    //    delete[] landCoverLWRadiation;
    //    delete[] sumMCRadiation;
    //    delete[] shadeFactor;
    //    delete[] shadeFactorMult;

    if(m_flushToDisk)
    {
      m_outputNetCDF->sync();
    }
  }

#endif
}

void RHEModel::closeOutputFiles()
{
  closeCSVOutputFile();
  closeOutputNetCDFFile();
}

void RHEModel::closeCSVOutputFile()
{

  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }

}

void RHEModel::closeOutputNetCDFFile()
{
#ifdef USE_NETCDF
  if (m_outputNetCDF)
  {
    m_outputNetCDF->sync();
    delete m_outputNetCDF;
    m_outputNetCDF = nullptr;
  }
#endif
}

QFileInfo RHEModel::relativePathToAbsolute(const QFileInfo &fileInfo)
{
  if (fileInfo.isRelative())
  {
    if (!m_inputFile.filePath().isEmpty() &&
        !m_inputFile.filePath().isNull() &&
        QFile::exists(m_inputFile.absoluteFilePath()))
    {
      QFileInfo absoluteFilePath = m_inputFile.absoluteDir().absoluteFilePath(fileInfo.filePath());

      if (absoluteFilePath.absoluteDir().exists())
      {
        return absoluteFilePath;
      }
    }
  }

  return fileInfo;
}

const unordered_map<string, int> RHEModel::m_inputFileFlags({
                                                              {"[OPTIONS]", 1},
                                                              {"[OUTPUTS]", 2},
                                                              {"[ELEMENTJUNCTIONS]", 3},
                                                              {"[ELEMENTS]", 4},
                                                              {"[BOUNDARY_CONDITIONS]", 5},
                                                              {"[HYDRAULICS]", 6},
                                                              {"[RADIATIVE_FLUXES]", 7},
                                                              {"[METEOROLOGY]", 8},
                                                              {"[TIMESERIES]", 9},
                                                              {"[OUTPUTVARIABLES]", 10}
                                                            });

const unordered_map<string, int> RHEModel::m_optionsFlags({
                                                            {"START_DATETIME", 1},
                                                            {"END_DATETIME", 2},
                                                            {"REPORT_INTERVAL", 3},
                                                            {"TIME_STEP", 4},
                                                            {"ALBEDO", 5},
                                                            {"WATER_EMISSIVITY", 6},
                                                            {"ATM_EMISSIVITY_COEFF", 7},
                                                            {"ATM_LW_REFLECTION", 8},
                                                            {"STEFANBOLTZMANN_CONSTANT", 9},
                                                            {"VERBOSE", 10},
                                                            {"FLUSH_TO_DISK_FREQ", 11},
                                                            {"PRINT_FREQ", 12},
                                                            {"EXTINCTION_COEFF", 13},
                                                          });


const unordered_map<string, int> RHEModel::m_hydraulicVariableFlags({{"DEPTH", 1},
                                                                     {"WIDTH", 2},
                                                                    });

const unordered_map<string, int> RHEModel::m_meteorologicalVariableFlags({
                                                                           {"RELATIVE_HUMIDITY", 1},
                                                                           {"AIR_TEMPERATURE", 2},
                                                                           {"LANDCOVER_TEMPERATURE", 4},
                                                                         });

const QRegExp RHEModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");
