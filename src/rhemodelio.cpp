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

#include "rhemodel.h"
#include "element.h"
#include "elementjunction.h"
#include "radiativefluxtimeseriesbc.h"
#include "hydraulicstimeseriesbc.h"
#include "meteorologytimeseriesbc.h"
#include "elementtimeseriesbc.h"
#include "temporal/timedata.h"
#include "threadsafenetcdf/threadsafencfile.h"
#include "threadsafenetcdf/threadsafencdim.h"
#include "threadsafenetcdf/threadsafencatt.h"
#include "timeseries.h"

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
                readSuccess = readInputFileBoundaryConditionsTag(line, error);
                break;
              case 6:
                readSuccess = readInputFileUniformHydraulicsTag(line, error);
                break;
              case 7:
                readSuccess = readInputFileNonUniformHydraulicsTag(line, error);
                break;
              case 8:
                readSuccess = readInputFileUniformRadiativeFluxesTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileNonUniformRadiativeFluxesTag(line, error);
                break;
              case 10:
                readSuccess = readInputFileUniformMeteorologyTag(line, error);
                break;
              case 11:
                readSuccess = readInputFileNonUniformMeteorologyTag(line, error);
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

    ThreadSafeNcVar junctionX =  m_outputNetCDF->addVar("x", NcType::nc_DOUBLE, junctionDim);
    junctionX.putAtt("long_name", "Junction X-Coordinate");
    junctionX.putAtt("units", "m");
    m_outNetCDFVariables["x"] = junctionX;

    ThreadSafeNcVar junctionY =  m_outputNetCDF->addVar("y", NcType::nc_DOUBLE, junctionDim);
    junctionY.putAtt("long_name", "Junction Y-Zoordinate");
    junctionY.putAtt("units", "m");
    m_outNetCDFVariables["y"] = junctionY;

    ThreadSafeNcVar junctionZ =  m_outputNetCDF->addVar("z", NcType::nc_DOUBLE, junctionDim);
    junctionZ.putAtt("long_name", "junction z-coordinate");
    junctionZ.putAtt("units", "m");
    m_outNetCDFVariables["z"] = junctionZ;

    double *vertx = new double[m_elementJunctions.size()];
    double *verty = new double[m_elementJunctions.size()];
    double *vertz = new double[m_elementJunctions.size()];
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

    //    ThreadSafeNcVar elementsVar =  m_outputNetCDF->addVar("elements", NcType::nc_DOUBLE, elementsDim);
    //    elementsVar.putAtt("long_name", "Distance");
    //    elementsVar.putAtt("units", "m");
    //    m_outNetCDFVariables["elements"] = elementsVar;

    int *fromJunctions = new int[m_elements.size()];
    int *toJunctions = new int[m_elements.size()];
    char **elementIds = new char *[m_elements.size()];
    //    double *els = new double[m_elements.size()];

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

    //hydraulics variables
    ThreadSafeNcVar depthVar =  m_outputNetCDF->addVar("depth", "double",
                                                       std::vector<std::string>({"time", "elements"}));
    depthVar.putAtt("long_name", "Flow Depth");
    depthVar.putAtt("units", "m");
    m_outNetCDFVariables["depth"] = depthVar;

    ThreadSafeNcVar widthVar =  m_outputNetCDF->addVar("width", "double",
                                                       std::vector<std::string>({"time", "elements"}));
    widthVar.putAtt("long_name", "Flow Top Width");
    widthVar.putAtt("units", "m");
    m_outNetCDFVariables["width"] = widthVar;

    ThreadSafeNcVar temperatureVar =  m_outputNetCDF->addVar("channel_temperature", "double",
                                                             std::vector<std::string>({"time", "elements"}));
    temperatureVar.putAtt("long_name", "Channel Temperature");
    temperatureVar.putAtt("units", "°C");
    m_outNetCDFVariables["channel_temperature"] = temperatureVar;


    ThreadSafeNcVar incomingSWSolarRadiationVar =  m_outputNetCDF->addVar("incoming_shortwave_solar_radiation", "double",
                                                                          std::vector<std::string>({"time", "elements"}));
    incomingSWSolarRadiationVar.putAtt("long_name", "Incoming Shortwave Solar Radiation");
    incomingSWSolarRadiationVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["incoming_shortwave_solar_radiation"] = incomingSWSolarRadiationVar;


    ThreadSafeNcVar netSWSolarRadiationVar =  m_outputNetCDF->addVar("net_shortwave_solar_radiation", "double",
                                                                     std::vector<std::string>({"time", "elements"}));
    netSWSolarRadiationVar.putAtt("long_name", "Net Shortwave Solar Radiation");
    netSWSolarRadiationVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["net_shortwave_solar_radiation"] = netSWSolarRadiationVar;


    ThreadSafeNcVar backLWRadiationVar =  m_outputNetCDF->addVar("back_longwave_radiation", "double",
                                                                 std::vector<std::string>({"time", "elements"}));
    backLWRadiationVar.putAtt("long_name", "Back Longwave Radiation");
    backLWRadiationVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["back_longwave_radiation"] = backLWRadiationVar;

    ThreadSafeNcVar sedNetSWSolarRadiationVar =  m_outputNetCDF->addVar("sediment_net_shortwave_solar_radiation", "double",
                                                                        std::vector<std::string>({"time", "elements"}));
    sedNetSWSolarRadiationVar.putAtt("long_name", "Net Shortwave Solar Radiation Reaching Sediment");
    sedNetSWSolarRadiationVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["sediment_net_shortwave_solar_radiation"] = sedNetSWSolarRadiationVar;


    ThreadSafeNcVar atmosphericLWRadiationVar =  m_outputNetCDF->addVar("atmospheric_longwave_radiation", "double",
                                                                        std::vector<std::string>({"time", "elements"}));
    atmosphericLWRadiationVar.putAtt("long_name", "Atmospheric Longwave Radiation");
    atmosphericLWRadiationVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["atmospheric_longwave_radiation"] = atmosphericLWRadiationVar;


    ThreadSafeNcVar landCoverLWRadiationVar =  m_outputNetCDF->addVar("landcover_longwave_radiation", "double",
                                                                      std::vector<std::string>({"time", "elements"}));
    landCoverLWRadiationVar.putAtt("long_name", "Landcover Longwave Radiation");
    landCoverLWRadiationVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["landcover_longwave_radiation"] = landCoverLWRadiationVar;


    ThreadSafeNcVar sumMCRadiationVar =  m_outputNetCDF->addVar("net_main_channel_radiation", "double",
                                                                std::vector<std::string>({"time", "elements"}));
    sumMCRadiationVar.putAtt("long_name", "Net Radiation Reaching Main Channel");
    sumMCRadiationVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["net_main_channel_radiation"] = sumMCRadiationVar;


    ThreadSafeNcVar shadeVar =  m_outputNetCDF->addVar("shade_factor", "double",
                                                       std::vector<std::string>({"time", "elements"}));
    shadeVar.putAtt("long_name", "Shade Factor");
    shadeVar.putAtt("units", "");
    m_outNetCDFVariables["shade_factor"] = shadeVar;


    ThreadSafeNcVar shadeMultVar =  m_outputNetCDF->addVar("shade_factor_multiplier", "double",
                                                       std::vector<std::string>({"time", "elements"}));
    shadeMultVar.putAtt("long_name", "Shade Factor Multiplier");
    shadeMultVar.putAtt("units", "");
    m_outNetCDFVariables["shade_factor_multiplier"] = shadeMultVar;

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
            m_verbose = QString::compare(options[1], "No", Qt::CaseInsensitive);
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

bool RHEModel::readInputFileBoundaryConditionsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];
    auto it = m_elementsById.find(id.toStdString());

    if (it != m_elementsById.end())
    {
      Element *element = it->second;

      bool found = false;
      QString type = columns[2];
      QString variable = columns[1].trimmed();

      if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
      {
        if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
        {
          bool valueOk ;
          double value = columns[3].toDouble(&valueOk);

          if (valueOk)
          {
            element->channelTemperature = value;
            found = true;
          }
          else
          {
            errorMessage = "Temperature BC value is invalid";
            return false;
          }
        }
        else
        {
          //          for (size_t i = 0; i < m_solutes.size(); i++)
          //          {
          //            std::string solute = m_solutes[i];

          //            if (!solute.compare(variable.toStdString()))
          //            {
          //              bool ok;
          //              double value = columns[3].toDouble(&ok);

          //              if (ok)
          //              {
          //                element->channelTemperature = value;
          //                found = true;
          //              }
          //            }
          //          }
        }

        if (!found)
        {
          return false;
        }
      }
      else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
      {
        int variableIndex = -2;

        if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
        {
          variableIndex = -1;
        }
        else
        {
          //          for (size_t i = 0; i < m_solutes.size(); i++)
          //          {
          //            std::string solute = m_solutes[i];

          //            if (!solute.compare(variable.toStdString()))
          //            {
          //              variableIndex = i;
          //              break;
          //            }
          //          }
        }

        if (variableIndex > -2)
        {
          QString filePath = columns[3];

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                ElementTimeSeriesBC *elementBC = new ElementTimeSeriesBC(element, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  elementBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(elementBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified BC filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified BC filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified BC filepath does not exist";
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

bool RHEModel::readInputFileUniformHydraulicsTag(const QString &line, QString &errorMessage)
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
            UniformHydraulicsTimeSeriesBC *uniformHydraulicsBC = new UniformHydraulicsTimeSeriesBC(fromElement, toElement,
                                                                                                   variableIndex, this);
            uniformHydraulicsBC->addValue(m_startDateTime, value);
            uniformHydraulicsBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(uniformHydraulicsBC);
          }
          else
          {
            errorMessage = "Uniform hydraulics value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "FILE", Qt::CaseInsensitive))
        {
          QString filePath = varValue;

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                UniformHydraulicsTimeSeriesBC *uniformHydraulicsBC = new UniformHydraulicsTimeSeriesBC(fromElement, toElement,
                                                                                                       variableIndex, this);
                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  uniformHydraulicsBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(uniformHydraulicsBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified uniform hydraulics filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified uniform hydraulics filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified uniform hydraulics filepath does not exist";
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

bool RHEModel::readInputFileNonUniformHydraulicsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 2)
  {
    auto it = m_hydraulicVariableFlags.find(columns[0].toStdString());

    if (it != m_hydraulicVariableFlags.end())
    {
      int variableIndex = it->second;

      QString filePath = columns[1];

      if (!filePath.isEmpty() && !filePath.isNull())
      {
        QFileInfo fileInfo(filePath);

        if (fileInfo.isRelative())
          fileInfo = relativePathToAbsolute(fileInfo);

        if (QFile::exists(fileInfo.absoluteFilePath()))
        {
          std::map<double, std::vector<double>> timeSeries;
          std::vector<std::string> headers;

          if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
          {
            for (size_t i = 0; i < headers.size(); i++)
            {
              auto eit = m_elementsById.find(headers[i]);

              if (eit != m_elementsById.end())
              {
                Element *element = eit->second;
                HydraulicsTimeSeriesBC *hydraulicsTimeSeries = new HydraulicsTimeSeriesBC(element, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[i];
                  hydraulicsTimeSeries->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(hydraulicsTimeSeries);
              }
              else
              {
                errorMessage = "Specified time varying hydraulic file is invalid";
                return false;
              }
            }
          }
          else
          {
            errorMessage = "Specified time varying hydraulic file is invalid";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified time varying hydraulic file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying hydraulic file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying hydraulic file is invalid";
      return false;
    }
  }
  else
  {
    errorMessage = "Specified time varying hydraulic file is invalid";
    return false;
  }

  return true;
}

bool RHEModel::readInputFileUniformRadiativeFluxesTag(const QString &line, QString &errorMessage)
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
          UniformRadiativeFluxTimeSeriesBC *uniformRadiationFluxBC = new UniformRadiativeFluxTimeSeriesBC(fromElement, toElement, this);
          uniformRadiationFluxBC->addValue(m_startDateTime, value);
          uniformRadiationFluxBC->addValue(m_endDateTime, value);
          m_boundaryConditions.push_back(uniformRadiationFluxBC);
        }
        else
        {
          errorMessage = "Radiation BC value is invalid";
          return false;
        }
      }
      else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
      {
        QString filePath = varValue;

        if (!filePath.isEmpty() && !filePath.isNull())
        {
          QFileInfo fileInfo(filePath);

          if (fileInfo.isRelative())
            fileInfo = relativePathToAbsolute(fileInfo);

          if (QFile::exists(fileInfo.absoluteFilePath()))
          {
            std::map<double, std::vector<double>> timeSeries;
            std::vector<std::string> headers;

            if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
            {
              UniformRadiativeFluxTimeSeriesBC *uniformRadiationFluxBC = new UniformRadiativeFluxTimeSeriesBC(fromElement, toElement, this);

              for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
              {
                double dateTime = it->first;
                double value = it->second[0];
                uniformRadiationFluxBC->addValue(dateTime, value);
              }

              m_boundaryConditions.push_back(uniformRadiationFluxBC);

              timeSeries.clear();
              headers.clear();
            }
            else
            {
              errorMessage = "Specified BC filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified BC filepath does not exist";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified BC filepath does not exist";
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

bool RHEModel::readInputFileNonUniformRadiativeFluxesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 1)
  {
    QString filePath = columns[0];

    if (!filePath.isEmpty() && !filePath.isNull())
    {
      QFileInfo fileInfo(filePath);

      if (fileInfo.isRelative())
        fileInfo = relativePathToAbsolute(fileInfo);

      if (QFile::exists(fileInfo.absoluteFilePath()))
      {
        std::map<double, std::vector<double>> timeSeries;
        std::vector<std::string> headers;

        if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
        {
          for (size_t i = 0; i < headers.size(); i++)
          {
            auto eit = m_elementsById.find(headers[i]);

            if (eit != m_elementsById.end())
            {
              Element *element = eit->second;
              RadiativeFluxTimeSeriesBC *radiativeFluxTimeSeries = new RadiativeFluxTimeSeriesBC(element, this);

              for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
              {
                double dateTime = it->first;
                double value = it->second[i];
                radiativeFluxTimeSeries->addValue(dateTime, value);
              }

              m_boundaryConditions.push_back(radiativeFluxTimeSeries);
            }
            else
            {
              errorMessage = "Specified time varying radiative flux file is invalid";
              return false;
            }
          }
        }
        else
        {
          errorMessage = "Specified time varying radiative flux file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying radiative flux file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying radiative flux file is invalid";
      return false;
    }

  }
  else
  {
    errorMessage = "Specified time varying radiative flux file is invalid";
    return false;
  }

  return true;
}

bool RHEModel::readInputFileUniformMeteorologyTag(const QString &line, QString &errorMessage)
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
            UniformMeteorologyTimeSeriesBC *uniformMeteorologyBC = new UniformMeteorologyTimeSeriesBC(fromElement, toElement,
                                                                                                      variableIndex, this);
            uniformMeteorologyBC->addValue(m_startDateTime, value);
            uniformMeteorologyBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(uniformMeteorologyBC);
          }
          else
          {
            errorMessage = "Uniform meteorology value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "FILE", Qt::CaseInsensitive))
        {
          QString filePath = varValue;

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                UniformMeteorologyTimeSeriesBC *uniformMeteorologyBC = new UniformMeteorologyTimeSeriesBC(fromElement, toElement,
                                                                                                          variableIndex, this);
                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  uniformMeteorologyBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(uniformMeteorologyBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified uniform meteorology filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified uniform meteorology filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified uniform meteorology filepath does not exist";
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

bool RHEModel::readInputFileNonUniformMeteorologyTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 2)
  {
    auto it = m_meteorologicalVariableFlags.find(columns[0].toStdString());

    if (it != m_meteorologicalVariableFlags.end())
    {
      int variableIndex = it->second;

      QString filePath = columns[1];

      if (!filePath.isEmpty() && !filePath.isNull())
      {
        QFileInfo fileInfo(filePath);

        if (fileInfo.isRelative())
          fileInfo = relativePathToAbsolute(fileInfo);

        if (QFile::exists(fileInfo.absoluteFilePath()))
        {
          std::map<double, std::vector<double>> timeSeries;
          std::vector<std::string> headers;

          if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
          {
            for (size_t i = 0; i < headers.size(); i++)
            {
              auto eit = m_elementsById.find(headers[i]);

              if (eit != m_elementsById.end())
              {
                Element *element = eit->second;
                MeteorologyTimeSeriesBC *meteorologyTimeSeries = new MeteorologyTimeSeriesBC(element, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[i];
                  meteorologyTimeSeries->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(meteorologyTimeSeries);
              }
              else
              {
                errorMessage = "Specified time varying hydraulic file is invalid";
                return false;
              }
            }
          }
          else
          {
            errorMessage = "Specified time varying hydraulic file is invalid";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified time varying hydraulic file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying hydraulic file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying hydraulic file is invalid";
      return false;
    }
  }
  else
  {
    errorMessage = "Specified time varying hydraulic file is invalid";
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

    double *depth = new double[m_elements.size()];
    double *width = new double[m_elements.size()];
    double *temperature = new double[m_elements.size()];
    double *incomingSWSolarRadiation = new double[m_elements.size()];
    double *netSWSolarRadiation = new double[m_elements.size()];
    double *backwaterLWRadiation = new double[m_elements.size()];
    double *sedNetSWSolarRadiation = new double[m_elements.size()];
    double *atmosphericLWRadiation = new double[m_elements.size()];
    double *landCoverLWRadiation = new double[m_elements.size()];
    double *sumMCRadiation = new double[m_elements.size()];
    double *shadeFactor = new double[m_elements.size()];
    double *shadeFactorMult = new double[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];
      depth[i] = element->channelDepth;
      width[i] = element->channelWidth;
      temperature[i] = element->channelTemperature;
      incomingSWSolarRadiation[i] = element->incomingSWSolarRadiation;
      netSWSolarRadiation[i] = element->netSWSolarRadiation;
      backwaterLWRadiation[i] = element->backLWRadiation;
      sedNetSWSolarRadiation[i] = element->sedNetSWSolarRadiation;
      atmosphericLWRadiation[i] = element->atmosphericLWRadiation;
      landCoverLWRadiation[i] = element->landCoverLWRadiation;
      sumMCRadiation[i] = element->netMCRadiation;
      shadeFactor[i] = element->shadeFactor;
      shadeFactorMult[i] = element->shadeFactorMultiplier;
    }

    m_outNetCDFVariables["depth"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), depth);

    m_outNetCDFVariables["width"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), width);

    m_outNetCDFVariables["channel_temperature"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), temperature);

    m_outNetCDFVariables["net_main_channel_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), sumMCRadiation);

    m_outNetCDFVariables["incoming_shortwave_solar_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), incomingSWSolarRadiation);

    m_outNetCDFVariables["net_shortwave_solar_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), netSWSolarRadiation);

    m_outNetCDFVariables["back_longwave_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), backwaterLWRadiation);

    m_outNetCDFVariables["sediment_net_shortwave_solar_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), sedNetSWSolarRadiation);

    m_outNetCDFVariables["atmospheric_longwave_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), atmosphericLWRadiation);

    m_outNetCDFVariables["landcover_longwave_radiation"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), landCoverLWRadiation);

    m_outNetCDFVariables["shade_factor"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), shadeFactor);

    m_outNetCDFVariables["shade_factor_multiplier"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), shadeFactorMult);

    delete[] depth;
    delete[] width;
    delete[] temperature;
    delete[] incomingSWSolarRadiation;
    delete[] netSWSolarRadiation;
    delete[] backwaterLWRadiation;
    delete[] sedNetSWSolarRadiation;
    delete[] atmosphericLWRadiation;
    delete[] landCoverLWRadiation;
    delete[] sumMCRadiation;
    delete[] shadeFactor;
    delete[] shadeFactorMult;

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
                                                              {"[UNIFORM_HYDRAULICS]", 6},
                                                              {"[NON_UNIFORM_HYDRAULICS]", 7},
                                                              {"[UNIFORM_RADIATIVE_FLUXES]", 8},
                                                              {"[NON_UNIFORM_RADIATIVE_FLUXES]", 9},
                                                              {"[UNIFORM_METEOROLOGY]", 10},
                                                              {"[NON_UNIFORM_METEOROLOGY]", 11}
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
                                                                         });

const QRegExp RHEModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");
