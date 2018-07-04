/*!
*  \file    rheproject.h
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
*  \todo
*  \warning
*/

#ifndef RHEMODEL_H
#define RHEMODEL_H

#include "rhecomponent_global.h"
#include "spatial/network.h"
#include "odesolver.h"

#include <vector>
#include <string>
#include <set>
#include <QFileInfo>
#include <QTextStream>
#include <unordered_map>

#ifdef USE_NETCDF
#include <netcdf>
#endif

class RHEComponent;
struct Element;
struct ElementJunction;
class Edge;
class RHEModel;
class IBoundaryCondition;
class ThreadSafeNcFile;

struct RHECOMPONENT_EXPORT SolverUserData
{
    RHEModel *model = nullptr;
    int variableIndex = -1;
};

typedef void (*RetrieveCouplingData)(RHEModel *model, double dateTime);

class RHECOMPONENT_EXPORT RHEModel : public QObject
{
    Q_OBJECT

    friend struct ElementJunction;
    friend struct Element;
    friend class UniformRadiativeFluxTimeSeriesBC;
    friend class NonPointSrcTimeSeriesBC;
    friend class UniformHydraulicsTimeSeriesBC;
    friend class UniformMeteorologyTimeSeriesBC;

  public:

    /*!
     * \brief RHEModel - Constructor for the Computational engine for the Stream Temperature Model.
     */
    RHEModel(RHEComponent *component);

    /*!
     * \brief ~RHEModel - Destructor for the Computational engine for the Stream Temperature Model.
     */
    ~RHEModel();

    /*!
     * \brief maxTimeStep
     * \return
     */
    double maxTimeStep() const;

    /*!
     * \brief setMaxTimeStep
     * \param timeStep
     */
    void setMaxTimeStep(double timeStep);

    /*!
     * \brief setTimeStepRelaxationFactor - Adaptive time-step relaxation factor
     * \param tStepRelaxFactor
     */
    void setTimeStepRelaxationFactor(double tStepRelaxFactor);

    /*!
     * \brief currentTimeStep - Current time step in seconds
     * \return
     */
    double currentTimeStep() const;

    /*!
     * \brief startDateTime - Start date and time specified as a modified julian date
     * \return
     */
    double startDateTime() const;

    /*!
     * \brief setStartDate - Sets the value of the start date and time
     * \param dateTime - Start date and time specified as modified julian date
     */
    void setStartDateTime(double dateTime);

    /*!
     * \brief endDateTime - End datetime specified as
     * \return
     */
    double endDateTime() const;

    /*!
     * \brief setEndDateTime Sets the value of the end datetime
     * \param dateTime
     */
    void setEndDateTime(double dateTime);

    /*!
     * \brief outputInterval
     * \return
     */
    double outputInterval() const;

    /*!
     * \brief setOutputInterval
     * \param interval
     */
    void setOutputInterval(double interval);

    /*!
     * \brief currentDateTime
     * \return
     */
    double currentDateTime() const;

    /*!
     * \brief waterDensity
     * \return
     */
    double waterDensity() const;

    /*!
     * \brief setWaterDensity
     * \param value
     */
    void setWaterDensity(double value);

    /*!
     * \brief specificHeatCapacityWater
     * \return
     */
    double specificHeatCapacityWater() const;

    /*!
     * \brief setSpecificHeatCapacityWater
     * \param value
     */
    void setSpecificHeatCapacityWater(double value);

    /*!
     * \brief verbose
     * \return
     */
    bool verbose() const;

    /*!
     * \brief setVerbose
     * \param verbose
     */
    void setVerbose(bool verbose);

    /*!
     * \brief printFrequency
     * \return
     */
    int printFrequency() const;

    /*!
     * \brief setPrintFrequency
     * \param printFreq
     */
    void setPrintFrequency(int printFreq);

    /*!
     * \brief flushToDiskFrequency
     * \return
     */
    int flushToDiskFrequency() const;

    /*!
     * \brief setFlushToDiskFrequency
     * \param diskFlushFrequency
     */
    void setFlushToDiskFrequency(int diskFlushFrequency);

    /*!
     * \brief numElementJunctions
     * \return
     */
    int numElementJunctions() const;

    /*!
     * \brief addControlVolumeNode
     * \param id
     * \param x
     * \param y
     * \param z
     * \return
     */
    ElementJunction *addElementJunction(const std::string &id, double x = 0, double y = 0, double z = 0);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(const std::string &id);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(int id);

    /*!
     * \brief getElementJunction
     * \param id
     * \return
     */
    ElementJunction *getElementJunction(const std::string &id) ;

    /*!
     * \brief getElementJunction
     * \param index
     * \return
     */
    ElementJunction *getElementJunction(int index) ;

    /*!
     * \brief numElements
     * \return
     */
    int numElements() const;

    /*!
     * \brief addElement
     * \param id
     * \param fromElement
     * \param toElement
     * \return
     */
    Element *addElement(const std::string &id, ElementJunction *upStream, ElementJunction *downStream);

    /*!
     * \brief removeElement
     * \param id
     */
    void deleteElement(const std::string &id);

    /*!
     * \brief removeElement
     * \param index
     */
    void deleteElement(int index);

    /*!
     * \brief getElement
     * \param id
     * \return
     */
    Element *getElement(const std::string &id);

    /*!
     * \brief getElement
     * \param index
     * \return
     */
    Element *getElement(int index);

    /*!
     * \brief inputFile
     * \return
     */
    QFileInfo inputFile() const;

    /*!
     * \brief setInputFile
     * \param inputFile
     */
    void setInputFile(const QFileInfo &inputFile);

    /*!
     * \brief outputCSVFile
     * \return
     */
    QFileInfo outputCSVFile() const;

    /*!
     * \brief setOutputCSVFile
     * \param outputFile
     */
    void setOutputCSVFile(const QFileInfo &outputFile);

    /*!
     * \brief outputNetCDFFile
     * \return
     */
    QFileInfo outputNetCDFFile() const;

    /*!
     * \brief setOutputNetCDFFile
     * \param outputNetCDFFile
     */
    void setOutputNetCDFFile(const QFileInfo &outputNetCDFFile);

    /*!
     * \brief retrieveCouplingDataFunction
     * \return
     */
    RetrieveCouplingData retrieveCouplingDataFunction() const;

    /*!
     * \brief setRetrieveCouplingDataFunction
     */
    void setRetrieveCouplingDataFunction(RetrieveCouplingData retrieveCouplingDataFunction);

    /*!
     * \brief initialize
     * \param errors
     * \return
     */
    bool initialize(std::list<std::string> &errors);

    /*!
     * \brief update
     */
    void update();

    /*!
     * \brief finalize
     * \param errors
     * \return
     */
    bool finalize(std::list<std::string> &errors);

    /*!
     * \brief printStatus
     */
    void printStatus();

    /*!
     * \brief saveAs
     * \param filePath
     */
    void saveAs(const QFileInfo &filePath);

  private:

    /*!
     * \brief initializeInputFiles
     * \param errors
     * \return
     */
    bool initializeInputFiles(std::list<std::string> &errors);

    /*!
     * \brief initializeTimeVariables
     * \param errors
     * \return
     */
    bool initializeTimeVariables(std::list<std::string> &errors);

    /*!
     * \brief initializeElements
     * \param errors
     * \return
     */
    bool initializeElements(std::list<std::string> &errors);

    /*!
     * \brief initializeSolver
     * \param errors
     * \return
     */
    bool initializeSolver(std::list<std::string> & errors);

    /*!
     * \brief intializeOutputFiles
     * \param errors
     * \return
     */
    bool initializeOutputFiles(std::list<std::string> &errors);

    /*!
     * \brief initializeCSVOutputFile
     * \param errors
     * \return
     */
    bool initializeCSVOutputFile(std::list<std::string> &errors);

    /*!
     * \brief initializeNetCDFOutputFile
     * \param errors
     * \return
     */
    bool initializeNetCDFOutputFile(std::list<std::string> &errors);

    /*!
     * \brief initializeBoundaryConditions
     * \param errors
     * \return
     */
    bool initializeBoundaryConditions(std::list<std::string> &errors);

    /*!
     * \brief prepareForNextTimeStep
     */
    void prepareForNextTimeStep();

    /*!
     * \brief applyInitialConditions
     */
    void applyInitialConditions();

    /*!
     * \brief applyBoundaryConditions
     */
    void applyBoundaryConditions(double dateTime);

    /*!
     * \brief computeTimeStep
     * \return
     */
    double computeTimeStep();

    /*!
     * \brief computeLongDispersion
     */
    void computeLongDispersion();

    /*!
     * \brief solveHeat
     * \param timeStep
     */
    void solveHeatTransport(double timeStep);

    /*!
     * \brief computeDTDt
     * \param model
     * \param variableIndex
     * \param t
     * \param y
     * \param dydt
     */
    static void computeDTDt(double t, double y[], double dydt[], void *userData);

    /*!
     * \brief computeSoluteDYDt
     * \param t
     * \param y
     * \param dydt
     * \param userData
     */
    static void computeDSoluteDt(double t, double y[], double dydt[], void *userData);

    /*!
     * \brief solveJunctionHeatContinuity Solve
     * \param timeStep
     */
    void solveJunctionHeatContinuity(double timeStep);

    /*!
     * \brief readInputFileOptionTag
     * \param line
     */
    bool readInputFileOptionTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileOutputTag
     * \param line
     */
    bool readInputFileOutputTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileSolutesTag
     * \param line
     */
    bool readInputFileSolutesTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementJunctionsTag
     * \param line
     */
    bool readInputFileElementJunctionsTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementsTag
     * \param line
     */
    bool readInputFileElementsTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileBoundaryConditionsTag
     * \param line
     */
    bool readInputFileBoundaryConditionsTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileUniformHydraulicsTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileUniformHydraulicsTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileNonUniformHydraulicsTag
     * \param line
     */
    bool readInputFileNonUniformHydraulicsTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileRadiativeFluxesTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileUniformRadiativeFluxesTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileNonUniformRadiativeFluxesTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileNonUniformRadiativeFluxesTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileUniformMeteorologyTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileUniformMeteorologyTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileNonUniformMeteorologyTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileNonUniformMeteorologyTag(const QString &line, QString &errorMessage);

    /*!
     * \brief writeOutput
     */
    void writeOutput();

    /*!
     * \brief writeCSVOutput
     */
    void writeCSVOutput();

    /*!
     * \brief writeNetCDFOutput
     */
    void writeNetCDFOutput();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputFiles();

    /*!
     * \brief closeCSVOutputFile
     */
    void closeCSVOutputFile();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputNetCDFFile();

    /*!
     * \brief relativePathtoAbsolute
     * \param inputFile
     * \return
     */
    QFileInfo relativePathToAbsolute(const QFileInfo& fileInfo);


    /*!
     * \brief findProfile
     * \param from
     * \param to
     * \param m_profile
     * \return
     */
    bool findProfile(Element *from, Element *to, std::list<Element*> &profile);

  private:

    //Time variables
    double m_timeStep, //seconds
    m_startDateTime, //Modified Julian Day
    m_endDateTime, //Modified Julian Day
    m_currentDateTime, //Modified Julian Day
    m_prevDateTime, //Modified Julian Day
    m_maxTimeStep, //seconds
    m_outputInterval, //seconds
    m_nextOutputTime;//Julian Day

    int  m_printFrequency, //Number of timesteps before printing
    m_currentPrintCount, // Number of timesteps
    m_flushToDiskFrequency, // Number of times to write output stored in memory to disk
    m_currentflushToDiskCount; //Number of timesteps that have been stored in memory so far since the last flush to disk


    bool m_verbose, //Print simulation information to console
    m_flushToDisk;

    //Element junctions
    std::vector<ElementJunction*> m_elementJunctions;
    std::unordered_map<std::string, ElementJunction*> m_elementJunctionsById; //added for fast lookup using identifiers instead of indexes.

    //1D Computational elements
    std::vector<Element*> m_elements;
    std::unordered_map<std::string, Element*> m_elementsById; //added for fast lookup using identifiers instead of indexes.

    //Boundary conditions list
    std::vector<IBoundaryCondition*> m_boundaryConditions;

    //Global water properties
    double m_waterDensity, //kg/m^3
    m_cp,// 4187.0; // J/kg/C
    m_extinctionCoefficient, //extinction coefficient m-1
    m_albedo, //
    m_atmLWReflection,
    m_emissWater,
    m_stefanBoltzmannConst,
    m_atmEmissCoeff,
    m_totalIncomingSolarRadiation, //(kJ/m^2)
    m_totalNetSWSolarRadiation, //(kJ/m^2)
    m_totalBackLWRadiation, //(kJ/m^2)
    m_totalSedNetSWSolarRadiation, //(kJ/m^2)
    m_totalAtmosphericLWRadiation, //(kJ/m^2)
    m_totalLandCoverLWRadiation;//(kJ/m^2)

    //File input and output
    QFileInfo m_inputFile, //Input filepath
    m_outputCSVFileInfo, //Output CSV filepath
    m_outputNetCDFFileInfo; //Output NetCDF filepath

#ifdef USE_NETCDF
    ThreadSafeNcFile *m_outputNetCDF = nullptr; //NetCDF output file object
#endif

    QTextStream m_outputCSVStream; //Output CSV filestream
    static const std::unordered_map<std::string, int> m_inputFileFlags, //Input file flags
    m_optionsFlags, //Input file flags
    m_hydraulicVariableFlags, //Hydraulic variable flags
    m_meteorologicalVariableFlags; //Meteorology variables

    static const QRegExp m_dateTimeDelim;
    QRegExp m_delimiters; //Regex delimiter

    //Function to retrieve and apply coupling data
    RetrieveCouplingData m_retrieveCouplingDataFunction;

    //Parent component
    RHEComponent *m_component;
};

#endif // RHEMODEL_H
