/*!
 *  \file    rhecomponent.h
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

#ifndef RHECOMPONENT_H
#define RHECOMPONENT_H

#include "rhecomponent_global.h"
#include "rhecomponentinfo.h"
#include "temporal/abstracttimemodelcomponent.h"

#include <unordered_map>

class RHEModel;
class HCGeometry;
struct Element;
class ElementInput;
class ElementOutput;
class Dimension;
class Unit;

class RHECOMPONENT_EXPORT RHEComponent : public AbstractTimeModelComponent,
    public virtual HydroCouple::ICloneableModelComponent
{
    Q_OBJECT

    Q_INTERFACES(HydroCouple::ICloneableModelComponent)


  public:

    /*!
     * \brief RHEComponent constructor
     * \param id Unique identifier for this component instance.
     * \param modelComponentInfo the parent ModelComponentInfo that generated this component instance.
     */
    RHEComponent(const QString &id, RHEComponentInfo* modelComponentInfo = nullptr);

    /*!
     * \brief ~RHEComponent destructor
     */
    virtual ~RHEComponent() override;

    /*!
     * \brief validate validates this component model instance
     * \return Returns a list of error messages.
     */
    QList<QString> validate() override;

    /*!
     * \brief prepare Prepares the model component instance.
     */
    void prepare() override;

    /*!
     * \brief update
     * \param requiredOutputs
     */
    void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

    /*!
     * \brief finish
     */
    void finish() override;

    /*!
     * \brief modelInstance
     * \return
     */
    RHEModel *modelInstance() const;

    /*!
     * \brief parent
     * \return
     */
    HydroCouple::ICloneableModelComponent* parent() const override;

    /*!
     * \brief clone
     * \return
     */
    HydroCouple::ICloneableModelComponent* clone() override;

    /*!
     * \brief clones
     * \return
     */
    QList<HydroCouple::ICloneableModelComponent*> clones() const override;


  protected:

    bool removeClone(RHEComponent *component);

    /*!
     * \brief intializeFailureCleanUp
     */
    void initializeFailureCleanUp() override;

  private:

    /*!
     * \brief createArguments
     */
    void createArguments() override;

    /*!
     * \brief createInputFileArguments
     */
    void createInputFileArguments();

    /*!
     * \brief initializeArguments
     * \param message
     * \return
     */
    bool initializeArguments(QString &message) override;

    /*!
     * \brief initializeInputFilesArguments
     * \param message
     * \return
     */
    bool initializeInputFilesArguments(QString &message);

    /*!
     * \brief createGeometriesMap
     */
    void createGeometries();

    /*!
     * \brief createInputs
     */
    void createInputs() override;

    /*!
     * \brief createDepthInput
     */
    void createDepthInput();

    /*!
     * \brief createTopWidthInput
     */
    void createTopWidthInput();

    /*!
     * \brief createTemperatureInput
     */
    void createTemperatureInput();

    /*!
     * \brief createLandCoverTemperatureInput
     */
    void createLandCoverTemperatureInput();

    /*!
     * \brief createSkyViewFactorInput
     */
    void createSkyViewFactorInput();

    /*!
     * \brief createShadeFactorInput
     */
    void createShadeFactorInput();

    /*!
     * \brief createShadeFactorMultiplierInput
     */
    void createShadeFactorMultiplierInput();

    /*!
     * \brief createOutputs
     */
    void createOutputs() override;

    /*!
     * \brief createRadiationOutputs
     */
    void createRadiationOutputs();

  private:

    IdBasedArgumentString *m_inputFilesArgument;

    ElementInput *m_depthInput,
    *m_topWidthInput,
    *m_channelTemperatureInput,
    *m_skyviewFactorInput,
    *m_shadeFactorInput,
    *m_shadeFactorMultiplierInput;


    ElementOutput *m_netSWSolarRadiationOutput,
    *m_backLWRadiationOutput,
    *m_sedNetSWSolarRadiationOutput,
    *m_atmosphericLWRadiationOutput,
    *m_landCoverLWRadiatioOutput;

    Unit *m_radiationFluxUnit,
    *m_heatFluxUnit,
    *m_temperatureUnit;

    Dimension *m_timeDimension,
    *m_geometryDimension;

    std::vector<QSharedPointer<HCGeometry>> m_elementGeometries;
    std::vector<QSharedPointer<HCGeometry>> m_elementJunctionGeometries;
    RHEModel *m_modelInstance;

    RHEComponent *m_parent;
    QList<HydroCouple::ICloneableModelComponent*> m_clones;
};

#endif //RHECOMPONENT_H
