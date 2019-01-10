#ifndef METEOROLOGYTIMESERIESBC_H
#define METEOROLOGYTIMESERIESBC_H

#include "iboundarycondition.h"
#include "rhecomponent_global.h"

#include <QObject>
#include <QSharedPointer>

struct Element;
class TimeSeries;
class RHEModel;
class DataCursor;

class RHECOMPONENT_EXPORT MeteorologyBC: public QObject,
    public virtual IBoundaryCondition
{
  public:

    MeteorologyBC(Element *startElement, Element *endElement, int variableIndex, RHEModel *model);

    virtual ~MeteorologyBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

    QSharedPointer<TimeSeries> timeSeries() const;

    void setTimeSeries(const QSharedPointer<TimeSeries> &timeseries);

  private:
    std::vector<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    int m_variableIndex;
    DataCursor *m_dataCursor;
    QSharedPointer<TimeSeries> m_timeSeries;
    RHEModel *m_model;
};

#endif // METEOROLOGYTIMESERIESBC_H
