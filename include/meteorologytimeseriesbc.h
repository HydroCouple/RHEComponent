#ifndef METEOROLOGYTIMESERIESBC_H
#define METEOROLOGYTIMESERIESBC_H


#include "abstracttimeseriesbc.h"

class RHECOMPONENT_EXPORT MeteorologyTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    MeteorologyTimeSeriesBC(Element *element, int variableIndex, RHEModel *model);

    virtual ~MeteorologyTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *element() const;

    void setElement(Element *element);

  private:

    Element *m_element;
    int m_variableIndex;
};


class RHECOMPONENT_EXPORT UniformMeteorologyTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    UniformMeteorologyTimeSeriesBC(Element *startElement, Element *endElement, int variableIndex, RHEModel *model);

    virtual ~UniformMeteorologyTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

  private:
    std::list<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    int m_variableIndex;
};

#endif // METEOROLOGYTIMESERIESBC_H
