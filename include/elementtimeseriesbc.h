#ifndef ELEMENTTIMESERIESBC_H
#define ELEMENTTIMESERIESBC_H

#include "abstracttimeseriesbc.h"


class RHECOMPONENT_EXPORT ElementTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    ElementTimeSeriesBC(Element *element, int variableIndex, RHEModel *model);

    virtual ~ElementTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *element() const;

    void setElement(Element *element);

  private:

    Element *m_element;
    int m_variableIndex;
};


#endif // ELEMENTTIMESERIESBC_H
