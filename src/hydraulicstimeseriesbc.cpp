
#include "element.h"
#include "rhemodel.h"
#include "hydraulicstimeseriesbc.h"

HydraulicsTimeSeriesBC::HydraulicsTimeSeriesBC(Element *element, int variableIndex, RHEModel *model)
  :AbstractTimeSeriesBC(model),
   m_element(element),
   m_variableIndex(variableIndex)
{

}

HydraulicsTimeSeriesBC::~HydraulicsTimeSeriesBC()
{

}

void HydraulicsTimeSeriesBC::findAssociatedGeometries()
{

}

void HydraulicsTimeSeriesBC::prepare()
{

}

void HydraulicsTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableIndex)
    {
      case 1:
        m_element->channelDepth = value;
        break;
      case 2:
        m_element->channelWidth = value;
        break;
    }
  }
}

Element *HydraulicsTimeSeriesBC::element() const
{
  return m_element;
}

void HydraulicsTimeSeriesBC::setElement(Element *element)
{
  m_element = element;
}


UniformHydraulicsTimeSeriesBC::UniformHydraulicsTimeSeriesBC(Element *startElement, Element *endElement, int variableIndex, RHEModel *model)
  :AbstractTimeSeriesBC(model),
   m_startElement(startElement),
   m_endElement(endElement),
   m_variableIndex(variableIndex)
{

}

UniformHydraulicsTimeSeriesBC::~UniformHydraulicsTimeSeriesBC()
{

}

void UniformHydraulicsTimeSeriesBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void UniformHydraulicsTimeSeriesBC::prepare()
{

}

void UniformHydraulicsTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    for(Element *element : m_profile)
    {
      switch (m_variableIndex)
      {
        case 1:
          element->channelDepth = value;
          break;
        case 2:
          element->channelWidth = value;
          break;
      }
    }
  }
}

Element *UniformHydraulicsTimeSeriesBC::startElement() const
{
  return m_startElement;
}

void UniformHydraulicsTimeSeriesBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *UniformHydraulicsTimeSeriesBC::endElement() const
{
  return m_endElement;
}

void UniformHydraulicsTimeSeriesBC::setEndElement(Element *element)
{
  m_endElement = element;
}
