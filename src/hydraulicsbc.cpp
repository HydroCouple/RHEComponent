
#include "stdafx.h"
#include "element.h"
#include "rhemodel.h"
#include "hydraulicsbc.h"
#include "core/datacursor.h"
#include "temporal/timeseries.h"

HydraulicsBC::HydraulicsBC(Element *startElement,
                           Element *endElement,
                           int variableIndex,
                           RHEModel *model)
  : QObject(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_variableIndex(variableIndex),
    m_model(model)
{
  m_dataCursor = new DataCursor();
}

HydraulicsBC::~HydraulicsBC()
{
  delete m_dataCursor;
}

void HydraulicsBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void HydraulicsBC::prepare()
{

}

void HydraulicsBC::applyBoundaryConditions(double dateTime)
{
  double value = 0;

  if(m_timeSeries->numColumns() == (int) m_profile.size())
  {
    switch (m_variableIndex)
    {
      case 1:
        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
          {
            m_profile[i]->channelDepth = value;
          }
        }
        break;
      case 2:
        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
          {
            m_profile[i]->channelWidth = value;
          }
        }
        break;
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      switch (m_variableIndex)
      {
        case 1:
          for(size_t i = 0 ; i < m_profile.size(); i++)
          {
            m_profile[i]->channelDepth = value;
          }
          break;
        case 2:
          for(size_t i = 0 ; i < m_profile.size(); i++)
          {
            m_profile[i]->channelWidth = value;
          }
          break;
      }
    }
  }
}

void HydraulicsBC::clear()
{

}

Element *HydraulicsBC::startElement() const
{
  return m_startElement;
}

void HydraulicsBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *HydraulicsBC::endElement() const
{
  return m_endElement;
}

void HydraulicsBC::setEndElement(Element *element)
{
  m_endElement = element;
}

QSharedPointer<TimeSeries> HydraulicsBC::timeSeries() const
{
  return m_timeSeries;
}

void HydraulicsBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}
