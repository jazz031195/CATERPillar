#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "slidergroup.h"

SlidersGroup::SlidersGroup(const QString &title, QWidget *parent)
    : QGroupBox(title, parent)
{
    dial = new QDial;
    dial->setFocusPolicy(Qt::StrongFocus);

}


void SlidersGroup::setMinimum(int value)
{
    dial->setMinimum(value);
}

void SlidersGroup::setMaximum(int value)
{
    dial->setMaximum(value);
}

void SlidersGroup::invertAppearance(bool invert)
{
    dial->setInvertedAppearance(invert);
}

void SlidersGroup::invertKeyBindings(bool invert)
{
    dial->setInvertedControls(invert);
}


