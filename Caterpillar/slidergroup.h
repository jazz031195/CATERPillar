#ifndef SLIDERGROUP_H
#define SLIDERGROUP_H

#include <QMainWindow>
#include <QSlider>
#include <QLabel>
#include <QString>
#include <QGroupBox>
#include <QScrollBar>
#include <QDial>
#include <QBoxLayout>


class SlidersGroup : public QGroupBox
{
    Q_OBJECT

public:
    SlidersGroup(const QString &title, QWidget *parent = nullptr);

signals:
    void valueChanged(int value);

public slots:

    void setMinimum(int value);
    void setMaximum(int value);
    void invertAppearance(bool invert);
    void invertKeyBindings(bool invert);

private:
    QDial *dial;
    QBoxLayout *slidersLayout;
};

#endif // SLIDERGROUP_H
