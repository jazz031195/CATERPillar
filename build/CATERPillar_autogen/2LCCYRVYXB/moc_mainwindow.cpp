/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../GUI/mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Window_t {
    QByteArrayData data[28];
    char stringdata0[458];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Window_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Window_t qt_meta_stringdata_Window = {
    {
QT_MOC_LITERAL(0, 0, 6), // "Window"
QT_MOC_LITERAL(1, 7, 19), // "onSaveButtonClicked"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 30), // "onSelectDirectoryButtonClicked"
QT_MOC_LITERAL(4, 59, 19), // "SelectSWCFileButton"
QT_MOC_LITERAL(5, 79, 9), // "PlotCells"
QT_MOC_LITERAL(6, 89, 10), // "axons_plot"
QT_MOC_LITERAL(7, 100, 15), // "glial_pop1_plot"
QT_MOC_LITERAL(8, 116, 15), // "glial_pop2_plot"
QT_MOC_LITERAL(9, 132, 18), // "blood_vessels_plot"
QT_MOC_LITERAL(10, 151, 21), // "ReadGlialCellsFromSWC"
QT_MOC_LITERAL(11, 173, 8), // "filePath"
QT_MOC_LITERAL(12, 182, 16), // "ReadAxonsFromSWC"
QT_MOC_LITERAL(13, 199, 16), // "ReadAxonsFromCSV"
QT_MOC_LITERAL(14, 216, 8), // "fileName"
QT_MOC_LITERAL(15, 225, 17), // "ReadAxonsFromFile"
QT_MOC_LITERAL(16, 243, 22), // "ReadGlialCellsFromFile"
QT_MOC_LITERAL(17, 266, 21), // "ReadGlialCellsFromCSV"
QT_MOC_LITERAL(18, 288, 24), // "ReadBloodVesselsFromFile"
QT_MOC_LITERAL(19, 313, 22), // "generateMonteCarloConf"
QT_MOC_LITERAL(20, 336, 15), // "runMCSimulation"
QT_MOC_LITERAL(21, 352, 16), // "onGrowthProgress"
QT_MOC_LITERAL(22, 369, 15), // "completed_depth"
QT_MOC_LITERAL(23, 385, 11), // "total_depth"
QT_MOC_LITERAL(24, 397, 18), // "onSwellingProgress"
QT_MOC_LITERAL(25, 416, 12), // "current_icvf"
QT_MOC_LITERAL(26, 429, 11), // "target_icvf"
QT_MOC_LITERAL(27, 441, 16) // "onGrowthFinished"

    },
    "Window\0onSaveButtonClicked\0\0"
    "onSelectDirectoryButtonClicked\0"
    "SelectSWCFileButton\0PlotCells\0axons_plot\0"
    "glial_pop1_plot\0glial_pop2_plot\0"
    "blood_vessels_plot\0ReadGlialCellsFromSWC\0"
    "filePath\0ReadAxonsFromSWC\0ReadAxonsFromCSV\0"
    "fileName\0ReadAxonsFromFile\0"
    "ReadGlialCellsFromFile\0ReadGlialCellsFromCSV\0"
    "ReadBloodVesselsFromFile\0"
    "generateMonteCarloConf\0runMCSimulation\0"
    "onGrowthProgress\0completed_depth\0"
    "total_depth\0onSwellingProgress\0"
    "current_icvf\0target_icvf\0onGrowthFinished"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Window[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   94,    2, 0x08 /* Private */,
       3,    0,   95,    2, 0x08 /* Private */,
       4,    0,   96,    2, 0x08 /* Private */,
       5,    4,   97,    2, 0x08 /* Private */,
      10,    1,  106,    2, 0x08 /* Private */,
      12,    1,  109,    2, 0x08 /* Private */,
      13,    1,  112,    2, 0x08 /* Private */,
      15,    1,  115,    2, 0x08 /* Private */,
      16,    1,  118,    2, 0x08 /* Private */,
      17,    1,  121,    2, 0x08 /* Private */,
      18,    1,  124,    2, 0x08 /* Private */,
      19,    0,  127,    2, 0x08 /* Private */,
      20,    0,  128,    2, 0x08 /* Private */,
      21,    2,  129,    2, 0x08 /* Private */,
      24,    2,  134,    2, 0x08 /* Private */,
      27,    0,  139,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool, QMetaType::Bool, QMetaType::Bool, QMetaType::Bool,    6,    7,    8,    9,
    QMetaType::Void, QMetaType::QString,   11,
    QMetaType::Void, QMetaType::QString,   11,
    QMetaType::Void, QMetaType::QString,   14,
    QMetaType::Void, QMetaType::QString,   14,
    QMetaType::Void, QMetaType::QString,   14,
    QMetaType::Void, QMetaType::QString,   14,
    QMetaType::Void, QMetaType::QString,   14,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Double, QMetaType::Double,   22,   23,
    QMetaType::Void, QMetaType::Double, QMetaType::Double,   25,   26,
    QMetaType::Void,

       0        // eod
};

void Window::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Window *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->onSaveButtonClicked(); break;
        case 1: _t->onSelectDirectoryButtonClicked(); break;
        case 2: _t->SelectSWCFileButton(); break;
        case 3: _t->PlotCells((*reinterpret_cast< const bool(*)>(_a[1])),(*reinterpret_cast< const bool(*)>(_a[2])),(*reinterpret_cast< const bool(*)>(_a[3])),(*reinterpret_cast< const bool(*)>(_a[4]))); break;
        case 4: _t->ReadGlialCellsFromSWC((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: _t->ReadAxonsFromSWC((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: _t->ReadAxonsFromCSV((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: _t->ReadAxonsFromFile((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 8: _t->ReadGlialCellsFromFile((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 9: _t->ReadGlialCellsFromCSV((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 10: _t->ReadBloodVesselsFromFile((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 11: _t->generateMonteCarloConf(); break;
        case 12: _t->runMCSimulation(); break;
        case 13: _t->onGrowthProgress((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 14: _t->onSwellingProgress((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 15: _t->onGrowthFinished(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject Window::staticMetaObject = { {
    QMetaObject::SuperData::link<QWidget::staticMetaObject>(),
    qt_meta_stringdata_Window.data,
    qt_meta_data_Window,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *Window::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Window::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Window.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int Window::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 16)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 16;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
