/****************************************************************************
** Meta object code from reading C++ file 'ScatterDataModifier.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../Caterpillar/ScatterDataModifier.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ScatterDataModifier.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_ScatterDataModifier_t {
    QByteArrayData data[31];
    char stringdata0[442];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ScatterDataModifier_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ScatterDataModifier_t qt_meta_stringdata_ScatterDataModifier = {
    {
QT_MOC_LITERAL(0, 0, 19), // "ScatterDataModifier"
QT_MOC_LITERAL(1, 20, 24), // "backgroundEnabledChanged"
QT_MOC_LITERAL(2, 45, 0), // ""
QT_MOC_LITERAL(3, 46, 7), // "enabled"
QT_MOC_LITERAL(4, 54, 18), // "gridEnabledChanged"
QT_MOC_LITERAL(5, 73, 11), // "fontChanged"
QT_MOC_LITERAL(6, 85, 4), // "font"
QT_MOC_LITERAL(7, 90, 20), // "shadowQualityChanged"
QT_MOC_LITERAL(8, 111, 7), // "quality"
QT_MOC_LITERAL(9, 119, 13), // "setSmoothDots"
QT_MOC_LITERAL(10, 133, 6), // "smooth"
QT_MOC_LITERAL(11, 140, 11), // "changeTheme"
QT_MOC_LITERAL(12, 152, 5), // "theme"
QT_MOC_LITERAL(13, 158, 18), // "changePresetCamera"
QT_MOC_LITERAL(14, 177, 16), // "changeLabelStyle"
QT_MOC_LITERAL(15, 194, 10), // "changeFont"
QT_MOC_LITERAL(16, 205, 28), // "shadowQualityUpdatedByVisual"
QT_MOC_LITERAL(17, 234, 52), // "QtDataVisualization::QAbstrac..."
QT_MOC_LITERAL(18, 287, 2), // "sq"
QT_MOC_LITERAL(19, 290, 19), // "changeShadowQuality"
QT_MOC_LITERAL(20, 310, 20), // "setBackgroundEnabled"
QT_MOC_LITERAL(21, 331, 14), // "setGridEnabled"
QT_MOC_LITERAL(22, 346, 18), // "setVoxelParameters"
QT_MOC_LITERAL(23, 365, 11), // "voxel_size_"
QT_MOC_LITERAL(24, 377, 33), // "std::vector<std::vector<doubl..."
QT_MOC_LITERAL(25, 411, 2), // "X_"
QT_MOC_LITERAL(26, 414, 2), // "Y_"
QT_MOC_LITERAL(27, 417, 2), // "Z_"
QT_MOC_LITERAL(28, 420, 2), // "R_"
QT_MOC_LITERAL(29, 423, 11), // "modifyScale"
QT_MOC_LITERAL(30, 435, 6) // "scale_"

    },
    "ScatterDataModifier\0backgroundEnabledChanged\0"
    "\0enabled\0gridEnabledChanged\0fontChanged\0"
    "font\0shadowQualityChanged\0quality\0"
    "setSmoothDots\0smooth\0changeTheme\0theme\0"
    "changePresetCamera\0changeLabelStyle\0"
    "changeFont\0shadowQualityUpdatedByVisual\0"
    "QtDataVisualization::QAbstract3DGraph::ShadowQuality\0"
    "sq\0changeShadowQuality\0setBackgroundEnabled\0"
    "setGridEnabled\0setVoxelParameters\0"
    "voxel_size_\0std::vector<std::vector<double> >\0"
    "X_\0Y_\0Z_\0R_\0modifyScale\0scale_"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ScatterDataModifier[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      15,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   89,    2, 0x06 /* Public */,
       4,    1,   92,    2, 0x06 /* Public */,
       5,    1,   95,    2, 0x06 /* Public */,
       7,    1,   98,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       9,    1,  101,    2, 0x0a /* Public */,
      11,    1,  104,    2, 0x0a /* Public */,
      13,    0,  107,    2, 0x0a /* Public */,
      14,    0,  108,    2, 0x0a /* Public */,
      15,    1,  109,    2, 0x0a /* Public */,
      16,    1,  112,    2, 0x0a /* Public */,
      19,    1,  115,    2, 0x0a /* Public */,
      20,    1,  118,    2, 0x0a /* Public */,
      21,    1,  121,    2, 0x0a /* Public */,
      22,    5,  124,    2, 0x0a /* Public */,
      29,    1,  135,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Bool,    3,
    QMetaType::Void, QMetaType::Bool,    3,
    QMetaType::Void, QMetaType::QFont,    6,
    QMetaType::Void, QMetaType::Int,    8,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,   10,
    QMetaType::Void, QMetaType::Int,   12,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QFont,    6,
    QMetaType::Void, 0x80000000 | 17,   18,
    QMetaType::Void, QMetaType::Int,    8,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int, 0x80000000 | 24, 0x80000000 | 24, 0x80000000 | 24, 0x80000000 | 24,   23,   25,   26,   27,   28,
    QMetaType::Void, QMetaType::Double,   30,

       0        // eod
};

void ScatterDataModifier::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<ScatterDataModifier *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->backgroundEnabledChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: _t->gridEnabledChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: _t->fontChanged((*reinterpret_cast< const QFont(*)>(_a[1]))); break;
        case 3: _t->shadowQualityChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->setSmoothDots((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->changeTheme((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: _t->changePresetCamera(); break;
        case 7: _t->changeLabelStyle(); break;
        case 8: _t->changeFont((*reinterpret_cast< const QFont(*)>(_a[1]))); break;
        case 9: _t->shadowQualityUpdatedByVisual((*reinterpret_cast< QtDataVisualization::QAbstract3DGraph::ShadowQuality(*)>(_a[1]))); break;
        case 10: _t->changeShadowQuality((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->setBackgroundEnabled((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: _t->setGridEnabled((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: _t->setVoxelParameters((*reinterpret_cast< const int(*)>(_a[1])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[2])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[3])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[4])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[5]))); break;
        case 14: _t->modifyScale((*reinterpret_cast< const double(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (ScatterDataModifier::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ScatterDataModifier::backgroundEnabledChanged)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (ScatterDataModifier::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ScatterDataModifier::gridEnabledChanged)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (ScatterDataModifier::*)(const QFont & );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ScatterDataModifier::fontChanged)) {
                *result = 2;
                return;
            }
        }
        {
            using _t = void (ScatterDataModifier::*)(int );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ScatterDataModifier::shadowQualityChanged)) {
                *result = 3;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject ScatterDataModifier::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_meta_stringdata_ScatterDataModifier.data,
    qt_meta_data_ScatterDataModifier,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *ScatterDataModifier::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ScatterDataModifier::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ScatterDataModifier.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int ScatterDataModifier::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 15)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 15;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 15)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 15;
    }
    return _id;
}

// SIGNAL 0
void ScatterDataModifier::backgroundEnabledChanged(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void ScatterDataModifier::gridEnabledChanged(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void ScatterDataModifier::fontChanged(const QFont & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void ScatterDataModifier::shadowQualityChanged(int _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
