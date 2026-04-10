/****************************************************************************
** Meta object code from reading C++ file 'ScatterDataModifier.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../GUI/ScatterDataModifier.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ScatterDataModifier.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_ScatterDataModifier_t {
    QByteArrayData data[13];
    char stringdata0[178];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ScatterDataModifier_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ScatterDataModifier_t qt_meta_stringdata_ScatterDataModifier = {
    {
QT_MOC_LITERAL(0, 0, 19), // "ScatterDataModifier"
QT_MOC_LITERAL(1, 20, 18), // "setVoxelParameters"
QT_MOC_LITERAL(2, 39, 0), // ""
QT_MOC_LITERAL(3, 40, 11), // "voxel_size_"
QT_MOC_LITERAL(4, 52, 33), // "std::vector<std::vector<doubl..."
QT_MOC_LITERAL(5, 86, 8), // "X_axons_"
QT_MOC_LITERAL(6, 95, 8), // "Y_axons_"
QT_MOC_LITERAL(7, 104, 8), // "Z_axons_"
QT_MOC_LITERAL(8, 113, 8), // "R_axons_"
QT_MOC_LITERAL(9, 122, 13), // "X_astrocytes_"
QT_MOC_LITERAL(10, 136, 13), // "Y_astrocytes_"
QT_MOC_LITERAL(11, 150, 13), // "Z_astrocytes_"
QT_MOC_LITERAL(12, 164, 13) // "R_astrocytes_"

    },
    "ScatterDataModifier\0setVoxelParameters\0"
    "\0voxel_size_\0std::vector<std::vector<double> >\0"
    "X_axons_\0Y_axons_\0Z_axons_\0R_axons_\0"
    "X_astrocytes_\0Y_astrocytes_\0Z_astrocytes_\0"
    "R_astrocytes_"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ScatterDataModifier[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    9,   19,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int, 0x80000000 | 4, 0x80000000 | 4, 0x80000000 | 4, 0x80000000 | 4, 0x80000000 | 4, 0x80000000 | 4, 0x80000000 | 4, 0x80000000 | 4,    3,    5,    6,    7,    8,    9,   10,   11,   12,

       0        // eod
};

void ScatterDataModifier::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<ScatterDataModifier *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->setVoxelParameters((*reinterpret_cast< const int(*)>(_a[1])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[2])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[3])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[4])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[5])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[6])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[7])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[8])),(*reinterpret_cast< const std::vector<std::vector<double> >(*)>(_a[9]))); break;
        default: ;
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
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 1;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
