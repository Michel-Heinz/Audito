! Copyright (C) 2023 Michel Heinz

elemental type(_F_TYPE) function _CAT3(_F,_PowObjScalar_,_KIND)(obj, value) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value
#ifdef _FGH
    integer(i8) :: a
    real(r8) :: out(SIZE(obj%g),SIZE(obj%g))

    a = SIZE(obj%g)
    out = SPREAD(obj%g, 2, a) * SPREAD(obj%g, 1, a)
    new%h = value * ((value-1) * obj%f ** (value-2) * out + obj%f ** (value-1) * obj%h)
#endif
#ifdef _FGL
    new%l = value * ((value-1) * obj%f ** (value-2) * DOT_PRODUCT(obj%g, obj%g) + obj%f ** (value-1) * obj%l)
#endif

    new%f = obj%f ** value
    new%g = obj%f ** (value-1) * value * obj%g
end function _CAT3(_F,_PowObjScalar_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_PowScalarObj_,_KIND)(value, obj) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value
    real(r8) :: f_new, val
#ifdef _FGH
    integer(i8) :: a
    real(r8) :: out(SIZE(obj%g),SIZE(obj%g))

    a = SIZE(obj%g)
    out = SPREAD(obj%g, 2, a) * SPREAD(obj%g, 1, a)
#endif

    f_new = value ** obj%f
    new%f = f_new
#ifdef _INT
    val = REAL(value, r8)
    new%g = obj%g * f_new * log(val)
#ifdef _FGH
    new%h = (obj%h * f_new + f_new * out) * log(val)
#endif
#ifdef _FGL
    new%l = (obj%l * f_new + f_new * DOT_PRODUCT(obj%g, obj%g)) * log(val)
#endif
#else
    new%g = obj%g * f_new * log(value)
#ifdef _FGH
    new%h = (obj%h * f_new + f_new * out) * log(value)
#endif
#ifdef _FGL
    new%l = (obj%l * f_new + f_new * DOT_PRODUCT(obj%g, obj%g)) * log(value)
#endif
#endif
end function _CAT3(_F,_PowScalarObj_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_DivObjScalar_,_KIND)(obj, value) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value

    new%f = obj%f / value
    new%g = obj%g / value
#ifdef _FGH
    new%h = obj%h / value
#endif
#ifdef _FGL
    new%l = obj%l / value
#endif
end function _CAT3(_F,_DivObjScalar_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_DivScalarObj_,_KIND)(value, obj) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value
    type(_F_TYPE) :: temp

    temp = _CAT2(_F,_PowObjScalar_i1(obj, -1_i1))
    new%f = value * temp%f
    new%g = value * temp%g
#ifdef _FGH
    new%h = value * temp%h
#endif
#ifdef _FGL
    new%l = value * temp%l
#endif
end function _CAT3(_F,_DivScalarObj_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_AddObjScalar_,_KIND)(obj, value) result(new)
    _TYPE, intent(in) :: value
    class(_F_TYPE), intent(in) :: obj

    new%f = obj%f + value
    new%g = obj%g
#ifdef _FGH
    new%h = obj%h
#endif
#ifdef _FGL
    new%l = obj%l
#endif
end function _CAT3(_F,_AddObjScalar_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_AddScalarObj_,_KIND)(value, obj) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value

    new%f = value + obj%f
    new%g = obj%g
#ifdef _FGH
    new%h = obj%h
#endif
#ifdef _FGL
    new%l = obj%l
#endif
end function _CAT3(_F,_AddScalarObj_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_SubObjScalar_,_KIND)(obj, value) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value

    new%f = obj%f - value
    new%g = obj%g
#ifdef _FGH
    new%h = obj%h
#endif
#ifdef _FGL
    new%l = obj%l
#endif
end function _CAT3(_F,_SubObjScalar_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_SubScalarObj_,_KIND)(value, obj) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value

    new%f = value - obj%f
    new%g = - obj%g
#ifdef _FGH
    new%h = obj%h
#endif
#ifdef _FGL
    new%l = obj%l
#endif
end function _CAT3(_F,_SubScalarObj_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_MulObjScalar_,_KIND)(obj, value) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value

    new%f = obj%f * value
    new%g = obj%g * value
#ifdef _FGH
    new%h = obj%h * value
#endif
#ifdef _FGL
    new%l = obj%l * value
#endif
end function _CAT3(_F,_MulObjScalar_,_KIND)


elemental type(_F_TYPE) function _CAT3(_F,_MulScalarObj_,_KIND)(value, obj) result(new)
    class(_F_TYPE), intent(in) :: obj
    _TYPE, intent(in) :: value

    new%f = value * obj%f
    new%g = value * obj%g
#ifdef _FGH
    new%h = value * obj%h
#endif
#ifdef _FGL
    new%l = value * obj%l
#endif
end function _CAT3(_F,_MulScalarObj_,_KIND)

#undef _F_TYPE
#undef _KIND
#undef _TYPE
#undef _F
#undef _INT
