#include "KnobsTF1.h"

#include "TF1.h"

using namespace ant::calibration::gui;
using namespace ant::calibration::gui::KnobsTF1;


ParameterKnob::ParameterKnob(const std::string& Name, TF1* Func, int par, IndicatorProperties::Type_t type, Color_t color, double LineW):
    IndicatorKnob(Name, IndicatorProperties(type,color,LineW)),
    func(Func),
    parameter_index(par)
{
}

double ParameterKnob::get() const
{
    return func->GetParameter(parameter_index);
}

void ParameterKnob::set(double a)
{
    func->SetParameter(parameter_index,a);
}


FixedParameterKnob::FixedParameterKnob(const std::string& Name, TF1* Func, int par, IndicatorProperties::Type_t type, Color_t color, double LineW):
    IndicatorKnob(Name, IndicatorProperties(type,color,LineW)),
    func(Func),
    parameter_index(par)
{
}

double FixedParameterKnob::get() const
{
    return func->GetParameter(parameter_index);
}

void FixedParameterKnob::set(double a)
{
    (void)a;
}



RangedParameterKnob::RangedParameterKnob(const std::string& Name, TF1* Func, int par,
                                         RangedParameterKnob::ConstraintType constraint_type_,
                                         IndicatorProperties::Type_t gui_type, Color_t color, double LineW) :
ParameterKnob(Name, Func, par, gui_type, color, LineW),
  constraint_type(constraint_type_)
{
}

void RangedParameterKnob::set(double a)
{
    double min, max;
    func->GetRange(min, max);

    switch(constraint_type) {
    case ConstraintType::lowerThanMin:
        if(min<a)
            return;
        break;
    case ConstraintType::lowerThanMax:
        if(max<a)
            return;
        break;
    case ConstraintType::higherThanMax:
        if(max>a)
            return;
        break;
    case ConstraintType::higherThanMin:
        if(min>a)
            return;
        break;
    default:
        return;
    }
    ParameterKnob::set(a);
}

TransformedParameterKnob::TransformedParameterKnob(const std::string& Name, TF1* Func, int par,
                                                   transformation_t trafo,
                                                   transformation_t trafo_inverse,
                                                   IndicatorProperties::Type_t type, Color_t color, double LineW) :
    ParameterKnob(Name, Func, par, type, color, LineW),
    transformation(trafo),
    transformation_inverse(trafo_inverse)
{

}

double TransformedParameterKnob::get() const
{
    return transformation(func->GetParameter(parameter_index), func);
}

void TransformedParameterKnob::set(double a)
{
    func->SetParameter(parameter_index, transformation_inverse(a, func));
}


ReferenceParameterKnob::ReferenceParameterKnob(const std::string& Name, TF1* Func,
                                               int par, int reference,
                                               IndicatorProperties::Type_t type, Color_t color, double LineW):
    TransformedParameterKnob(Name,Func,par,
                             [reference] (double a, TF1* f) { return a + f->GetParameter(reference); },
                             [reference] (double a, TF1* f) { return a - f->GetParameter(reference); },
                             type,color,LineW)
{
}

RangeKnob::RangeKnob(const std::string& Name, TF1* Func, RangeEndType Type, Color_t color, double LineW)
    :IndicatorKnob(Name,{IndicatorProperties::Type_t::slider_vertical,color,LineW}),
      func(Func),
      type(Type)
{
}

double RangeKnob::RangeKnob::get() const
{
    double min, max;
    func->GetRange(min,max);

    if(type==RangeEndType::lower) {
        return min;
    }

    return max;
}

void RangeKnob::RangeKnob::set(double a)
{
    double min, max;
    func->GetRange(min,max);

    switch (type) {

    case RangeEndType::lower:
        if(a < max)
            func->SetRange(a,max);
        break;

    case RangeEndType::upper:
            if(a>min)
          func->SetRange(min,a);
        break;

    default:
        break;
    }
}

FixedRangeKnob::FixedRangeKnob(const std::string& Name, TF1* Func, RangeEndType Type, Color_t color, double LineW):
    RangeKnob(Name, Func, Type, color, LineW)
{
}

double FixedRangeKnob::FixedRangeKnob::get() const
{
    double min, max;
    func->GetRange(min,max);

    if(type==RangeEndType::lower) {
        return min;
    }

    return max;
}

void FixedRangeKnob::FixedRangeKnob::set(double a)
{
    (void)a;
}

ParLimitKnob::ParLimitKnob(const std::string& Name, TF1* Func, int par, LimitEndType LimType, IndicatorProperties::Type_t type, Color_t color, double LineW):
    IndicatorKnob(Name, IndicatorProperties(type,color,LineW)),
    func(Func),
    parameter_index(par),
    limtype(LimType)
{
}

double ParLimitKnob::get() const
{
    double low, high;
    func->GetParLimits(parameter_index, low, high);
    if(limtype==LimitEndType::lower) {
        return low;
    }
    return high;
}

void ParLimitKnob::set(double a)
{
    double parval = func->GetParameter(parameter_index);
    double high, low;
    func->GetParLimits(parameter_index, low, high);
    if(limtype==LimitEndType::lower) {
        if(high==0){
            func->SetParLimits(parameter_index,a,2*parval-a);
            return;
        }
        func->SetParLimits(parameter_index,a,high);
        return;
    }
    if(limtype==LimitEndType::upper) {
        if(low==0){
            func->SetParLimits(parameter_index,2*parval-a,a);
            return;
        }
        func->SetParLimits(parameter_index,low,a);
        return;
    }
}

