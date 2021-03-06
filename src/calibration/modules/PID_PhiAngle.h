#pragma once

#include "Calibration.h"
#include "base/interval.h"

class TGraph;

namespace ant {

namespace expconfig {
namespace detector {
struct PID;
}}

namespace calibration {

class DataManager;

namespace gui {
class FitGaus;
}


class PID_PhiAngle : public Calibration::Module
{
public:
    PID_PhiAngle(const std::shared_ptr<expconfig::detector::PID>&  pid,
            const std::shared_ptr<DataManager>& calmgr
            );
    virtual ~PID_PhiAngle();

    class TheGUI : public gui::CalibModule_traits {
    protected:
        std::shared_ptr<DataManager> calibrationManager;
        std::shared_ptr<expconfig::detector::PID> pid_detector;
        class _FitGauss;
        std::shared_ptr<_FitGauss> func;

        gui::CalCanvas* canvas;

        TH1*  h_projection = nullptr;
        TGraph* h_result;

        std::vector<double> angles;
        std::vector<double> previousAngles;
        std::map< unsigned, std::vector<double> > fitParameters;

        double phi_offset = std::numeric_limits<double>::quiet_NaN();

        bool IgnorePreviousFitParameters = false;
    public:
        TheGUI(const std::string& basename,
               const std::shared_ptr<DataManager>& calmgr,
               const std::shared_ptr<expconfig::detector::PID>& pid
               );
        virtual ~TheGUI();

        virtual std::shared_ptr<TH1> GetHistogram(const WrapTFile& file) const override;
        virtual unsigned GetNumberOfChannels() const override;
        virtual void InitGUI(gui::ManagerWindow_traits& window) override;

        virtual void StartSlice(const interval<TID>& range) override;
        virtual DoFitReturn_t DoFit(const TH1& hist, unsigned channel) override;
        virtual void DisplayFit() override;
        virtual void StoreFit(unsigned channel) override;
        virtual bool FinishSlice() override;
        virtual void StoreFinishSlice(const interval<TID>& range) override;
    }; // TheGUI

    virtual void GetGUIs(std::list<std::unique_ptr<calibration::gui::CalibModule_traits> >& guis, ant::OptionsPtr options) override;

    // Updateable_traits interface
    virtual std::list<Loader_t> GetLoaders() override;

protected:
    std::shared_ptr<expconfig::detector::PID> pid_detector;
    std::shared_ptr<DataManager> calibrationManager;



};

}}
