#include "Setup.h"

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_EtaPrime :
        public Setup
{
public:
    virtual std::string GetName() const override {
        return "Setup_2014_EtaPrime";
    }

    Setup_2014_EtaPrime() {

        // setup the detectors of interest
        const auto trigger = std::make_shared<detector::Trigger>();
        const bool cherenkovInstalled = false;
        AddDetector(trigger);
        AddDetector<detector::EPT_2014>(GetBeamEnergy());
        AddDetector<detector::CB>();
        AddDetector<detector::PID_2014>();
        AddDetector<detector::TAPS_2013>(cherenkovInstalled, false); // no Cherenkov, don't use sensitive channels
        AddDetector<detector::TAPSVeto_2014>(cherenkovInstalled); // no Cherenkov


        // then calibrations need some rawvalues to "physical" values converters
        // they can be quite different (especially for the COMPASS TCS system), but most of them simply decode the bytes
        // to 16bit signed values
        /// \todo check if 16bit signed is correct for all those detectors
        const auto convert_MultiHit16bit = make_shared<calibration::converter::MultiHit16bit>();
        const auto convert_CATCH_Tagger = make_shared<calibration::converter::CATCH_TDC>(trigger->Reference_CATCH_TaggerCrate);
        const auto convert_CATCH_CB = make_shared<calibration::converter::CATCH_TDC>(trigger->Reference_CATCH_CBCrate);
        const auto convert_GeSiCa_SADC = make_shared<calibration::converter::GeSiCa_SADC>();

        // the order of the calibrations is important
        // add both CATCH converters first,
        // since they need to scan the detector read for their reference hit
        AddCalibration(convert_CATCH_Tagger);
        AddCalibration(convert_CATCH_CB);

        // then we add the others, and link it to the CATCH converters
        AddCalibration<calibration::Timing>(Detector_t::Type_t::EPT,
                                            convert_CATCH_Tagger,
                                            -325 // default offset in ns
                                            );
        AddCalibration<calibration::Timing>(Detector_t::Type_t::CB,
                                            convert_CATCH_CB,
                                            -325,      // default offset in ns
                                            interval<double>{-100, 100} // default time window cut in ns
                                            );
        AddCalibration<calibration::Timing>(Detector_t::Type_t::PID,
                                            convert_CATCH_CB,
                                            -325,
                                            interval<double>{-500, 500} // default time window cut in ns
                                            );
        AddCalibration<calibration::Timing>(Detector_t::Type_t::TAPS,
                                            convert_MultiHit16bit,
                                            -300, /// \todo different default for PbWO
                                            interval<double>{-500, 500},
                                            -0.100 /// \todo give measured time gains for BaF2
                                            );
        AddCalibration<calibration::Timing>(Detector_t::Type_t::TAPSVeto,
                                            convert_MultiHit16bit,
                                            160,
                                            interval<double>{-1000, 1000}, /// \todo make this window smaller...
                                            -0.05 // default gain
                                            );

        AddCalibration<calibration::Integral>(Detector_t::Type_t::CB,
                                              convert_GeSiCa_SADC,
                                              0,    // default pedestal in raw
                                              0.07, // default gain
                                              2     // default threshold in MeV
                                              );

        AddCalibration<calibration::Integral>(Detector_t::Type_t::PID,
                                              convert_MultiHit16bit,
                                              100,    // default pedestal in raw
                                              0.014,  // default gain
                                              0.001   // default threshold in MeV
                                              );

        AddCalibration<calibration::Integral>(Detector_t::Type_t::TAPS,
                                              convert_MultiHit16bit,
                                              100,   // default pedestal in raw
                                              0.30,  // default gain
                                              1      // default threshold in MeV
                                              );

        AddCalibration<calibration::Integral>(Detector_t::Type_t::TAPSVeto,
                                              convert_MultiHit16bit,
                                              100,     // default pedestal in raw
                                              0.010, // default gain
                                              0.1    // default threshold in MeV
                                              );
    }

    virtual double GetBeamEnergy() const override {
        return 1604.0;
    }

    virtual cluster_thresholds_t GetClusterThresholds() const override {
        return {
            {Detector_t::Type_t::CB,   15}, // in MeV
            {Detector_t::Type_t::TAPS, 20}, // in MeV
        };
    }

    bool Matches(const THeaderInfo& header) const override {
        if(!Setup::Matches(header))
            return false;
        /// \todo Make beamtime match stricter than just detectors
        return true;

    }

    void BuildMappings(std::vector<hit_mapping_t>& hit_mappings,
                       std::vector<scaler_mapping_t>& scaler_mappings) const
    {
        Setup::BuildMappings(hit_mappings, scaler_mappings);
        // you may tweak the mapping at this location here
        // for example, ignore elements
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_EtaPrime)

}}} // namespace ant::expconfig::setup
