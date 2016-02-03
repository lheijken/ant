#include "Logger.h"

#include "TError.h"
#include <sstream>

// setup the logger, will be compiled as a little library
INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace ant::logger;

void SetupLogger(int argc, char* argv[]) {
    START_EASYLOGGINGPP(argc, argv);
    SetupLogger();
}

void SetupLogger() {
    el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
    el::Configurations loggerConf;
    loggerConf.setToDefault();
    loggerConf.setGlobally(el::ConfigurationType::Format, "%datetime [%level] %fbase:%line : %msg");
    loggerConf.setGlobally(el::ConfigurationType::ToFile, "false");
    loggerConf.set(el::Level::Verbose,  el::ConfigurationType::Format, "%datetime [%level-%vlevel] %fbase:%line : %msg");
    el::Loggers::reconfigureLogger("default", loggerConf);

    // set ROOT error handler
    SetErrorHandler([] (
                    int level, Bool_t, const char *location,
                    const char *msg) {
        // Bool_t is abort, not used for now...

        if (level < gErrorIgnoreLevel)
            return;

        stringstream ss;
        ss << "ROOT<" << location << ">: " << msg;
        if(level < kInfo) { LOG(INFO) << ss.str(); }
        else if(level == kInfo) { LOG(INFO) << ss.str(); }
        else if(level == kWarning) { LOG(WARNING) << ss.str(); }
        else {
            LOG(ERROR) << ss.str();
            if(DebugInfo::nProcessedEvents>=0)
                LOG(DEBUG) << "nProcessEvents=" << DebugInfo::nProcessedEvents;
            if(DebugInfo::nUnpackedBuffers>=0)
                LOG(DEBUG) << "nUnpackedBuffers=" << DebugInfo::nUnpackedBuffers;
        }

    });
}

long long DebugInfo::nProcessedEvents = -1;
int DebugInfo::nUnpackedBuffers = -1;

