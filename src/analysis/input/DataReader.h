#pragma once

#include <memory>

namespace ant {

class Event;
class ReadTFiles;

namespace input {

/**
 * @brief Abstract base class for data input modules
 *
 * Data input modules read MC/Detector data from somewhere.
 * Examples:
 *  * goat file reader
 *  * new ant data foramt reader
 */
class DataReader {
protected:
    std::shared_ptr<ReadTFiles> files;
public:

    DataReader(const std::shared_ptr<ReadTFiles>& rootfiles);
    virtual ~DataReader();
    virtual std::shared_ptr<Event> ReadNextEvent() =0;

    virtual bool hasData() const =0;

    virtual long long EventsRead() const =0;

    /**
     * @brief Get total number of events available
     * @return number, -1 if unknown/does not apply
     */
    virtual long long TotalEvents() const { return -1; }

    class Exception : public std::runtime_error {
      using std::runtime_error::runtime_error; // use base class constructor
    };
};


}

}
