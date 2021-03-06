#pragma once

#include "tree/TID.h"

#include "Reconstruct_traits.h"

#include <queue>
#include <list>
#include <memory>

namespace ant {


namespace reconstruct {

class UpdateableManager {
public:
    /**
     * @brief UpdateableManager initializes the manager
     * @param updateables list of updateable items to be managed
     */
    UpdateableManager(const std::list< std::shared_ptr<Updateable_traits> >& updateables_);

    /**
     * @brief UpdateParameters make the managed items ready for given currentPoint
     * @param currentPoint the time point
     */
    void UpdateParameters(const TID& currentPoint);

private:
    struct queue_item_t {
        TID NextChangePoint;
        Updateable_traits::Loader_t Item;
        queue_item_t(const TID& nextChangePoint,
                     Updateable_traits::Loader_t item) :
            NextChangePoint(nextChangePoint),
            Item(item)
        {}
        bool operator<(const queue_item_t& other) const {
            // invert ordering such that item with earliest change point
            // comes first in priority queue
            return other.NextChangePoint < NextChangePoint;
        }
    };

    std::priority_queue<queue_item_t> queue;

    std::list< std::shared_ptr<Updateable_traits> > updateables;
    TID lastFlagsSeen;

    void DoQueueLoad(const TID& currPoint,
                     Updateable_traits::Loader_t loader);
};


}} // namespace ant::reconstruct
