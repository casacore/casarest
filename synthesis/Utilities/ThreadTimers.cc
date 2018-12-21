#include <synthesis/Utilities/ThreadTimers.h>
  
namespace casacore{
  DT operator- (const Timers & tLater, const Timers & tEarlier)
  {
    return DT (tLater.elapsed() - tEarlier.elapsed(),
               tLater.cpu() - tEarlier.cpu(), 0);
  }
};
