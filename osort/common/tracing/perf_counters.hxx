#ifndef F
#error Define F before including this file
#endif

// F(iscomputed, name, description, expression)
F(false, readCount, "Read pages", 0)
F(false, writeCount, "Write pages", 0)
F(false, accessCount, "Access count", 0)
F(false, swapCount, "Swap count", 0)
F(false, tlbHitCount, "TLB hit count", 0)
F(true, tlbHitRate, "TLB hit rate", static_cast<double>(tlbHitCount) / static_cast<double>(accessCount))

#undef F