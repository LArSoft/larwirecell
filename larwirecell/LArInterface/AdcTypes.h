// AdcTypes.h

// David Adams
// December 2015
//
// Integer and floating type for Adc signals.

#ifndef AdcTypes_H
#define AdcTypes_H

#include <vector>

typedef short AdcCount;
typedef std::vector<AdcCount> AdcCountVector;

typedef float AdcSignal;
typedef std::vector<AdcSignal> AdcSignalVector;
typedef std::vector<AdcSignalVector> AdcSignalVectorVector;

typedef float AdcPedestal;

typedef std::vector<bool> AdcFilterVector;

typedef short AdcFlag;
typedef std::vector<AdcFlag> AdcFlagVector;

const AdcFlag AdcGood = 0;         // ADC sample is fine
const AdcFlag AdcUnderflow = 1;    // ADC sample is underflow
const AdcFlag AdcOverflow = 2;     // ADC sample is overflow
const AdcFlag AdcStuckOff = 3;     // ADC sample has bits stuck at 0
const AdcFlag AdcStuckOn = 4;      // ADC sample has bits stuck at 1
const AdcFlag AdcSetFixed = 5;     // ADC sample is set to a fixed value
const AdcFlag AdcInterpolated = 6; // ADC sample is interpolated
const AdcFlag AdcExtrapolated = 7; // ADC sample is extrapolated

typedef unsigned int AdcChannel;
typedef std::vector<AdcChannel> AdcChannelVector;

#endif
