#pragma once
#include "../config.h"

namespace genesis {
#if defined(_SINGLE) || defined(_MIXED)
using WorkingReal = float;
#else
using WorkingReal = double;
#endif
#if defined(_SINGLE)
using IntegrationReal = float;
#else
using IntegrationReal = double;
#endif
}
