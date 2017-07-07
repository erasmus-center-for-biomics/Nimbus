// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <stdio.h>

// additional headers from the standard library
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

// headers from the threadingutils library
#include "Signal.h"
#include "TQueue.h"

// headers from the libnimbus library
#include "Read.h"
#include "SAMrecord.h"
#include "Amplicon.h"
#include "Alignment.h"
#include "io.h"
#include "AlignmentBuilder.h"
#include "AmpliconAlignment.h"

#if _POSIX_VERSION >= 200112L
#include <unistd.h>
#endif

#define LIMIT 5000
#define THREADSLEEPTIME 1 
