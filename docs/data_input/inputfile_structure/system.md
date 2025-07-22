# System

## Overview

This entry contains RHM global settings.

## Description

- `"LogLevel"` (string, [<span style="color:red">TOV, COOL, RH</span>]) **:** Defines verbosity of the programs. Choose from ["Error", "Info", "Debug", "Trace"], where "Error" only prints the halting reason, "Info" provides checkpoint information, "Debug" signalizes about potential faults during the main program loop and "Trace" delves into main loop procedures. Defaults to "Error".
- `"LogPath"` (string, [<span style="color:red">TOV, COOL, RH</span>]) **:** Path to the log file. Choose from ["stdout", "stderr"] or specify a file name(log will be appended). If not specified, logs to "stderr".