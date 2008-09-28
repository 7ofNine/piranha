#! /bin/sh
# Detect cache size in kilobytes using /proc/cpuinfo.
echo `cat /proc/cpuinfo | grep "cache size" | head -n 1 | awk -F ':' '{print $2}' | awk '{print $1}'`
